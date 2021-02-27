"""Ce fichier python contient les utilitaires du programme
"""
import sys, os
import subprocess

def readConfig(configFile):
	"""Lit le fichier de configuration
		-Args:
			_configfile: fichier de configuration
	"""
	process,goread,GenerateVectors = True,True,False
	with open(configFile, 'r') as conf:
		for line in conf:
			if line.startswith('#'):
				continue
			if line.startswith('GenerateVectors'):
				if line.strip().split('=')[1] == 'YES':
					GenerateVectors = True
			if line.startswith('PathToPdbFile'):
				PathToPdbFile = line.strip().split('=')[1]
				if not PathToPdbFile:
					print("ERR:Path to pdb structure file must be indicated")
					sys.exit()
			if line.startswith('FolderVectors'):
				FolderVectors = line.strip().split('=')[1]
				if not FolderVectors:
					print("ERR:Path to a folder used to store normal mode file must be indicated")
					sys.exit()
			if line.startswith('PathToAxis'):
				PathToAxis = line.strip().split('=')[1]
				if PathToAxis == 'NONE':
					print("Axis of protein missing, computed if volume constraints")
					PathToAxis = None
					# sys.exit()
			if line.startswith('EigenFile'):
				EigenFile = line.strip().split('=')[1]
				if EigenFile == 'NONE':
					process = False
			if line.startswith('ModesFile'):
				ModesFile = line.strip().split('=')[1]
				if ModesFile == 'NONE':
					goread = False
			if line.startswith('Type'):
				Type = line.strip().split('=')[1]
			if line.startswith('Contrainte'):
				Contrainte = line.strip().split('=')[1]
			if line.startswith('Collectivity'):
				Coll = line.strip().split('=')[1]
			if line.startswith('Selection'):
				Sel = line.strip().split('=')[1]
			if line.startswith('Temperature'):
				Tmp = line.strip().split('=')[1]
	return Config(PathToPdbFile,GenerateVectors,FolderVectors,PathToAxis,EigenFile,process,ModesFile,goread,Type,Contrainte,Coll,Sel,Tmp)

class Config:
	"""Garde en mémoire la configuration donnée par l'utilisateur
	"""
	def __init__(self,pdbFile,GenerV,folderV,pathToAxis,EigenFile,process,ModesFile,goread,Type,Contrainte,Coll,Sel,Tmp):
		self.pdb = pdbFile
		self.pdbName = self.pdb.split('.')[1].split('/')[2]
		self.generateNMA = GenerV
		self.eigenFolder = folderV
		self.axis = pathToAxis
		self.efile = EigenFile
		self.process = process
		self.mfile = ModesFile
		self.goread = goread
		self.type = Type
		self.contr = float(Contrainte)
		self.coll = float(Coll) #seuil de collectivité
		self.selection = Sel # selection de residu pour le calcul de SAS
		self.temp = int(Tmp)

	def check_axis(self):
		"""Create axis file from pdb if axis missing and needed
		"""
		if self.type == 'Volume' and self.axis == None:
			axe, H = determineSymetrie(self.pdb)
			self.axis = writeSymetryAxis(self.pdbName,axe,H)


############################################################################################################################
# writers/ readers

def ecrit_eigenVectors(vector, value, EigenFile, wmode, num):
    """Ecrit dans un fichier texte les vecteurs propres et valeurs propores calculés par 
    le module mdtraj.
        -Args:
            _vector: le vecteur propre
            _value: la valeur propre
            _EigenFile: chemin + nom du fichier 
            _wmode: mode d'écriture du fichier, append ou write
            _num: indice du mode
    """
    if wmode == 0:
        wmode = 'w'
    else:
        wmode = 'a'

    with open(EigenFile, wmode) as eig:
        eig.write('eigen vector {} : eigen value {}\n'.format(num,value))
        for elem in vector:
            eig.write('{}\n'.format(elem))


def write_PDB(pdbOri,trajectoire,p,numode=0,typ='ATOM'):
	"""Ecrit un fichier PDB. Modulaire
        -Args:
            _pdbOri: fichier pdb d'origine
            _trajectoire: la trajectoire selon un ou plusieur modes
            _p: mode d'ecriture 'w' ou 'a' 
            _numode:si mode normal seul alors entier correspondant a son indice. Sinon X pour les combinaisons
            _typ: type d'atome a prendre en compte, indique HETATM si molecule inorganique
	"""
	keep =[]
	nameF = 'out/Traj{}Mode{}'.format(numode,pdbOri.split('/')[-1])
	with open(pdbOri,'r') as fil:
		for line in fil:
			if line.startswith(typ):
				keep.append(line[0:27])
	with open(nameF,p) as out:
		if len(trajectoire) != len(keep)*3:
			for i,model in enumerate(trajectoire):
				out.write('MODEL        {}\n'.format(i))
				for atm in range(int(len(model)/3)):
					crd = 3*atm +2
					out.write('{}{:>11.3f} {:>7.3f} {:>7.3f}\n'.format(keep[atm],model[crd-2],model[crd-1],model[crd]))
				out.write('ENDMDL\n')
		else:
			out.write('MODEL        {}\n'.format(0))
			for atm in range(int(len(trajectoire)/3)):
				crd = 3*atm +2
				out.write('{}{:>11.3f} {:>7.3f} {:>7.3f}\n'.format(keep[atm],trajectoire[crd-2],trajectoire[crd-1],trajectoire[crd]))
			out.write('ENDMDL\n')

def writeSelectedModes(vector,value,name,C):
	"""Ecrit dans un fichier formaté les vecteurs propres correspondant aux modes 
	de plus haute collectivité
        -Args:
            _vector: vecteur propre
			_value: valeur propre
			_name: indice du vecteur
			_C: configuration
	"""
	print("ECRITURE.......")
	if name == 0:
		wmode = 'w'
	else:
		wmode = 'a'
	with open(C.mfile, wmode) as eig:
		eig.write('eigen vector:{:.5f} eigen value:{:.5f} frequence:{:.5f} amplitude:{:.5f} collectivite:{:.5f}\n'.format(name,float(value),float(vector[1]),float(vector[2]),float(vector[3])))
		for elem in vector[0]:
			eig.write('{}\n'.format(elem))

def readEigenVs(file):
	"""Lit un fichier comportant des vecteurs propres genéré par nma.py
	(fichier primaire de Konbi-Mod)
        -Args:
            _file: le fichier de vecteur propre au format txt
	"""
	Vectors = {}
	eigV = []
	with open(file,'r')as fil:
		for line in fil:
			if line.startswith('eigen'):
				if eigV:
					Vectors[eigv] = eigV
					eigV = []
				eigv = line.split()[6]
			else:
				eigV.append(float(line.strip()))	
	Vectors[eigv] = eigV
	return Vectors

def read_coord_PDB(nomfichier):
    """
    Lit le fichier pdb et ajoute les atomes de notre systeme dans une liste.
    Argument(s) : nom du fichier PDB
    Sortie : liste d'atomes
    """
    liste_atomes=[]
    with open(nomfichier,'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                numAtom = int(line[7:11])
                x = float(line[31:38])#*0.1#float(line[29:37])
                y = float(line[38:46])#*0.1#float(line[40:46])
                z = float(line[46:54])#*0.1# float(line[50:56])
                chain = line[21:22]
                resname = line[17:20]
                resNumber = line[23:27]
                xyz = [x,y,z]
                ty = ' '+str(line[13:14])#str(line[15:16])
                atome=Atom(numAtom,xyz,ty,chain,resname,resNumber)
                liste_atomes.append(atome)
    return(liste_atomes)

def readSelectedModes(eigenFile):
	"""Lit le fichier de vecteurs propres formaté produit par anime.py 
	"""
	with open(eigenFile, 'r') as eig:
		eigenV = {}
		for line in eig:
			if line.startswith('eigen'):
				val = float(line[line.index('value')+6:line.index('value')+13])
				frq = float(line[line.index('frequence')+10:line.index('frequence')+17])
				amp = float(line[line.index('amplitude')+10:line.index('amplitude')+17])
				col = float(line[line.index('collectivite')+13:line.index('collectivite')+20])
				eigenV[val] = [[],frq,amp,col]
			else:
				eigenV[val][0].append(float(line.strip()))
	return eigenV

############################################################################################################################
# système, commandes et creation de repertoires

def execute(cmd): 
    """Execute une commande donnée
    	-Args:
    		_cmd: la commande en chaine de char
    """
    try:
        subprocess.call(cmd, shell=True)
        print(cmd)

    except Exception as e:
        sys.stderr.write("Error: %s\n" % repr(e));

def verif_env(Config):
	"""Verifie que les repertoires sont présents et qu'ils sont nommées 
	en accord avec les noms utilisés dans le programme
		-Args:
			_Config: instance de la classe Config, contient des informations sur les noms de repertoires
	"""
	if not os.path.isdir(Config.eigenFolder):
		cmd = 'mkdir {}'.format(Config.eigenFolder)
		execute(cmd)

def outFolder(pathpname):
	"""Create folder name if not exist in current dir, else pass
	"""
	if not os.path.exists(pathpname):
		execute('mkdir {}'.format(pathpname))

############################################################################################################################
# utilitaires
def getxyz(Latm):
	"""Take List of atm and return list of coord
	"""
	xyz = []
	for a in Latm:
		xyz.append(a.traj[0])
		xyz.append(a.traj[1])
		xyz.append(a.traj[2])
	return xyz
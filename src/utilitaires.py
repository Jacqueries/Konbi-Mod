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

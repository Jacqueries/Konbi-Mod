"""Ce fichier python contient les utilitaires du programme
la classe config et la classe volume. 
"""
import sys, os
import subprocess
import numpy as np
import NMAContrainte as N
import ast

############################################################################################################################
# Configuration

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
		self.selection = self.reWriteSelection(Sel) # selection de residu pour le calcul de SAS ou de num atomes pour le calcul de distances
		self.temp = int(Tmp)

	def check_axis(self):
		"""Create axis file from pdb if axis missing and needed
		"""
		if self.type == 'Volume' and self.axis == None:
			axe, H = determineSymetrie(self.pdb)
			self.axis = writeSymetryAxis(self.pdbName,axe,H)

	def reWriteSelection(self,selection):
		"""Reformate la selection si le type d'analyse est Distance entre atomes
		"""
		if self.type == 'Distance':
			return ast.literal_eval(selection)
		return selection


############################################################################################################################
# Axe de symmetrie

def distance(xyzr,xyzd):
	"""Renvoie les distances au carré entre deux set de coordonnées
	"""
	return ((xyzr[0] - xyzd[0])**2+(xyzr[1] - xyzd[1])**2+(xyzr[2] - xyzd[2])**2)

def DistMaxPaires(coord):
	"""Détermine pour chaque point le point le plus éloigné et l'associe en une paire
	"""
	Paires = []
	for d in range(len(coord)):
		distMax = -1000
		for a in range(d+1,len(coord)):
			dist = distance(coord[d],coord[a])
			if dist > distMax:
				eph = [coord[d],coord[a],dist]
				distMax = dist
		Paires.append(eph)
	return Paires

def centroid(points):
	"""Estime le centroid d'un nuage de points 
	"""
	x,y,z=[],[],[]
	for i in range(len(points)):
		x.append(points[i][0])
		y.append(points[i][1])
		z.append(points[i][2])

	xmed = statistics.median(x)
	ymed = statistics.median(y)
	zmed = statistics.median(z)
	return (xmed,ymed,zmed)

def determineLargPlan(centre,nuage):
	"""prend un centre et un nuage de point en entrée et renvoie la moyenne des distances au centre
	"""
	sumDist = 0
	for point in nuage:
		sumDist += distance(point,centre)
	return np.sqrt(sumDist)/len(nuage)

def determineSymetrie(pdbfile):
	"""Approxime l'axe de symétrie d'un nuage de points 3D
		Prend en argument un fichier pdb
	"""
	liste_atm = NMC.read_coord_PDB(pdbfile)
	coords = [] # conteneur : coordonnées
	count_atm = 0
	for atm in liste_atm:
		coords.append(atm.xyz)
		count_atm += 1
	paires = DistMaxPaires(coords) # Associe par paires les atomes les plus distants
	paires_triee = sorted(paires, key = lambda x: x[2]) # Trie par distance des paires
	nPaires = len(paires_triee)
	two_perc = int(0.02 * nPaires) # On selectionne deux pourcent des paires parmis les plus distantes  
	keepD = [] # ensemble du centre 1
	keepA = [] # ensemble du centre 2
	for i in range(two_perc):
		keepD.append(paires_triee[-(1+i)][0])
		keepA.append(paires_triee[-(1+i)][1])
	cD = centroid(keepD) # calcule le centre 1
	cA = centroid(keepA) # calcule le centre 1
	# ces deux ensembles de points constituent l'axe
	lD = determineLargPlan(cD,keepD)
	lA = determineLargPlan(cA,keepA)
	hauteur = (lD+lA)/2
	return (cD,cA,hauteur)

def writeSymetryAxis(pdbName,axis,H):
	outFolder('./Struct/Axis')
	with open('./Struct/Axis/Axis_{}'.format(pdbName), 'w') as out:
		for point in axis:
			for coord in point:
				out.write('{}\n'.format(coord))
			out.write('#\n')
		out.write('H={}\n'.format(H))

############################################################################################################################
# Volume

class Volume:
	"""Classe Volume qui gère certaines opération du calcul de volume de canal
	Definit un pavé droit autour du canal le long d'un axe calculé et lu dans un fichier
	Le volume du canal est calculé en retirant au volume du pavé droit les volumes de 
	Van der Waals des atomes se trouvant à l'intérieur 
	"""
	def __init__(self,s1,s2,interval=2):
		self.s1 = s1
		self.s2 = s2
		self.volume = None
		self.voxels = None
		self.inter = interval

	def isInVol(self,pt):
		"""Check if pt is in Volume
		"""
		u = self.s1[0] - self.s2[3]
		v = self.s1[0] - self.s1[1]
		w = self.s1[0] - self.s1[3]
		if np.dot(u,self.s1[0]) > np.dot(u,pt) > np.dot(u,self.s2[3]) and \
			np.dot(v,self.s1[0]) > np.dot(v,pt) > np.dot(v,self.s1[1]) and \
			np.dot(w,self.s1[0]) > np.dot(w,pt) > np.dot(w,self.s1[3]):
			return True
		return False

	def calcVolume(self):
		"""calcule le volume du pavé droit
		"""
		u = self.s1[0] - self.s2[3]
		v = self.s1[0] - self.s1[1]
		w = self.s1[0] - self.s1[3]
		self.volume =  get_norme(u)*get_norme(v)*get_norme(w)
	
	def buildVoxels(self):
		"""Construit un nuage de points a l'interieur du volume
		"""
		v1 = np.array(vector(self.s1[0],self.s1[1]))
		v2 = np.array(vector(self.s1[0],self.s1[3]))
		v3 = np.array(vector(self.s1[0],self.s2[3]))
		it1 = int(get_norme(v1)/self.inter)
		it3 = int(get_norme(v3)/self.inter)
		v1,v2,v3 = v1/it1,v2/it1,v3/it3
		xyz = []
		for i in range(it3):
			pt = np.array(self.s1[0])+np.array(v3)*i
			for j in range(it1):
				pt = np.array(pt)+np.array(v2)
				for k in range(it1):
					pt = np.array(pt)+np.array(v1)
					xyz.append(pt)
				pt = pt - it1*np.array(v1)
		self.voxels = xyz

#####Utilitaires et calculs pour le volume#####

def move_on_axe(A,B,progress):
	"""Se deplace le long d'un vecteur AB
	"""
	x = A[0] + (B[0]-A[0])*progress
	y = A[1] + (B[1]-A[1])*progress
	z = A[2] + (B[2]-A[2])*progress
	return np.array([x,y,z])

def is_equal(n,m):
	"""Verifie l'égalité entre deux float
	"""
	return n+m < 1e-10 and n+m >= 0 or n+m > -1e-10 and n+m <= 0

def get_norme(AB):
	"""Renvoie la norme d'un vecteur AB
	"""
	return np.sqrt(AB[0]**2+AB[1]**2+AB[2]**2)

def vector(A,B):
	"""renvoie le vecteur AB
	"""
	return np.array([B[0]-A[0],B[1]-A[1],B[2]-A[2]])

def point_du_plan(P,A):
	"""Renvoie un point du plan
	"""
	return np.array([0,0,-P[3]/P[2]])

def is_on_the_plane(v,P):
    return is_equal(v[0]*P[0] + v[1]*P[1] + v[2]*P[2], P[3])

def plan(AB,B):
	"""Renvoie le plan colineaire a AB contenant B
	"""
	a = AB[0]
	b = AB[1]
	c = AB[2]
	d = -(a*B[0]+b*B[1]+c*B[2])
	return np.array([a,b,c,d])

def base_du_plan(A,B):
	"""Le plan P, le vecteur normal vNorm
	construit la base du plan v1,v2
	"""
	vNorm = vector(A,B)
	P = plan(vNorm,B)
	p1 = point_du_plan(P,B)
	v1 = vector(B,p1)
	v2 = np.cross(vNorm, v1)
	return (v1,v2)

def coins(v1,v2,Pt,h):
	"""Prend v1,v2 les bases d'un plan, un point Pt qui appartient à ce plan et h une longeur
	Renvoie les quatres points d'une surface (2*h)**2 centré sur Pt
	"""
	nv1 = get_norme(v1)
	nv2 = get_norme(v2)
	r1 = h/nv1
	r2 = h/nv2
	c0 = Pt + v1*r1 + v2*r2
	c1 = Pt - v1*r1 + v2*r2
	c2 = Pt - v1*r1 - v2*r2
	c3 = Pt + v1*r1 - v2*r2
	return np.array([c0,c1,c2,c3])

def surface(axe,h):
	"""Determine les deux surfaces faisant partie de deux plans orthogonaux à l'axe
	le point axe[0] appartient à s1 et axe[1] appartient à s2
	"""
	A = axe[0]
	B = axe[1]
	v1,v2 = base_du_plan(A,B)
	s1 = coins(v1,v2,B,h)
	v1,v2 = base_du_plan(B,A)
	s2 = coins(v1,v2,A,h)
	return (np.array([s1,s2]))

def verif(s1,s2,A,B):
	"""Vérifie que les surfaces s1 et s2 ont été correctement crées
	cad perpendiculairement à l'axe
	"""
	vNorm = vector(A,B)
	P1 = plan(vNorm,B)
	vNorm2 = vector(B,A)
	P2 = plan(vNorm2,A)
	for i in range(len(s1)):
		a = is_on_the_plane(s1[i],P2)
		b = is_on_the_plane(s2[i],P1)
		if not a or not b :
			print("Error is not on plane")

def get_axe(pathtoaxe):
	"""lecture du fichier axe dans ./Struct/Axes/
	retourne l'axe, et la hauteur du plan du volume
	"""
	ax = np.array([[0.,0.,0.],[0.,0.,0.]])
	with open(pathtoaxe,'r') as axe:
		p,i = 0,0
		for line in axe:
			if line.startswith('H='):
				H = float(line.strip().split('=')[1])

				continue
			if line.startswith('#'):
				p = 1
				i = 0
				continue
			ax[p][i] = line.strip()
			i+=1
	return (ax,H)

def build_example(file,axe,h): #build_example('Exemple.pdb',get_axe())
	"""Construit l'exemple du volume dans un fichier pdb
	Peut être superposé avec le fichier pdb pour visualiser l'encadrement du canal
	"""	
	s1,s2 = surface(axe,h)
	lon = get_norme(vector(axe[0],axe[1]))
	step = 1
	nstep = int(lon/step)+1
	verif(s1,s2,axe[1],axe[0])
	cx = axe[1]
	sx0 = s1[0]
	sx1 = s1[1]
	sx2 = s1[2]
	sx3 = s1[3]
	natm = 0
	with open (file,'w') as ex:
		for i in range(nstep):
			prg = (i+1)/nstep
			ex.write("ATOM{:>7d}{:>6s}{:>3s}{:>5d}{:>10.3f}{:>10.3f}{:>10.3f}{:>3d}{:>3d}\n".format(natm,'O ',"TIP3",1,cx[0],cx[1],cx[2],1,0 ))
			natm+=1
			ex.write("ATOM{:>7d}{:>6s}{:>3s}{:>5d}{:>10.3f}{:>10.3f}{:>10.3f}{:>3d}{:>3d}\n".format(natm,'O ',"TIP3",1,sx0[0],sx0[1],sx0[2],1,0 ))
			natm+=1
			ex.write("ATOM{:>7d}{:>6s}{:>3s}{:>5d}{:>10.3f}{:>10.3f}{:>10.3f}{:>3d}{:>3d}\n".format(natm,'O ',"TIP3",1,sx1[0],sx1[1],sx1[2],1,0 ))
			natm+=1
			ex.write("ATOM{:>7d}{:>6s}{:>3s}{:>5d}{:>10.3f}{:>10.3f}{:>10.3f}{:>3d}{:>3d}\n".format(natm,'O ',"TIP3",1,sx2[0],sx2[1],sx2[2],1,0 ))
			natm+=1
			ex.write("ATOM{:>7d}{:>6s}{:>3s}{:>5d}{:>10.3f}{:>10.3f}{:>10.3f}{:>3d}{:>3d}\n".format(natm,'O ',"TIP3",1,sx3[0],sx3[1],sx3[2],1,0 ))
			natm+=1
			cx = move_on_axe(axe[1],axe[0],prg)
			sx0 = move_on_axe(s1[0],s2[3],prg)
			sx1 = move_on_axe(s1[1],s2[2],prg)
			sx2 = move_on_axe(s1[2],s2[1],prg)
			sx3 = move_on_axe(s1[3],s2[0],prg)
		natm+=1

def BuildOBB(axe,h):
	"""Construit le volume R qui contient le canal
	Retourne R
	"""
	s1,s2 = surface(axe,h)
	verif(s1,s2,axe[1],axe[0])
	R = Volume(s1,s2)
	R.calcVolume()
	# R.buildVoxels() Methode d'estimation du volume trop lente pour l'instant
	return R


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

def writeWeights(wei):
	with open('out/Weights.txt','w') as w:
		for k,elem in enumerate(wei.combinaisonMax):
			w.write('Weights : {}\n'.format(k))
			for e in elem:
				w.write('{}\n'.format(e))

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
                atome=N.Atom(numAtom,xyz,ty,chain,resname,resNumber)
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
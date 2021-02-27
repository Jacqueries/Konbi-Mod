""" Programme principal, optimisation des poids associées aux modes normaux sous contrainte
Est lancé depuis Konbi-Mod
L'optimisation des poids se fait par un processus de Monte Carlo.
A la fin des itérations, le programme crée une trajectoire qui satisfait au mieux la contrainte
selon une combinaison linéaire des modes normaux de basse fréquence de la protéine
"""
import sys, os
from contextlib import contextmanager
import numpy as np
import copy
import anime as A
import weights as W
import utilitaires as U
import freesasa
import scipy.spatial.distance

############################################################################################################################
# classe Atome
class Atom:
    def __init__(self,numAtom,xyz,ty,chain,resname,resNumber):
        self.num = numAtom
        # self.numRes = numRes
        self.xyz = np.array(xyz) #origine
        self.traj = np.array(xyz) #déplacement
        self.ty = ty #type d'atome (C,N,O,H,S)
        self.chain = chain
        self.resname = resname
        self.resN = resNumber
        self.volume = None
        self.radii = self.affectradii(ty)#1.42 #angstrom
        self.neighbors = []
        self.Hnei = [] #neighbors for hessian construct
        self.Hdis = [] #distance associated with neighbors Hnei
        self.mass = self.affectmass(ty)

    def computeVolume(self):
    	self.volume = (4/3)*np.pi*self.radii**3

    def ajVoisin(self,numvoisin):
    	if not numvoisin in self.neighbors:
    		self.neighbors.append(numvoisin)

    def move_2(self,dpc):
    	self.traj = np.array(dpc)+self.xyz

    def move_W(self,dpc,t):
    	if t == 0:
    		self.traj = self.xyz
    	self.traj = np.array(dpc)+self.traj

    def AddHneiHdis(self,nei,dist):
    	"""add nei to Hnei and the distance 
    	"""
    	self.Hnei.append(nei)
    	self.Hdis.append(dist)

    def affectradii(self,ty):
    	dRadii = {'H':1.2,'C':1.7,'N':1.55,'O':1.52,
    	'S':1.8}
    	if ty not in dRadii.keys():
    		return 1.42
    	else:
    		return dRadii[ty]
    def affectmass(self,ty):
    	dmass = {'H':1.,'C':12.,'N':14.,'O':16.,
    	'S':32.}
    	if ty not in dmass.keys():
    		return 10.
    	else:
    		return dmass[ty]

############################################################################################################################
#utilitaires

def trieEigen(Eivector):
	"""Trie les vecteurs par collectivity 
	renvoie le dictionnaire avec comme clef le rang de collectivité
	permet de donner de l'importance (poids) aux collectivités fortes par construction
	"""
	sortedKeys = list({k: v for k, v in sorted(Eivector.items(), key=lambda item: item[1][1],reverse=True)}.keys())
	for i,k in enumerate(sortedKeys):
		Eivector[i] = Eivector[k]
		del Eivector[k]
	return Eivector

def visual(Latm,w,eigenvectors,pdb,C):
	"""Ecrit dans un fichier pdb chaque frame appartenant aux meilleures combinaisons de poids
	trouvé lors du processus de Monte carlo. Les poids correspondants a chacune de ces 
	frames seront aussi sauvegardés dans un fichier texte
		-Args:
			_Latm: liste d'atomes
			_w:instance de la classe weights
			_eigenvectors: vecteur propre
			_pdb: fichier pdb d'origine
			_C: configuration
	"""
	U.outFolder('out')
	for i,elem in enumerate(w.combinaisonMax):
		if i == 0:
			wmode = 'w'
		else:
			wmode = 'a'
		Latm = ReturnDisplacedLatm(Latm,eigenvectors,len(Latm),w.timeCollect[i],elem,C)
		U.write_PDB(pdb,U.getxyz(Latm),wmode,numode='visual',typ='ATOM')
	U.writeWeights(w)

def trajMaxCombinaison(Latm,eigenvectors,w,C):
	"""Ecrit dans un fichier la trajectoire sur 100 frames de la combinaison de modes qui 
	satisfait au mieux la contrainte posée
		-Args:
			_Latm: liste d'atomes
			_eigenvectors: vecteur propre
			_w:instance de la classe weights
			_C: configuration
	"""
	for t in range(100):
		if t == 0:
			wmode = 'w'
		else:
			wmode = 'a'
		Latm = ReturnDisplacedLatmW(Latm,eigenvectors,len(Latm),t,w.combinaisonMax[-1],C)
		U.write_PDB(C.pdb,U.getxyz(Latm),wmode,'X',typ='ATOM')


############################################################################################################################
# trajectoire

def ReturnDisplacedLatmW(Latm,modeLf,natm,t,weights,C):
	"""Modifie les coordonnées de la trajectoire pour l'écriture dans un fichier PDB
	Retourne une liste d'atome
		-Args:
			_Latm: liste d'atomes
			_modeLf: vecteur propre
			_natm: nombre d'atomes
			_t: pas de temps
			_weights: poids associés
			_C: configuration
	"""
	dpc = weightedDisplact(modeLf,t,natm,weights,C,Latm)
	for j in range(len(Latm)):
		i = int(j * 3)
		Latm[j].move_W([dpc[i],dpc[i+1],dpc[i+2]],t)
	return Latm

def ReturnDisplacedLatm(Latm,modeLf,natm,t,weights,C):
	"""Modifie les coordonnées selon la trajectoire
	Retourne une liste d'atome
		-Args:
			_Latm: liste d'atomes
			_modeLf: vecteur propre
			_natm: nombre d'atomes
			_t: pas de temps, ici vaudra toujours 0
			_weights: poids associés
			_C: configuration
	"""
	dpc = weightedDisplact(modeLf,t,natm,weights,C,Latm)
	for j in range(len(Latm)):
		i = int(j * 3)
		Latm[j].move_2([dpc[i],dpc[i+1],dpc[i+2]])
	return Latm

def weightedDisplact(Vect,t,natm,weights,C,Latm):
	"""Déplace les atomes selon une combinaison des modes
		-Args:
			_Vect: vecteur propre
			_t: pas de temps, ici vaudra toujours 0
			_natm: nombre d'atomes
			_weights: poids associés
			_C: configuration
			_Latm: liste d'atomes
	Référence : Dynamical Properties of the MscL of Escherichia coli: A Normal Mode Analysis
	H.Valadie,J.Lacapc,Y-H.Sanejouand and C.Etchebest
	"""
	dpc = []
	Tmp = C.temp # température indiquée dans le fichier de config
	for coor in range(natm*3):
		m = int(coor/3)
		dpc.append(0)
		for  k,mode in enumerate(Vect.keys()):
			amp = A.calc_amplitude(Vect[mode][1],Tmp) # calcul l'amplitude du mode
			dpc[coor] += weights[k]*Vect[mode][0][coor]*amp*np.cos(2*np.pi*Vect[mode][1]*t+np.pi)
		dpc[coor] = (1/(np.sqrt(Latm[m].mass)))*dpc[coor]
	return dpc

############################################################################################################################
# starters organisateurs
def fromLaunch(C):
	"""Lance le programme depuis Konbi-Mod.py
	nécessite la configuration C
	"""
	eigFile,pdb,axe = C.mfile, C.pdb, C.axis
	Latm = U.read_coord_PDB(pdb)
	eigenvectors = trieEigen(U.readSelectedModes(eigFile))
	Algo = C.type

	if Algo == 'Surface':
		selection = C.selection
		contrainte = C.contr
		Latm,w = goSurface(selection,Latm,contrainte,eigenvectors,len(Latm),C)
	if Algo == 'Volume':
		contrainte = C.contr
		Latm,w = goVolume(eigenvectors,Latm,contrainte,C)
	if w.combinaisonMax:
		visual(Latm,w,eigenvectors,pdb,C)
		trajMaxCombinaison(Latm,eigenvectors,w,C)
		writeMemory(w.memory)
	else:
		print('insufficient constraint, please review constraint before starting again')
		print('if you are using Type=Volume you can increase the H value in the Axis file')

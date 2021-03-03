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

import time
import multiprocessing

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
        self.mass = self.affectmass(ty)

    def computeVolume(self):
    	self.volume = (4/3)*np.pi*self.radii**3

    def ajVoisin(self,numvoisin):
    	if not numvoisin in self.neighbors:
    		self.neighbors.append(numvoisin)

    def move_2(self,dpc):
    	"""deplace les atomes par rapport à leur emplacement d'origine
    	"""
    	self.traj = np.array(dpc)+self.xyz

    def move_W(self,dpc,t):
    	"""deplace les atomes par rappport à leur emplacement précédent
    	"""
    	if t == 0:
    		self.traj = self.xyz
    	self.traj = np.array(dpc)+self.traj


    def affectradii(self,ty):
    	"""Affecte à un atome son rayon de van der waals
    	Source: PyMol
    	"""
    	dRadii = {'H':1.2,'C':1.7,'N':1.55,'O':1.52,
    	'S':1.8}
    	if ty not in dRadii.keys():
    		return 1.42
    	else:
    		return dRadii[ty]

    def affectmass(self,ty):
    	"""Affecte à un atome sa masse molaire
    	"""
    	dmass = {'H':1.,'C':12.,'N':14.,'O':16.,
    	'S':32.}
    	if ty not in dmass.keys():
    		return 10.
    	else:
    		return dmass[ty]

def dist_atomes(xyz1,xyz2):
    """
        Calcul la distance entre deux atomes
        Arguments:
            + atom1 : objet Atom
            + atom2 : objet Atom
        Retourne la distance entre ces atomes
    """
    return(np.sqrt((xyz1[0] - xyz2[0])**2 + (xyz1[1] - xyz2[1])**2 + (xyz1[2] - xyz2[2])**2))

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

def voisins(Latm,inside):
	"""Ajoute a un atome sa liste de voisins 
	Deux atomes sont considérés voisins si leurs volumes de van der valls se chevauchent
	"""
	for e in range(len(inside)):
		for k in range(e+1,len(inside)):
			i = inside[e]
			j = inside[k]
			dij = dist_atomes(Latm[i].traj,Latm[j].traj)

			if dij < Latm[i].radii + Latm[j].radii:
				Latm[i].ajVoisin(j)
				Latm[j].ajVoisin(i)
	return Latm

def retrieveData(iteration,metrique,C,mode):
	"""Ecrit les données de volume/ sasa/ ... au cours de l'optimisation
	"""
	file = '{}_{}.txt'.format(C.pdbName,C.temp)
	if mode:
		mode = 'a'
	else:
		mode = 'w'
	with open(file, mode) as fil:
		fil.write('{},{}\n'.format(iteration,metrique))

def writeMemory(memory):
	"""ecrit dans un fichier l'évolution des poids pour chaque mode a chaque itération
	"""
	with open('./out/memoryWeights.csv', 'w') as out:
		for i,mode in enumerate(memory):
			out.write('Mode{}'.format(i))
			for weight in mode:
				out.write(',{}'.format(weight))
			out.write('\n')

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
			dpc[coor] += weights[k]*Vect[mode][0][coor]*amp*np.cos(2*np.pi*Vect[mode][1]*t + np.pi )#+np.pi
		dpc[coor] = (1/(np.sqrt(Latm[m].mass)))*dpc[coor]
	return dpc

###UPDATE

def ReturnDisplacedLatmD(Latm,modeLf,natm,t,weights,C):
	"""Modifie les coordonnées selon la trajectoire
	Retourne une liste d'atome
		-Args:
			_Latm: liste d'atomes
			_modeLf: vecteur propre
			_natm: nombre d'atomes
			_t: pas de temps, ignoré
			_weights: poids associés
			_C: configuration
	"""
	for t in range(t+1):
		dpc = weightedDisplact(modeLf,t,natm,weights,C,Latm)
		for j in range(len(Latm)):
			i = int(j * 3)
			Latm[j].move_W([dpc[i],dpc[i+1],dpc[i+2]],t)
	return Latm

def weightedDisplactT(Vect,t,natm,weights,C,Latm):
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
	t1 = time.time()
	
	manager = multiprocessing.Manager()
	coor_dict = manager.dict()

	Tmp = C.temp # température indiquée dans le fichier de config

	nproc = multiprocessing.cpu_count()
	nbatch = int((natm*3)/(nproc*4))
	if natm%nproc != 0:
		nbatch += 1
	coor = 0
	for batch in range(nbatch):
		jobs = []
		count = 0
		while coor < (natm*3) and count < (nproc*4):
			# m = int(coor/3)
			p = multiprocessing.Process(target=parallelDisplacement, args=(coor_dict,coor,Vect,t,weights,Latm,Tmp))
			jobs.append(p)
			p.start()
			count += 1
			coor += 1

		for proc in jobs:
			proc.join()

	dpc = []
	for coor in range(natm*3):
		dpc.append(coor_dict[coor])
	t2 = time.time()
	print("Seconds = {:.2f}".format((t2-t1)))
	return dpc

def parallelDisplacement(coor_dict,coor,Vect,t,weights,Latm,Tmp):
	coor_dict[coor] = 0
	for  k,mode in enumerate(Vect.keys()):
		amp = A.calc_amplitude(Vect[mode][1],Tmp) # calcul l'amplitude du mode
		coor_dict[coor] += weights[k]*Vect[mode][0][coor]*amp*np.cos(2*np.pi*Vect[mode][1]*t)#+np.pi
	coor_dict[coor] = (1/(np.sqrt(Latm[int(coor/3)].mass)))*coor_dict[coor]

############################################################################################################################
# Monte Carlo Volume

def double_overlap(pos1, pos2, r1, r2):
    """
    double_overlap(pos1, pos2, r1, r2)
    Calculate the overlap volume of two spheres of radius r1, r2, at positions
        pos1, pos2
    """
    d = sum((pos1 - pos2) ** 2) ** 0.5
    # check they overlap
    if d >= (r1 + r2):
        return 0
    # check if one entirely holds the other
    if r1 > (d + r2):  # 2 is entirely contained in one
        return 4. / 3. * np.pi * r2 ** 3
    if r2 > (d + r1):  # 1 is entirely contained in one
        return 4. / 3. * np.pi * r1 ** 3

    vol = (np.pi * (r1 + r2 - d) ** 2 * (d ** 2 + (2 * d * r1 - 3 * r1 ** 2 +2 * d * r2 - 3 * r2 ** 2)+ 6 * r1 * r2)) / (12 * d)
    return vol

def volOverlap(Latm):
	"""Estime le volume commun de VdW
	"""
	sumOverlap = 0
	watchList = []
	for atm in Latm:
		for v in atm.neighbors:
			if v in watchList:
				continue
			watchList.append(v)
			sumOverlap += double_overlap(atm.traj,Latm[v].traj,atm.radii,Latm[v].radii)
	return sumOverlap

def CalcVolCanal(Latm,R):
	"""Estime le volume occupé par les atomes pour une configuration xyz dans le volume R
		retourne le volume du canal
		(estimation qui peut être affinée mais difficile mathématiquement)
	"""
	inside = []
	nstep = 1000
	volVdW = 0
	for i in range(len(Latm)):
		if R.isInVol(Latm[i].traj):
			inside.append(i)
			Latm[i].computeVolume()
			volVdW += Latm[i].volume
	Latm = voisins(Latm,inside)
	sumOverlap = volOverlap(Latm)
	cuttedVol = 0

	vpoche = R.volume- (volVdW- (sumOverlap+ cuttedVol))
	return vpoche

def goVolume(eigenvectors,Latm,contrainte,C):
	"""Lance l'optimisation des poids pour maximiser le volume
		-Args:
			_eigenvectors: les vecteurs propres selectionés
			_Latm: la liste d'atomes
			_contrainte: le ratio d'augmentation ou de diminution de la sasa pour la selection
			_C: configuration
	"""
	natm = len(Latm)
	axe,h = U.get_axe(C.axis)
	U.build_example('./Struct/OBB.pdb',axe,h) # construit une représentation cubique du volume 
	R = U.BuildOBB(axe,h) # construit l'objet Volume
	volCanali = CalcVolCanal(Latm,R) # calule le volume initial du canal
	contrainte = contrainte*volCanali # définition de la contrainte
	modeE = 0
	aj = 1 # permet d'aller vers une augmentation ou une diminution du volume
	if volCanali > contrainte:
		aj = - 1 # si contrainte est inférieure au volume initial alors on va vers une diminution
	volCanalmax,volDefault = volCanali,volCanali
	print("Volume initial : {}".format(volCanali))
	w = W.Weights(len(eigenvectors)) # initialisation des poids (aléatoire)
	d,prgs = 0,0 # d : pas de temps utilisé, prgs : décompte des itérations
	Svol = []
	while  contrainte*aj > volCanali*aj and prgs < C.niter: # tant que la contrainte n'est pas satisfaite et que le nombre d'iteration nest pas atteint
		w.reajustWeights(prgs) # reajustement des poids
		Latm = ReturnDisplacedLatm(Latm,eigenvectors,natm,d,w.weights,C) # deplacement des atomes selon les modes pondérés
		volCanali = CalcVolCanal(Latm,R) # calcul de la surface accessible au solvent après déplacement
		print(volCanali)
		if volCanali*aj > volCanalmax*aj: # si augmentation du volume et post saturation
			w.reajustLimits(prgs) # reajustement des bornes des poids (entre -1 et 1)
			w.saveCombinaisonMax(d) # sauvegarder cette combinaison 
			volCanalmax = volCanali
			Svol = []
			print("Etape {}, évolution du volume du canal : {:.3f}".format(prgs,aj*(volCanali-volDefault)))
			# retrieveData(prgs,(volCanali-volDefault),C,modeE)
			# modeE+=1
		else :
			Svol,w = W.watch(volCanali,Svol,w)
			w.precState() # retourne au vecteur de poids précédent si pas d'amélioration
		prgs+=1
	return [Latm,w]


############################################################################################################################
# Monte Carlo Surface accessible au solvant

def calcSASA(Latm,selection):
	"""Calcule la surface accessible au solvent (SAS) des acides aminés de la selecion
	Retourne la SAS pour une sélection donnée
	"""
	freesasa.setVerbosity(1)
	structure = freesasa.Structure()
	for a in Latm:
		structure.addAtom(a.ty,a.resname,a.resN,a.chain,a.traj[0],a.traj[1],a.traj[2])
	result = freesasa.calc(structure)
	selections = freesasa.selectArea((selection,'all, resn ala'), structure, result)
	return selections[selection.split()[0][:-1]]

def goSurface(selection,Latm,contrainte,eigenvectors,natm,C):
	"""Lance l'optimisation des poids pour maximiser la suface accessible au solvent
	de la selection selon la contrainte donnée par l'utilisateur
		-Args:
			_selection: la selection d'acide aminés pour le calcul de surface accessible 
			au solvant, abrégé sasa ou sas
			_Latm: la liste d'atomes
			_contrainte: le ratio d'augmentation ou de diminution de la sasa pour la selection
			_eigenvectors: les vecteurs propres selectionés
			_natm:le nombre d'atomes
			_C: configuration

	"""
	saInit = calcSASA(Latm,selection) # calcule la sas de la selection
	contrainte = contrainte*saInit # définition de la contrainte
	modeE = 0
	aj = 1 # permet d'aller vers une augmententation ou diminution la sas
	if saInit > contrainte:
		aj = - 1
	SASAmax, saDefault = saInit, saInit
	print("Surface initiale de {} : {}".format(selection,saInit))
	w = W.Weights(len(eigenvectors)) # initialisation des poids (aléatoire)
	prgs,d = 0,0 # d : pas de temps utilisé, prgs : décompte des itérations
	SASA = []
	while contrainte*aj > saInit*aj and prgs < C.niter: # tant que la contrainte n'est pas satisfaite et que le nombre d'iteration nest pas atteint
		w.reajustWeights(prgs) # reajustement des poids
		Latm = ReturnDisplacedLatm(Latm,eigenvectors,natm,d,w.weights,C) # deplacement des atomes selon les modes pondérés
		saInit = calcSASA(Latm,selection) # calcul de la surface accessible au solvent après déplacement

		if saInit*aj > SASAmax*aj : # si sas init est depasse
			w.reajustLimits(prgs) # reajustement des bornes des poids (entre -1 et 1)
			w.saveCombinaisonMax(d) # sauvegarder cette combinaison 
			SASAmax = saInit
			SASA = []
			print("Etape {}, évolution de la surface accessible au solvant : {:.3f}".format(prgs,aj*(saInit- saDefault)))
			# retrieveData(prgs,(saInit- saDefault),C,modeE)
			modeE+=1
		else :
			SASA,w = W.watch(saInit,SASA,w)		
			w.precState() # retourne au vecteur de poids précédent si pas d'amélioration	
		prgs+=1

	return [Latm,w]

############################################################################################################################
# Distance

def calcDistance(paires,Latm):
	"""Prend en entrée des paires d'atomes (liste de listes) contenant des numéros d'atomes
	calcule la distance pour chaque paire et renvoie une liste correspondante
	Les numéros d'atomes doivent correspondre à la numérotation du fichier pdb
		-Args:
			_paires: liste de liste de paires de numéros d'atomes
			_Latm: liste d'objet Atome
	"""
	liste_distance = []
	for paire in paires:
		liste_distance.append(dist_atomes(Latm[paire[0]-1].traj,Latm[paire[1]-1].traj))
	return liste_distance

def goDistance2(selection,Latm,contrainte,eigenvectors,natm,C):
	"""Lance l'optimisation des poids pour maximiser la distance
	de la selection selon la contrainte donnée par l'utilisateur
		-Args:
			_selection: la selection de numéro d'atomes pour lesquel on veut contraindre la distance
			sous forme de liste de listes de paires
			_Latm: la liste d'atomes
			_contrainte: le ratio d'augmentation ou de diminution de la distance moyenne de la selection
			_eigenvectors: les vecteurs propres selectionés
			_natm:le nombre d'atomes
			_C: configuration
	"""
	distInit = np.mean(calcDistance(selection,Latm)) # calcule la sas de la selection
	contrainte = contrainte*distInit # définition de la contrainte
	# modeE = 0
	aj = 1 # permet d'aller vers une augmententation ou diminution la distance
	if distInit > contrainte:
		aj = - 1
	distmax, distDefault = distInit, distInit
	print("Distance initiale des atomes de la selection {} : {:.3f}".format(selection,distInit))
	w = W.Weights(len(eigenvectors)) # initialisation des poids (aléatoire)
	prgs,d = 0,0 # d : pas de temps utilisé, prgs : décompte des itérations
	Distances = []
	while contrainte*aj > distInit*aj and prgs < C.niter: # tant que la contrainte n'est pas satisfaite et que le nombre d'iteration nest pas atteint
		w.reajustWeights(prgs) # reajustement des poids
		t1 = time.time()
		Latm = ReturnDisplacedLatm(Latm,eigenvectors,natm,d,w.weights,C) # deplacement des atomes selon les modes pondérés
		t2 = time.time()
		distInit = np.mean(calcDistance(selection,Latm)) # calcul de la surface accessible au solvent après déplacement
		t3 = time.time()
		if distInit*aj > distmax*aj: # si sas init est depasse et que iterations post saturation
			w.reajustLimits(prgs) # reajustement des bornes des poids (entre -1 et 1)
			w.saveCombinaisonMax(d) # sauvegarder cette combinaison 
			distmax = distInit
			Distances = []
			print("Etape {}, évolution de la distance : {:.3f}".format(prgs,aj*(distInit- distDefault)))
			# retrieveData(prgs,(distInit- distDefault),C,modeE)
			# modeE+=1

		else :
			Distances,w = W.watch(distInit,Distances,w)		
			w.precState() # retourne au vecteur de poids précédent si pas d'amélioration	
		# print("Seconds = {:.2f} , {:.2f}, {:.2f}".format((t2-t1),(t3-t2),(((t3-t2)-(t2-t1)))))
		prgs+=1
	return [Latm,w]

###Update

def goDistance(selection,Latm,RatioC,eigenvectors,natm,C):
	"""Lance l'optimisation des poids pour maximiser la distance
	de la selection selon la contrainte donnée par l'utilisateur
		-Args:
			_selection: la selection de numéro d'atomes pour lesquel on veut contraindre la distance
			sous forme de liste de listes de paires
			_Latm: la liste d'atomes
			_RatioC: le ratio d'augmentation ou de diminution de la distance moyenne de la selection
			_eigenvectors: les vecteurs propres selectionés
			_natm:le nombre d'atomes
			_C: configuration
	"""
	distInit = calcDistance(selection,Latm) # calcule la sas de la selection
	contrainte = []
	for i in range(len(distInit)): 
		contrainte.append(RatioC*distInit[i]) # définition de la contrainte
	# modeE = 0
	aj = 1 # permet d'aller vers une augmententation ou diminution la distance
	if distInit[0] > contrainte[0]:
		aj = - 1
	Satisfait = False
	distmax, distDefault = copy.deepcopy(distInit), copy.deepcopy(distInit)
	print("Distance initiale des atomes de la selection {} : {}".format(selection,distInit))
	w = W.Weights(len(eigenvectors)) # initialisation des poids (aléatoire)
	prgs,d = 0,0 # d : pas de temps utilisé, prgs : décompte des itérations
	Distances = []

	while not Satisfait and prgs < C.niter: # tant que la contrainte n'est pas satisfaite et que le nombre d'iteration nest pas atteint
		w.reajustWeights(prgs) # reajustement des poids

		Latm = ReturnDisplacedLatmD(Latm,eigenvectors,natm,d,w.weights,C) # deplacement des atomes selon les modes pondérés
		distInit = calcDistance(selection,Latm) # calcul de la surface accessible au solvent après déplacement
		progress = True
		for i,distance in enumerate(distInit): # pour chaque distance entre une paire d'atome
			if distance*aj < distmax[i]*aj: # si une des distane ne progresse pas vers la contrainte
				progress = False # break
				break
		
		if progress: # si tt les distances sont depasse 
			
			w.reajustLimits(prgs) # reajustement des bornes des poids (entre -1 et 1)
			w.saveCombinaisonMax(d) # sauvegarder cette combinaison 
			distmax = copy.deepcopy(distInit)
			Distances = []
			print("Etape {}, évolution de la distance : {}".format(prgs,[aj*(distInit[i]- distDefault[i]) for i in range(len(distInit)) ]))
			# retrieveData(prgs,(distInit- distDefault),C,modeE)
			# modeE+=1
			d+=1
			if d > 20:
				d=10

		else :
			Distances,w = W.watch(distInit,Distances,w)		
			w.precState() # retourne au vecteur de poids précédent si pas d'amélioration	
		# print("Seconds = {:.2f} , {:.2f}, {:.2f}".format((t2-t1),(t3-t2),(((t3-t2)-(t2-t1)))))
		Satisfait = True
		for i in range(len(distInit)): 
			if contrainte[i]*aj > distInit[i]*aj:
				Satisfait = False
				break
		prgs+=1

	return [Latm,w]

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
	if Algo == 'Distance':
		selection = C.selection
		contrainte = C.contr
		Latm,w = goDistance(selection,Latm,contrainte,eigenvectors,len(Latm),C)

	if w.combinaisonMax:
		visual(Latm,w,eigenvectors,pdb,C)
		trajMaxCombinaison(Latm,eigenvectors,w,C)
		writeMemory(w.memory)
	else:
		print('insufficient constraint, please review constraint before starting again')
		print('if you are using Type=Volume you can increase the H value in the Axis file')

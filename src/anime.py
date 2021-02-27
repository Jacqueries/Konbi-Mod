"""Ecrit dans un fichier les vecteurs de basses fréquence filtrées par collectivité
Peut être utilisé seul pour écrire dans des fichiers .pdb les trajectoires selon certains modes
Usage :
	python Animate.py fichier_vecteurs_propres fichier_pdb température seuil nombre_de_frames

fichier_vecteurs_propres: fichier de vecteurs propres et valeurs propres généré par nma.py
fichier_pdb: un fichier au format pdb 
température: un entier. en Kelvin
seuil: seuil de collectivité au dessus duquel on conserve un mode. entre 0 et 1
nombre_de_frames: nombre de frames pour la trajectoire de sortie.

	exemple : python Animate.py eigenvecteur.txt pdbFile.pdb 100 0.2 50

La trajectoire sera écrite au format pdb et peut être visualisée avec le programme VMD par exemple
"""
import sys,os
import mdtraj as md
import numpy as np
import copy
import utilitaires as U

############################################################################################################################
# utilitaires et organisateurs propres a ce fichier

def clean(Vect,natm):
	"""Verification des vecteurs propres, duplicats etc..
		-Args:
			_Vect: liste de vecteurs propres associés a des valeurs propres
			_natm: nombre d'atomes dans le système
	"""
	todel = []
	for key in Vect.keys():
		if len(Vect[key]) != (natm*3):
			print('Vecteur invalide ...')
			print('  nombre d\'éléments du vecteur : {}'.format(len(Vect[key])))
			print('  nombre attendu : {} '.format(natm*3))
			todel.append(key)
		flag = False
		if float(key) < 0:
			print('Fréquence du vecteur négative, supprimé...')
			todel.append(key)
			continue
		for elem in Vect[key]:
			if elem != 0.0: # lié a un bug interne d'écriture du fichier
				break
			else :
				flag = True
		if flag:
			todel.append(key)
	for dl in todel:
		del Vect[dl]
	return Vect

def wanted_mode(Vect,nmode,index):
	"""Permet d'extraire des modes spécifiques de la liste de tout les modes
		-Args:
			_Vect: liste de vecteurs propres associés a des valeurs propres
			_nmodes: nombre de modes désirés en sortie
			_index: index du mode désiré si nmodes = 1
	"""
	order = list(Vect.keys())
	order.sort()
	if index != -1:
		tmp = order[0]
		order[0] = order[index]
		order[index] = tmp
	nVect = {}
	n = 0
	while n < nmode:
		nVect[order[n]]=Vect[order[n]]
		n+=1
	return nVect

def initiate(xyz,Vector,natm,T,nrep,seuil):
	"""Lance les calculs de collectivité d'amplitude et de fréquence.
	Lance le calcul d'une trajectoire. Si collectivité < seuil retourne None
	Retourne une trajectoire du déplacement des coordonnées selon le mode contenu dans Vector
		-Args:
			_xyz: liste de coordonnées des atomes
			_Vector: Dict qui contient un vecteur propre et une valeur propre.   
			_natm: nombre d'atomes dans le système	
			_T: entier, température
			_nrep: nombre de frames pour la trajectoire
			_seuil: float entre 0 et 1. Seuil de collectivité des mouvements
	"""
	Vecteur = calcule_amp_frq(Vector,T)
	trajectoire,distMax = animate2(xyz,Vecteur,nrep,natm)
	key = list(Vecteur.keys())[0]
	coll = collectivity(distMax,Vecteur[key][2],natm)
	if coll < seuil:
		trajectoire = None
	else:
		print("Clé = {}\nCollectivité = {}\nFréquence = {}\nAmplitude = {}".format(key,coll,Vecteur[key][1],Vecteur[key][2]))

	return trajectoire


############################################################################################################################
# calculs

def distance(xyzr,xyzd):
	"""Renvoie les distances au carré entre deux set de coordonnées
	"""
	return ((xyzr[0] - xyzd[0])**2+(xyzr[1] - xyzd[1])**2+(xyzr[2] - xyzd[2])**2)

def calc_frequence(eigv):
	"""Calcule une fréquence associée à la valeur propre eigv
	Référence: Les modes normaux de basse fréquence des protéines Yves-Henri Sanejouand
	"""
	return (np.sqrt(eigv))/(2*np.pi)

def calc_amplitude(frq,tmp):
	"""Calcule l'amplitude associée à la fréquence à une température donnée
	Référence: Les modes normaux de basse fréquence des protéines Yves-Henri Sanejouand
	"""
	kb = 1.98*10**-3 # Kcal/(mol*K)
	return (np.sqrt(2*kb*tmp)/(2*np.pi*frq))

def calcule_amp_frq(Vecteur,Temperature):
	"""Reformate Vecteur pour qu'il contienne la fréquence du mode associé et l'amplitude à une tempérture donnée
		-Args:
			_Vecteur: Dict qui contient un vecteur propre et une valeur propre. 
			_Température: Entier température en kelvin 
	"""
	mode = list(Vecteur.keys())[0]
	frequence = calc_frequence(float(mode))
	amplitude = calc_amplitude(frequence,Temperature)
	Vecteur[mode] = [Vecteur[mode],frequence,amplitude]
	return Vecteur

def collectivity(distMax,amp,natm):
	"""Calcule la collectivité d'un mode
	avec distMax une liste d'amplitude des mouvements de chaque atome, 
	natm le nombre d'atomes et amp l'amplitude du mode
	Référence : Dynamical Properties of the MscL of Escherichia coli: A Normal Mode Analysis
	H.Valadie,J.Lacapc,Y-H.Sanejouand and C.Etchebest

	"""
	tot = sum(distMax)
	kappa = 0
	for i in range(natm):
		alpha = (distMax[i]/tot) # de sorte a ce que la somme des distMax soit égale à 1
		kappa += alpha*np.log(alpha)
	return (1/natm)*np.exp(-kappa)

def deplacement(tr,distMax,xy):
	"""Garde en mémoire les deplacement d'amplitude maximums des atomes et les renvoie dans distMax

		tr = trajectoire courante
		xy = coordonées xyz de référence
		distMax le vecteur de déplacement
	"""
	for i in range(0,len(xy),3):
		disti = distance([xy[i],xy[i+1],xy[i+2]],[tr[i],tr[i+1],tr[i+2]])
		if disti > distMax[int(i/3)]:
			distMax[int(i/3)]=disti
	return distMax

def addition_vecteurs(v1,v2,ope):
	"""Additionne ou soustrait deux vecteurs v1,v2
	ope peut être 1 ou -1
	"""
	if len(v1) != len(v2):
		sys.exit('Err: addition_vecteurs, vecteurs de tailles différentes')
	res = []
	for i in range(len(v1)):
		res.append(v1[i] + ope*v2[i])
	return res

############################################################################################################################
# trajectoire, animation des modes normaux par animate.py SEUL

def animate2(xyz,Vecteur,nrep,natm,way='default'):
	"""Deplace les atomes de la structure selon un mode sur nrep frames
		-Args:
			_xyz: liste de coordonnées des atomes
			_Vecteur: Dict qui contient un vecteur propre et une valeur propre, une fréquence et une amplitude   
			_nrep: nombre de frames pour la trajectoire
			_natm: nombre d'atomes dans le systeme
			_way: module la valeur retournée. Depend du programme qui tourne. 
			animate.py : way=default, Konbi-Mod: way=Free
	"""
	tab = list(range((nrep*2)+1)) # contient les frames de la trajectoire, au centre est l'état initial 
	tab[nrep] = xyz
	distMax = np.zeros(natm)
	for d in range(0,nrep):
		t = d
		tab[d+nrep+1] = addition_vecteurs(tab[d+nrep],deplacement_along_mode2(Vecteur,t,natm),1)
		distMax = deplacement(tab[d+nrep+1],distMax,xyz)
		if way == 'default':
			tab[nrep-d-1] = addition_vecteurs(tab[nrep-d],deplacement_along_mode2(Vecteur,-t,natm),1)
	if way == 'default':
		return [tab,distMax]
	else :
		return distMax

def deplacement_along_mode2(Vect,t,natm):
	"""Deplace toutes les coordonées du système selon un mode
		-Args:
			_Vect: Dict qui contient un vecteur propre et une valeur propre, une fréquence et une amplitude   
			_t: frame du déplacement
			_natm: nombre d'atomes du système
	Référence : Dynamical Properties of the MscL of Escherichia coli: A Normal Mode Analysis
	H.Valadie,J.Lacapc,Y-H.Sanejouand and C.Etchebest
	"""
	dpc = [] # liste des deplacements des coordonnées
	for coor in range(natm*3): # pour chaque coordonnée du système
		dpc.append(0)
		mode = list(Vect.keys())[0]
		# Vect[mode][0][coor] = ieme(coor) coordonnée du mode, Vect[mode][2] = amplitude du mode
		# Vect[mode][1] = fréquence du mode
		dpc[coor] = Vect[mode][0][coor]*Vect[mode][2]*np.cos(2*np.pi*Vect[mode][1]*t+np.pi)
		dpc[coor] = 1/(np.sqrt(1.55))*dpc[coor] 
		# pondération par la masse moyenne de l'atome. m(C,H,N,O)=1.55)
		# inéxact mais ce module ne devrait pas être employé pour visualiser les modes, seulement dans le cadre de Konbi-Mod 
	return dpc

############################################################################################################################
# starters et main(s)

def returnLowFreq(xyz,Vector,natm,T,nrep,threshold):
	"""Conserve les mode de collectivité supérieure au seuil
		-Args:
			_xyz: coordonées du système
			_Vector: vecteur propre
			_natm: nombre d'atomes du système
			_T: température indiquée dans le fichier de configuration
			_nrep: nombre de frames pour la trajectoire (pour le calcul de collectivité)
			_threshold: seuil de collectivité
	"""

	Vecteur = calcule_amp_frq(Vector,T)
	distMax = animate2(xyz,Vecteur,nrep,natm,'Free')
	key = list(Vecteur.keys())[0]
	amp = Vecteur[key][2]
	freq = Vecteur[key][1]
	co = collectivity(distMax,amp,natm)
	Vecteur[key] = [Vecteur[key][0],freq,amp,co]
	print("Clé = {}\nCollectivité = {}\nFréquence = {}\nAmplitude = {}".format(key,co,freq,amp))
	if co >= threshold:
		print('Mode conservé')
		return Vecteur
	else:
		print('Collectivité trop basse, mode ignoré')
		return False

def mainPlay(xyz,T,nrep,threshold,C):
	"""Ecrit dans un fichier les vecteurs de basses frequences et haute collectivité pour usage ultérieur
		-Args:
			_xyz: coordonées du système
			_T: température indiquée dans le fichier de configuration
			_nrep: nombre de frames pour la trajectoire (pour le calcul de collectivité)
			_threshold: seuil de collectivité
			_C: configuration
	"""
	natm = int(len(xyz)/3)
	vects = U.readEigenVs(C.efile)
	vects = clean(vects,natm)

	lowFreq = [] # porte mal son nom,, conserve les index des modes à ignorer
	for i,k in enumerate(vects.keys()):
		vect = {}
		vect[k] = copy.deepcopy(vects[k])
		vect = returnLowFreq(xyz,vect,natm,T,nrep,threshold) # contient les modes de collectivité > seuil
		if not vect:
			lowFreq.append(k)
		else :
			a = list(vect.values())[0]	
			vects[k] = [a[0],a[1],a[2],a[3]]

	for l in lowFreq:
		del vects[l]

	for i,key in enumerate(vects.keys()):
		U.writeSelectedModes(vects[key],key,i,C)

def fromLaunch(Config):
	"""Lance le programme depuis Konbi-Mod
		-Args:
			_Config: configuration
	"""
	T = Config.temp # température
	nrep = 50 # nombre de frames  
	thre = Config.coll
	pdb = md.load_pdb(Config.pdb)
	coord = pdb[0].xyz[0]
	xyz = [c*10 for l in coord for c in l] # conversion en angström
	mainPlay(xyz,T,nrep,thre,Config)

def main(file,pdbF,xyz,T,nrep,seuil):
	"""Hors programme principal (propre a anime.py). permet d'obtenir les trajectories de certains modes 
	et d'ensuite les visualiser
	Ecrit les trajectoires dans le format pdb
		-Args:
			_file: fichier contenant les vecteurs propres et valeurs propres générées par nma.py 
			(fichier primaire de Konbi-Mod.py, peut être récupéré dans le repertoire mes_Modes après une run)
			_pdbF: fichier pdb
			_xyz: liste de coordonnées des atomes
			_T: entier, température
			_seuil: float entre 0 et 1. Seuil de collectivité des mouvements
	"""
	natm = int(len(xyz)/3)
	vects = U.readEigenVs(file)
	vects = clean(vects,natm) # vérifications il existe parfois des duplicats ou des bugs

	for i,k in enumerate(vects.keys()):
		vect = wanted_mode(vects,1,i)
		trajectoire = initiate(xyz,vect,natm,T,nrep,seuil)
		if trajectoire != None:
			print("Writing.........")
			U.write_PDB(pdbF,trajectoire,'w',i)

if __name__ == '__main__':
	vecteur = sys.argv[1]
	pdbfile = sys.argv[2]
	temperature = int(sys.argv[3])
	seuil = float(sys.argv[4])
	nframes = int(sys.argv[5])
	pdb = md.load_pdb(pdbfile) # lecture du fichier pdb par le module mdtraj
	coord = pdb[0].xyz[0]
	xyz = [c*10 for l in coord for c in l] # conversion des coordonnées dans un autre format et passage des nm aux angström
	main(vecteur,pdbfile,xyz,temperature,nframes,seuil)


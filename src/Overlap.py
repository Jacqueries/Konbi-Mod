"""Overlap.py permet d'identifier les modes responsables de la transition entredeux états conformationels.
=====================================================
Usage:
	python3 Overlap.py structure_départ.pdb structure_arrivée.pdb vecteurs.txt

prend en entrée deux fichiers pdb et un fichier texte contenant les vecteurs propres
"""
import sys,os
import utilitaires as U
import NMAContrainte as N
import mdtraj as md
import numpy as np


def dist_between_two_strucOld(s1,s2):
	"""Calcule un vecteur qui contient les distances entre les mêmes atomes de deux structures
		-Args:
			_s1: liste des coordonnées de la première structure
			_s2: liste des coordonnées de la deuxième structure
	"""
	distances = []
	for i in range(0,len(s1),3):
		distances.append(N.dist_atomes( [s1[i],s1[i+1],s1[i+2]],[s2[i],s2[i+1],s2[i+2]]))
	return distances

def dist_between_two_struc(s1,s2):
	"""Calcule un vecteur qui contient les distances entre les mêmes atomes de deux structures
		-Args:
			_s1: liste des coordonnées de la première structure
			_s2: liste des coordonnées de la deuxième structure
	"""
	distances = []
	for i in range(len(s1)):
		diff = s1[i]-s2[i]
		distances.append(diff)
	return np.array(distances)

def verify_transform(pdb1,pdb2,lenVector):
	"""Verifie que c'est la même protéine dans les deux fichiers pdb puis revoie leur coordonnées
	sous forme de deux listes
		-Args:
			_pdb1: fichier au format pdb
			_pdb2: fichier au format pdb
	"""
	pdb_1 = md.load_pdb(pdb1) # lecture du fichier pdb par le module mdtraj
	pdb_2 = md.load_pdb(pdb2)
	if pdb_1.n_atoms != pdb_2.n_atoms or pdb_1.n_residues != pdb_2.n_residues:
		print("Pas le même nombre d'atomes/residues entre les deux structures !!")
		print("Nombre d'atomes : {} vs {}".format(pdb_1.n_atoms,pdb_2.n_atoms))
		print("Nombre de résidues : {} vs {}".format(pdb_1.n_residues,pdb_2.n_residues))
		sys.exit()
	if pdb_1.n_atoms*3 != lenVector:
		print("Nombre de coordonnées du système: {} vs nombre de coordonnées normales{}".format(pdb_1.n_atoms*3,lenVector))
		print("Pas le même nombre de coordonnées normales et de coordonnées atomiques !!")
		sys.exit()
	xyz1 = [c*10 for l in pdb_1[0].xyz[0] for c in l] # conversion des coordonnées dans un autre format et passage des nm aux angström
	xyz2 = [c*10 for l in pdb_2[0].xyz[0] for c in l]

	return (xyz1,xyz2)


def calcOverlap(modes,Vdistance):
	"""Prend un vecteur de distance entre deux structures en entrée, un dictionnaires de vecteurs propres
	et renvoi un vecteur de valeur qui correspond a l'overlap de chaque mode
		-Args:
			_modes: dictionnaire de modes normaux (vecteurs propres) de basse fréquence
			_Vdistance: Vecteur de distance entre toutes les coordonnées des atomes identiques de deux états conformationnels 
	"""
	Overlap = []
	for mode in modes.keys():
		Overlap.append(np.vdot(modes[mode],Vdistance))
	return Overlap

if __name__ == '__main__':

	if len(sys.argv) == 2 and sys.argv[1] == 'help':
		print(__doc__)
		sys.exit()

	pdbD = sys.argv[1]
	pdbA = sys.argv[2]
	eigenFile = sys.argv[3]

	vectors = U.readEigenVs(eigenFile)
	lenVector = len(list(vectors.values())[0])
	s1,s2 = verify_transform(pdbD,pdbA,lenVector)
	VectDistances = dist_between_two_struc(s1,s2)
	norme = np.linalg.norm(VectDistances)
	VectDistUnitaire = VectDistances/norme
	Overlap = calcOverlap(vectors,VectDistUnitaire)
	print(sorted(Overlap))
	print(Overlap)
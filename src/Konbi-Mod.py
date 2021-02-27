"""Konbi-Mod permet de trouver une combinaison de modes normaux de basses fréquences
d'une protéine qui satisfait une contrainte donnée par l'utilisateur.
"""

__auteur__ = ("Ilyas Grandguillaume")
__contact__ = ("ilyasgrand@live.fr")

import mdtraj as md
from nma import ANMA
import sys,os
import utilitaires as U

configFile = sys.argv[1]
Config = U.readConfig(configFile)
U.verif_env(Config)

if Config.generateNMA:

	Config.efile = '{}/EignenVector{}.txt'.format(Config.eigenFolder,Config.pdbName) # crée le fichier qui contiendra les vecteurs propres
	
	##################################################################################################################
	# calcul des modes normaux par le module mdtraj et ecriture dans Config.efile
	pdb = md.load_pdb(Config.pdb)
	anma = ANMA(mode=50, rmsd=0.06, n_steps=50, selection='all')
	anma_traj = anma.fit_transform(pdb,Config.efile)
	###################################################################################################################
	
	# Config.mfile = '{}\\ModesBassesFrequences{}.txt'.format(Config.eigenFolder,Config.pdbName)
	# A.fromLaunch(Config)
	# N.fromLaunch(Config)

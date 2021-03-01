# Konbi-Mod
Konbi-Mod permet de trouver une combinaison de modes normaux qui satisfait une contrainte, par exemple d'augmenter la surface accessible au solvant d'une cavité, d'ourvrir un canal, ou bien d'augmenter la distance entre certains atomes.

Ce programme à été crée dans le cadre d'un projet encadré par un professeur à l'université paris diderot, et à été développé en python 3.8.

# Prérequis:
	-python3 ou supérieur

	-GNU Compiler Collection gcc

	-module python mdtraj: https://mdtraj.org/1.9.4/installation.html
		$ conda install -c conda-forge mdtraj

	-module python nma : https://github.com/mdtraj/nma
		$ git clone https://github.com/mdtraj/nma.git && cd nma
		$ python setup.py install

	-module python networkx : https://networkx.org
		$ conda install networkx
	
	-module python freesasa : https://github.com/freesasa/freesasa-python
		$ pip install freesasa

# Utilisation

		$python configuration.txt

Konbi_Mod nécessite un fichier de configuration pour fonctionner. Un exemple de ce fichier est donné dans ce repertoire. Konbi-Mod nécessite aussi un répertoire qui contient ./src ./Struct ./Modes.\newline
Dans le fichier de configuration vous pouvez indiquer si vous voulez que le script calcule les modes de basses fréquences avec :  
		-GenerateVectors=YES  
Une fois calculés pour un système donné pensez, pour gagner du temps, à indiquer dans le fichier de configuration:  
		-GenerateVectors=NO  
Les vecteurs et valeurs propres calculés seront stockés dans un fichier primaire brut dans le répertoire Modes, sous le nom EigenVectorPDB.txt avec PDB le nom du fichier pdb sur lequel vous calculez les modes.  Un fichier secondaire contenant les modes de collectivité supérieure à un seuil sera crée à partir de EigenVectorPDB.txt, formaté différemment.  
Vous pouvez ensuite indiquer dans le fichier de configuration le chemin vers le fichier primaire avec:  
		-EigenFile=./Modes/EignenVectorPDB.txt  
et vers le secondaire avec:  
		ModesFile=./Modes/ModesBassesFrequencesPDB.txt  
Pour une run avec le même seuil de collecivité, seul le fichier ModesBassesFrequencesPDB.txt est nécessaire, pour des runs avec un seuil de collectivité différent, il est necessaire de le génerer à nouveau à partir du fichier primaire, pour cela indiquez  
		ModesFile=NONE

Le reste du fichier de configuration permet d'indiqueer au programme les elements suivants:  
	-Le chemin vers le fichier pdb sur lequel on désire travailler  
		PathToPdbFile=./Struct/PDB.pdb  
	-Le nom du repertoire ou le programme écrira EignenVectorPDB.txt et ModesBassesFrequencesPDB.txt  
		FolderVectors=Modes  
	-Le chemin vers l'axe de symmétrie de la protéine  
		PathToAxis=./Struct/Axis/Axis_PDB.txt  
Ceci n'est nécessaire que quand la contrainte est un volume à atteindre, sinon ce paramètre est ignoré.  
	-Le type de contrainte voulue: Distance, Volume ou Surface  
		Type=Distance  
	-La synthaxe de la selection qui dépend de la contrainte: pour une contrainte de distance entre  atomes, indiquer une liste de paires de numéros d'atomes du fichier .pdb, par exemple:  
		Selection=[[647,1691]]  
pour une contrainte de surface accessible au solvant, la synthaxe complète est disponible ici: https://freesasa.github.io/1.1/Selection.html . Par exemple pour les résidus 366 à 403 de la chaine A:  
		Selection=sel, resn 366-403  
pour une contrainte sur le volume, aucune selecion n'est necessaire.  
	-Le ratio d'augmentation ou de diminution de la contrainte. Par exemple 0.95 pour une diminution de 5%, 1.05 pour une augmentation de 5%  
		Contrainte=1.05  
	-Le seuil de collectivité de mouvements des modes pour la génération des modes normaux:  
		Collectivity=0.25  
	-La température en Kelvin pour le calcul d'amplitude des modes:  
		Temperature=300  
	-Le nombre d'itération lors de l'optimisation:  
		Nombre_iterations=500  

Dans le fichier de config les # sont des lignes commentaires. L'ordre des elements n'a pas d'importance.

Astuce : Pour un fichier pdb donné, faite une première run avec vos paramètres de collectivité, de température et les paramètres par défaut (dans configdefaut) et spécifiez uniquement GenerateVectors=YES. Après la génération des deux fichiers, indiquez GenerateVectors=NO et spécifiez l'emplacement de vos deux fichiers de vecteurs propres. Vous pouvez ensuite essayer différents types de contrainte, de selection, d'itération et de température.

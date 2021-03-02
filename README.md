# Konbi-Mod
Konbi-Mod permet de trouver une combinaison de modes normaux qui satisfait une contrainte, par exemple d'augmenter la surface accessible au solvant d'une cavité, d'ouvrir un canal, ou bien d'augmenter la distance entre certains atomes.

Ce programme à été écrit dans le cadre d'un projet encadré par un professeur à l'université paris diderot, et à été développé en python 3.8.

# Prérequis
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

		$python ./src/Konbi-Mod.py configuration.txt

Konbi_Mod nécessite un fichier de configuration en entrée pour fonctionner. Un exemple de ce fichier est donné dans ce répertoire. Konbi-Mod nécessite aussi un répertoire qui contient ./src ./Struct ./Modes.\newline
Dans le fichier de configuration vous pouvez indiquer si vous voulez que le script calcule les modes de basses fréquences avec :

	GenerateVectors=YES  

Une fois calculés pour un système donné pensez, pour gagner du temps, à indiquer dans le fichier de configuration:  

	GenerateVectors=NO  

Les vecteurs et valeurs propres calculés seront stockés dans un fichier primaire brut dans le répertoire Modes, sous le nom EigenVectorPDB.txt avec PDB le nom du fichier pdb sur lequel vous calculez les modes.  Un fichier secondaire contenant les modes de collectivité supérieure à un seuil sera créé à partir de EigenVectorPDB.txt, formaté différemment.  
Vous pouvez ensuite indiquer dans le fichier de configuration le chemin vers le fichier primaire avec:  

	EigenFile=./Modes/EignenVectorPDB.txt  

et vers le secondaire avec:  

	ModesFile=./Modes/ModesBassesFrequencesPDB.txt  

Pour une run avec le même seuil de collecivité, seul le fichier ModesBassesFrequencesPDB.txt est nécessaire, pour des runs avec un seuil de collectivité différent, il est nécessaire de le générer à nouveau à partir du fichier primaire, pour cela indiquez  

	ModesFile=NONE

Le reste du fichier de configuration permet d'indiquer au programme les éléments suivants:  

	-Le chemin vers le fichier pdb sur lequel on désire travailler  
		PathToPdbFile=./Struct/PDB.pdb  

	-Le nom du répertoire ou le programme écrira EignenVectorPDB.txt et ModesBassesFrequencesPDB.txt  
		FolderVectors=Modes  

	-Le chemin vers l'axe de symétrie de la protéine  
		PathToAxis=./Struct/Axis/Axis_PDB.txt  

Ceci n'est nécessaire que quand la contrainte est un volume à atteindre, sinon ce paramètre est ignoré. L'axe est généré automatiquement si rien n'est indiqué.

	-Le type de contrainte voulue: Distance, Volume ou Surface  
		Type=Distance  

	-La syntaxe de la sélection qui dépend de la contrainte: pour une contrainte de distance entre  atomes, indiquer une liste de paires de numéros d'atomes du fichier .pdb, par exemple:   
		Selection=[[647,1691]]  

Pour une contrainte de surface accessible au solvant, la syntaxe complète est disponible ici: https://freesasa.github.io/1.1/Selection.html . Par exemple pour les résidus 366 à 403 de la chaîne A:  
		
		Type=Selection
		Selection=sel, resn 366-403  

Pour une contrainte de volume, aucune sélection n'est nécessaire.  

	-Le ratio d'augmentation ou de diminution de la contrainte. Par exemple 0.95 pour une diminution de 5%, 1.05 pour une augmentation de 5%  
		Contrainte=1.05  

	-Le seuil de collectivité de mouvements des modes pour la génération des modes normaux:  
		Collectivity=0.25  

	-La température en Kelvin pour le calcul d'amplitude des modes:  
		Temperature=300  

	-Le nombre d'itération lors de l'optimisation:  
		Nombre_iterations=500  *
(le nombre d'itération doit être supérieur à deux fois le nombre de modes)  

Dans le fichier de config les # sont des lignes commentaires. L'ordre des elements n'a pas d'importance.

*Astuce :* Pour un fichier pdb donné, faites une première run avec vos paramètres de collectivité, de température et les paramètres par défaut (dans configTemplate.txt) et spécifiez uniquement GenerateVectors=YES. Après la génération des deux fichiers, indiquez GenerateVectors=NO et spécifiez l'emplacement de vos deux fichiers de vecteurs propres. Vous pouvez ensuite essayer différents types de contraintes, de sélection, d'itération et de température.

# Sortie

Konbi-Mod génère en premier lieu les fichiers contenant les modes de basses fréquences (50 modes, non filtré = fichier primaire) et les modes de basses fréquences au dessus du seuil de collectivité indiqué (n modes filtrés = fichier secondaire) dans le répertoire ./Modes.

Si la contrainte choisie est Type=Volume pour une protéine à canal, Konbi-Mod va générer un axe de symétrie qui va passer le long de ce canal. Cet axe sera écrit dans un fichier situé dans le répertoire ./Struct/Axis. Ensuite le pavé droit sera écrit dans un fichier dans ./Struct sous le nom OBB.pdb. Sa position pourra être visualisée en le chargeant avec la structure de la protéine dans un logiciel de visualisation moléculaire.   

Ensuite lors de l'optimisation le fichier memoryWeights.csv sera généré dans ./out. Il contient la valeur des poids à chaque itération.

Après l'optimisation deux fichier pdb seront générées, l'un contenant la trajectoire selon la combinaison de modes la plus favorable par rapport à la contrainte donnée, et l'autre contenant la frame *0* du système à chaque étape de l'optimisation ou la valeur de la contrainte à été améliorée (ces étapes sont notifiées sur le terminal pendant l'optimisation). Un fichier Weights.txt donne les poids attribués au modes à chacune de ces étapes. tous ces fichiers seront générées dans ./out.

# Exemple

		$python ./src/Konbi-Mod.py ./Exemples/configuration3pgk.txt

La protéine 3pgk à deux domaines articulés autour d'une charnière. Les déplacements le long de ses modes sont amples et c'est un modèle idéal pour éprouver Konbi-Mod. La contrainte définie est une augmentation de la surface accessible au solvant de 20% des acides aminés 366 à 403 qui constituent la région  charnière de la protéine.

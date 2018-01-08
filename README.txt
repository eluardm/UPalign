Matthias ELUARD
M2 bib Paris Diderot

UPalign est un outil d'alignement structural flexible basée sur l'alignement des Unités Protéiques.
Dépendances : TM-align, R

Pour chacunes des 2 protéines à aligner :
1 - Générer des Unités Protéiques avec Protein Peeling. http://www.dsimb.inserm.fr/dsimb_tools/peeling3/
2 - Copier les resulats dans un dossier portant le nom du fichier PDB dans le répertoire /results/
3 - Exécuter le script dpuis la racine du répertoire ./UPalign/
4 - Les resulats se trouvent dans les dossier "date+protA+protB" et "date+protB+protA" dans /results

Exemple de commande :
./src/upalign.py "./data/1a03.pdb" "./data/1cnp.pdb"




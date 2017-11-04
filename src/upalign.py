#!/usr/bin/python3
# coding: utf-8

DICOAA = {
	'GLU':'E',
	'ASP':'D',
	'ALA':'A',
	'ARG':'R',
	'ASN':'N',
	'CYS':'C',
	'GLN':'Q',
	'GLY':'G',
	'HIS':'H',
	'ILE':'I',
	'LEU':'L',
	'LYS':'K',
	'MET':'M',
	'PHE':'F',
	'PRO':'P',
	'SER':'S',
	'THR':'T',
	'TRP':'W',
	'TYR':'Y',
	'VAL':'V'
}

class Protein:
	pdb = []
	pdbTmp = []
	seq = ""

	def __init__(self, filePath):
		self.import_pdb(filePath)
		self.read_seq()

	def import_pdb(self, filePath):
		try:
			with open(filePath, 'r') as fichier:
				self.pdb = fichier.readlines()
		except IOError as erreur:
			print(erreur)
			exit(1)

	def read_seq(self):
		#Ne lit que la chaine A pour le moment.
		numRes = 0
		for ligne in self.pdb:
			if ligne[0:3] == 'TER':
				break
			if ligne[0:4] == 'ATOM' and int(ligne[22:26]) != numRes:
				numRes = int(ligne[22:26])
				self.seq += DICOAA[ligne[17:20]]

if __name__ == '__main__':
	protA = Protein("../data/4dmo/4dmo.pdb")
	protB = Protein("../data/2pqt/2pqt.pdb")



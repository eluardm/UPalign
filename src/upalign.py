#!/usr/bin/python3
# coding: utf-8

import os
from os.path import basename
import glob 
from copy import deepcopy
import re
import numpy as np

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
    id = ''
    path = ''
    pdb = []
    pdbTmp = []
    seq = ""

    def __init__(self, filePath):
        fileName, fileExtension = os.path.splitext(filePath)
        #if(fileExtension!='.pdb'):
        #    print("Erreur extention .pdb")
        #    exit(1)
        self.path = filePath
        self.id = basename(fileName)
        self.pathTmp = "../tmp/tmp_"+self.id+".pdb"
        self.import_pdb(filePath)
        self.read_seq()

    def import_pdb(self, filePath):
        '''Charge le fichier PDB en mémoire, ne conserve que les lignes ATOM, TER et END.
        Uniquement la chaine A'''
        try:
            with open(filePath, 'r') as fichier:
                coordtmp = []
                for ligne in fichier:
                    key = ligne[0:6].strip()
                    if key == 'ATOM':
                        self.pdb.append(ligne)
                        coordtmp.append([float(ligne[30:38]),float(ligne[38:46]),float(ligne[46:54])])
                    if key in ['TER','END']:
                        self.pdb.append(ligne)
                        break
                self.coord = np.asarray(coordtmp)
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

    def get_pdb(self):
        return "".join(self.pdb)

    def write_pdbtmp(self):
        with open("../tmp/tmp_"+self.id+".pdb","w") as fOut:
            fOut.write(self.get_pdb())

    def rotate(self, matrix):
        newCoord = np.zeros(self.coord.shape)
        newCoord[:,0]=matrix[0,0]+matrix[0,1]*self.coord[:,0]+matrix[0,2]*self.coord[:,1]+matrix[0,3]*self.coord[:,2]
        newCoord[:,1]=matrix[1,0]+matrix[1,1]*self.coord[:,0]+matrix[1,2]*self.coord[:,1]+matrix[1,3]*self.coord[:,2]
        newCoord[:,2]=matrix[2,0]+matrix[2,1]*self.coord[:,0]+matrix[2,2]*self.coord[:,1]+matrix[2,3]*self.coord[:,2]
        self.coord = deepcopy(newCoord)

class Alignement:
    #num = -1
    #UP = (-1,-1) #start et end

    def __init__(self, outAlign, i):
        #self.num = i
        self.list_alignedA = []
        self.list_alignedB = []
        lignes = outAlign.split('\n')
        #self.length = int(lignes[16].split(',')[0].split('=')[-1])
        #res = re.search("in_PU"+str(i)+"_(?P<start>\d+)-(?P<end>\d+)",lignes[11])
        #if res is not None:
        #    self.UP = int(res.group('start')),int(res.group('end'))
        #else:
        #    print("erreur sortie tm-align")
        #    exit()

        seqA = list(lignes[22])
        motif = list(lignes[23])
        seqB = list(lignes[24])

        tailleMin = 6
        i = 0
        lim = len(motif)
        while i < lim:
            length = 1
            if motif[i] in ['.',':']:
                sA = len(seqA[:i+1])-seqA[:i+1].count('-')
                #+self.UP[0]-1
                sB = len(seqB[:i+1])-seqB[:i+1].count('-')
                while(i+length < lim and motif[i+length] in ['.',':']):
                    length+=1
                eA = sA + length - 1
                eB = sB + length - 1
                if(length > tailleMin):
                    self.list_alignedA += range(sA,eA+1)
                    self.list_alignedB += range(sB,eB+1)
            i += length


def remove_aligned(atom, aligned):
    pdbtmp = []
    i = 0
    numResPrec = -1
    for ligne in atom:
        numRes = int(ligne[22:26])
        if(numRes != numResPrec):
            i+=1
        if i not in aligned:
            pdbtmp.append(ligne)
        numResPrec = numRes
    return pdbtmp



def recup_pdb_pu(id, i):
    file1 = glob.glob('../results/{}/*'.format(id))
    file2 = [fichier for fichier in file1 if(basename(fichier).startswith('in_PU'+str(i)) and not basename(fichier).endswith('.png'))]
    return file2

def import_matrix(fileMatrix):
    with open(fileMatrix) as fMat:
        matrix = np.zeros((3,4))
        for ligne in fMat:
            if ligne.startswith(' 1') or ligne.startswith(' 2') or ligne.startswith(' 3'):
                matrix[int(ligne.split()[0])-1,:] = float(ligne.split()[1]),float(ligne.split()[2]),float(ligne.split()[3]),float(ligne.split()[4])
    return deepcopy(matrix)

if __name__ == '__main__':
    protB = Protein("../data/1e0s.pdb")
    protA = Protein("../data/2j5x.pdb")
    #protA.write_pdbtmp()
    
    

    #Creer les rep pour stocker les fichiers de peeling
    #os.system("mkdir ../results/"+protA.id)
    #os.system("mkdir ../results/"+protB.id)
    i = 0
    limite = 10 #securité while

    while(i<limite):
        protB.write_pdbtmp()
        protBtmp = deepcopy(protB)

        listePU = recup_pdb_pu(protA.id, i)
        if(not listePU):
            #liste vide si true
            break

        for pu in listePU:

            sortie = os.popen("TMalign {} {} -m ../tmp/matrix.txt".format(pu,protBtmp.pathTmp), "r").read()
            print(sortie)
            a = Alignement(sortie, i)
            matRot = import_matrix("../tmp/matrix.txt")
            print("-------------------------------------")
            print(a.list_alignedA)
            print(a.list_alignedB)
            print("-------------------------------------")

            puTmp = Protein(pu)
            puTmp.rotate(matRot)
            protBtmp.pdb = remove_aligned(protBtmp.pdb, a.list_alignedB)

            protBtmp.write_pdbtmp()

        i+=1







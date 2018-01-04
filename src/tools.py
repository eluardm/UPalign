#!/usr/bin/python3
# coding: utf-8

import os
from os.path import basename
import glob 
from copy import deepcopy
import re
import numpy as np
from time import gmtime, strftime
import sys
import logging

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
    seq = ""

    def __init__(self, filePath):
        self.pdb = []
        self.numRes = []
        pdbTmp = []
        self.coord = None
        coordtmp = []
        fileName, fileExtension = os.path.splitext(filePath)
        #if(fileExtension!='.pdb'):
        #    print("Erreur extention .pdb")
        #    exit(1)
        self.path = filePath
        self.id = basename(fileName)
        self.pathTmp = "./tmp/tmp_"+self.id+".pdb"
        self.import_pdb(filePath)
        self.read_seq()

    def import_pdb(self, filePath):
        '''Charge le fichier PDB en mémoire, ne conserve que les lignes ATOM, TER et END.
        Uniquement la chaine A'''
        try:
            with open(filePath, 'r') as fichier:
                coordtmp = []
                num = 0
                for ligne in fichier:
                    key = ligne[0:6].strip()
                    if key == 'ATOM':
                        self.pdb.append(ligne)
                        num = int(ligne[22:26].strip())
                        if(num not in self.numRes):
                            self.numRes.append(num)
                        coordtmp.append([float(ligne[30:38]),float(ligne[38:46]),float(ligne[46:54])])
                    if key in ['TER','END']:
                        self.pdb.append(ligne)
                        break
                self.coord = np.asarray(coordtmp)
        except IOError as erreur:
            logging.error(erreur)
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
        with open("./tmp/tmp_"+self.id+".pdb","w") as fOut:
            fOut.write(self.get_pdb())

    def rotate(self, matrix):
        newCoord = np.zeros(self.coord.shape)
        newCoord[:,0]=matrix[0,0]+matrix[0,1]*self.coord[:,0]+matrix[0,2]*self.coord[:,1]+matrix[0,3]*self.coord[:,2]
        newCoord[:,1]=matrix[1,0]+matrix[1,1]*self.coord[:,0]+matrix[1,2]*self.coord[:,1]+matrix[1,3]*self.coord[:,2]
        newCoord[:,2]=matrix[2,0]+matrix[2,1]*self.coord[:,0]+matrix[2,2]*self.coord[:,1]+matrix[2,3]*self.coord[:,2]
        self.coord = deepcopy(newCoord)

    def remove_aligned(self, aligned):
        pdbtmp = []
        coordtmp = []
        i = 0
        numResPrec = -1
        for ind, ligne in enumerate(self.pdb):
            numRes = int(ligne[22:26])
            if(numRes != numResPrec):
                i+=1
            if i not in aligned:
                pdbtmp.append(ligne)
                if ligne[0:4] == 'ATOM':
                    coordtmp.append(self.coord[ind])
            numResPrec = numRes
        self.pdb = deepcopy(pdbtmp)
        self.coord = deepcopy(np.asarray(coordtmp))

class Alignement:

    def __init__(self, outAlign, i):
        self.list_alignedA = []
        self.list_alignedB = []
        lignes = outAlign.split('\n')
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
                sB = len(seqB[:i+1])-seqB[:i+1].count('-')
                while(i+length < lim and motif[i+length] in ['.',':']):
                    length+=1
                eA = sA + length - 1
                eB = sB + length - 1
                if(length > tailleMin):
                    self.list_alignedA += range(sA,eA+1)
                    self.list_alignedB += range(sB,eB+1)
            i += length
        self.long = len(self.list_alignedA)

def recup_pdb_pu(id, i):
    file1 = glob.glob('./results/{}/*'.format(id))
    file2 = [fichier for fichier in file1 if(basename(fichier).startswith('in_PU'+str(i)) and not basename(fichier).endswith('.png'))]
    return file2

def import_matrix(fileMatrix):
    with open(fileMatrix) as fMat:
        matrix = np.zeros((3,4))
        for ligne in fMat:
            if ligne.startswith(' 1') or ligne.startswith(' 2') or ligne.startswith(' 3'):
                matrix[int(ligne.split()[0])-1,:] = float(ligne.split()[1]),float(ligne.split()[2]),float(ligne.split()[3]),float(ligne.split()[4])
    return deepcopy(matrix)

def sumTermeTM(protA, protB, align):
    i = 0

    lTargetA = len(protA.seq)
    lTargetB = len(protB.seq)

    sumTmA = 0
    sumTmB = 0

    d0A = 1.24*(lTargetA-15)**(1.0/3)-1.8
    d0B = 1.24*(lTargetB-15)**(1.0/3)-1.8
    
    caA = get_coordCA(protA.pdb,protA.coord, align.list_alignedA)
    caB = get_coordCA(protB.pdb,protB.coord, align.list_alignedB)

    while(i < align.long):
        dist = dist3D(caA[i,:], caB[i,:])
        sumTmA += 1/(1+(dist/d0A)**2) 
        sumTmB += 1/(1+(dist/d0B)**2) 
        i += 1
    return sumTmA, sumTmB

def tm_score(somme, l):
    return(somme/l)

def get_coordCA(pdb,coord, res):
    i = 1
    ca = []
    cp = 0
    for inx, ligne in enumerate(pdb):
        if ligne[0:4] == 'ATOM' and ligne[12:16].strip() == 'CA':
            if(i in res):
                cp +=1
                ca.append(coord[inx])
            i += 1
    return deepcopy(np.asarray(ca))

def dist3D(coord1, coord2):
    return(np.linalg.norm(coord1-coord2))

def write_pdb_align(fileName, prot, lastAtom = 0, lastChain = 64, align = None):
    with open(fileName,"a") as fOut:
        lastAtomTmp = 0
        for ind, ligne in enumerate(prot.pdb):
            if(ligne[0:4] == 'ATOM'):
                if(align == None or prot.numRes.index(int(ligne[22:26]))+1 in align):
                    fOut.write(ligne[0:6]+\
                        "{:5d}".format(int(ligne[6:11])+lastAtom)+\
                        ligne[11:21]+\
                        chr(lastChain+1)+\
                        ligne[22:30]+\
                        "{:8.3f}{:8.3f}{:8.3f}".format(prot.coord[ind][0],prot.coord[ind][1],prot.coord[ind][2])+\
                        ligne[54:])
                    lastAtomTmp = int(ligne[6:11])
            else:
                fOut.write(ligne)
    return(lastAtomTmp, lastChain+1)

def file_plot_score(level, tmA, tmb, outPath):
    with open(outPath+"plotTMscore.dat","a") as fOut:
        if(level == 0):
            fOut.write("Level\tTM-scoreA\t TM-scoreB\n")
        fOut.write("{}\t{}\t{}\n".format(level, tmA, tmb))
            

def file_plot_align(aligned, level, numPu, xlim, numRes, numResInit, outPath):
    '''
    numResInit = numero des résidu dans le pdb initial car le num des residu dans le pdb des pu n'est pas le mème car résidus manquant ignorés.
    Le plot des résidus alignés doit se faire selon la numérotation initial.
    '''
    with open(outPath+"plotAligned"+str(level+1)+".dat","a") as fOut:
        if(fOut.tell() == 0):
            fOut.write("Residus\tNumPU\tLevel\n")
        for ele in aligned:
            #NumPU+1 car coloration selon numPU en R -> or, 0 = blanc
            fOut.write(str(numResInit[numRes[ele-1]-1])+"\t"+str(numPu+1) + "\t"+ str(level+1)+"\n")

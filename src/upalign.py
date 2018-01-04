#!/usr/bin/python3
# coding: utf-8

import os
from os.path import basename
import glob 
from copy import deepcopy
import re
import numpy as np
from time import gmtime, strftime, sleep
import sys
import logging
from tools import *

def upalign(fileA, fileB):

    protB = Protein(fileB)
    protA = Protein(fileA)

    #Creer les rep pour stocker les fichiers de peeling
    #os.system("mkdir ./results/"+protA.id)
    #os.system("mkdir ./results/"+protB.id)
    i = 0
    limite = 10 #securité while
    


    outPath = "./results/"+strftime("%H%M%S", gmtime())+protA.id+protB.id+"/"
    os.system("mkdir "+outPath)
    logging.basicConfig(filename=outPath+'output.log',level=logging.DEBUG)

    #TM-score initial
    sortie = os.popen("TMalign {} {}".format(fileA,fileB), "r").read()
    logging.info(sortie)
    sortie = sortie.split('\n')
    file_plot_score(0, float(sortie[17][9:17]), float(sortie[17][9:17]), outPath)

    while(i<limite):
        #Pourcentage de residu aligné au final
        pcA = 0
        pcB = 0
        
        listePU = recup_pdb_pu(protA.id, i)
        if(not listePU):
            #liste vide si true
            break

        protBtmp = deepcopy(protB)
        protBtmp.write_pdbtmp()
        lastAtom, lastChain = write_pdb_align(outPath+"level_"+str(i+1)+".pdb", protBtmp)

        (sumTmA, sumTmB) = (0,0)
        for numPu, pu in enumerate(listePU):

            sortie = os.popen("TMalign {} {} -m ./tmp/matrix.txt".format(pu,protBtmp.pathTmp), "r").read()
            logging.info(sortie)

            a = Alignement(sortie, i)

            matRot = import_matrix("./tmp/matrix.txt")

            puTmp = Protein(pu)
            puTmp.rotate(matRot)
            #addition element par element de tuple
            (sumTmA, sumTmB) = map(lambda x,y: x+y,(sumTmA, sumTmB),sumTermeTM(puTmp, protBtmp, a))
            lastAtom, lastChain = write_pdb_align(outPath+"level_"+str(i+1)+".pdb", puTmp, lastAtom, lastChain, align = a.list_alignedA)
            protBtmp.remove_aligned(a.list_alignedB)
            protBtmp.write_pdbtmp()
            pcA += len(a.list_alignedA)
            pcB += len(a.list_alignedB)

            file_plot_align(a.list_alignedA, i, numPu, len(protB.seq), puTmp.numRes, protA.numRes, outPath)

        TMscoreA = tm_score(sumTmA, len(protA.seq))
        TMscoreB = tm_score(sumTmB, len(protB.seq))
        file_plot_score(i+1, TMscoreA, TMscoreB, outPath)
        logging.info("level = "+str(i))
        logging.info("TMscoreA = "+str(TMscoreA))
        logging.info("TMscoreB = "+str(TMscoreB))
        logging.info("Nombre alignés = "+str(pcA))
        logging.info("Pourcentage alignés A = "+str(pcA/len(protA.seq)))
        logging.info("Pourcentage alignés B = "+str(pcB/len(protB.seq)))
        i+=1

    os.system("./src/plotAlign.R "+str(i)+" "+str(len(protB.seq))+" "+outPath)
    os.system("./src/plotTMscore.R "+outPath)

if __name__ == '__main__':

    if(len(sys.argv) == 3):
        if(os.path.isfile(sys.argv[1])):
            if(os.path.isfile(sys.argv[2])):
                upalign(sys.argv[1], sys.argv[2])
                sleep(1)
                upalign(sys.argv[2], sys.argv[1])
            else:
                logging.error("erreur : "+sys.argv[2]+" n'existe pas !\n")
                exit("erreur : "+sys.argv[2]+" n'existe pas !\n")
        else:
            logging.error("erreur : "+sys.argv[1]+" n'existe pas !\n")
            exit("erreur : "+sys.argv[1]+" n'existe pas !\n")
    else:
        logging.error("erreur : nombre d'argument non valide.\n")
        exit("erreur : nombre d'argument non valide.\n")


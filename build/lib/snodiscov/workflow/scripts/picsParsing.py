#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd
from pathlib import Path
import os 
import sys
import re


## Recup Srtingtie ID to include in Parsing results file

## Recup class code for localisation pics and gffcompare Parsing results file
def allpicsParsing(allpics, gtfPics, gffcompParsingOut, allpicsParsingOut):
    sys.stdout=open(gffcompParsingOut, 'w')
    dicoClassCode=defaultdict(list)
    print("Id","Localisation","strand", sep='\t')
    with open(allpics, 'r') as file:
        for li in file:
            li=li.rstrip('\n').split()
            chr=li[0]; star=li[1]; end=li[2]; strand=li[3]; classcode=li[5]; id=chr+':'+star+'-'+end
            # snoDiscov v1.0: key on chr:start (was start only -> collisions
            # between loci starting at the same position on different chromosomes)
            dicoClassCode[chr+':'+star]=[strand,classcode]
            print(id,classcode, strand, sep='\t')

    sys.stdout.close()
    sys.stdout = sys.__stdout__  # was open('/dev/stdout','w'): it truncated redirected logs

    sys.stdout=open(allpicsParsingOut, 'w')
    print('Id','gene_id','Chr','Star','End','Length','strand','Class_code', sep=('\t'))
    with open(gtfPics, 'r') as allfile:
        for li in allfile:
            li=li.rstrip('\n').split()
            if li[2]=='gene':
                chr=li[0]; star=li[3]; end=li[4]; strgid=re.sub("\"|;","",li[9])
                id=str(chr+':'+star+'-'+end); length=int(end)-int(star)
                print(id,strgid,chr,star,end,length,dicoClassCode[chr+':'+star][0],dicoClassCode[chr+':'+star][1], sep='\t')
                

    sys.stdout.close()
    sys.stdout = sys.__stdout__  # was open('/dev/stdout','w'): it truncated redirected logs


####################################
#----  Parclip REsults Parsing ----#
####################################

def parclipParsing(data, parclipParsingOut, test):
    cdfilename=[
        "NOP58_parclip",
        "NOP56_parclip",
        "FBL_parclip"
    ]
    sys.stdout=open(parclipParsingOut, 'w')
    print('Id',test, sep='\t')
    with open(data, 'r') as parclipfile:
        for li in parclipfile:
            li=li.strip('\n').split()
            chr=re.sub("chr","",li[6]); star=li[7]; end=li[8]; score=float(li[4])
            Id=str(chr+':'+star+'-'+end)
            if test in cdfilename:
                print(Id,"cd", sep='\t')
            else:
                print(Id,"haca", sep='\t')


    sys.stdout.close()
    sys.stdout = sys.__stdout__  # was open('/dev/stdout','w'): it truncated redirected logs


##################################
#---- Eclip Results Parsing -----#
##################################

def eclipParsing(data, eclipParsingOut):
    sys.stdout=open(eclipParsingOut, 'w')
    print('Id','DKC1_eclip', sep='\t')
    with open(data, 'r') as eclipfile:
        for li in eclipfile:
            li=li.rstrip('\n').split()
            chr=re.sub("chr","", li[6]); star=li[7]; end=li[8]
            Id=str(chr+':'+star+'-'+end)
            print(Id,'haca', sep='\t')

    sys.stdout.close()
    sys.stdout = sys.__stdout__  # was open('/dev/stdout','w'): it truncated redirected logs


##################################################
#------- Depletion Results DESeq2 Parsed --------#
##################################################

def deseq2Parsing(deseq2tab, snoGtf, tab_out, testType):
    dkc1Dep_dico=defaultdict(list) #Dico for retrieve gene id, log2foldChange and significative pvalue (<=.05)
    with open(deseq2tab, 'r') as dkc1file:
        header=dkc1file.readline()
        for li in dkc1file:
            li=li.rstrip('\n').split()
            Gene_id=re.sub("\"","",li[0]); lfc=li[2]; pvalue=li[5]
            if testType.startswith("dkc1_depletion"):
                if float(pvalue) <= 0.05 and float(lfc)<=-0.5:
                    res="haca" 
                    dkc1Dep_dico[Gene_id]=[res, lfc, pvalue]
                else:
                    res='negative'
                    dkc1Dep_dico[Gene_id]=[res, lfc, pvalue]
            else:
                if float(pvalue) <= 0.05 and float(lfc)<=-0.5:
                    res="cd" 
                    dkc1Dep_dico[Gene_id]=[res, lfc, pvalue]
                else:
                    res='negative'
                    dkc1Dep_dico[Gene_id]=[res, lfc, pvalue]

    sys.stdout=open(tab_out, 'w')
    print('Id','boxType', 'l2fc','pvalue', sep='\t')
    with open(snoGtf, 'r') as annotfile:
        for li in annotfile:
            li=li.rstrip('\n').split()
            if li[2]=='transcript':
                gene_id=re.sub("\"|;","",li[9]); star=li[3]; id=li[0]+':'+star+'-'+li[4]
                if gene_id in dkc1Dep_dico:
                #dkc1Annot_dico[gene_id]=id
                    print(id,dkc1Dep_dico[gene_id][0],dkc1Dep_dico[gene_id][1], dkc1Dep_dico[gene_id][2],sep='\t')

    sys.stdout.close()
    sys.stdout = sys.__stdout__  # was open('/dev/stdout','w'): it truncated redirected logs

#########################################################################################

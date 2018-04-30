#!/usr/bin/env python

# Main libray conaining all the functions to calculate features.
#
# Written by Antoine LU @gmail.com in December 2017 @UMR7216 Paris Diderot.

# -*- coding: utf-8 *-*

#lots de fonction a inclure pour l'execution du script test #TODO not the test but this should be the main script.

import math
import re
import pandas as pd

def getKozak(df, k, j):
    '''
    Retrieve the context of the Kozak of each sequence in a dataframe
    IMPORTANT: if no cdna start (codon START position) is given, doesn't
    work --> amelioration, chercher A/GccATGG
    Input:
    df: Pandas dataframe, read from a csv file
    k: default, k = 10: number of nucleotide before START-1 (A)
        and after START+1 (G) of the START codon (ATG)
        >> k---START---k is the Kozak sequence
    j: default, j = 20: number of nucleotide to select around
        the kozak sequence
    Output:
    kozak: which is the sequence of the context of the Kozak
    '''
    kozak = []
    for indice in range(len(df['tmp_seq'])):
        seq = df.iloc[indice, df.columns.get_loc('tmp_seq')]
        #Nan supprimes prealablement du dataset df
        if math.isnan(df.iloc[indice, df.columns.get_loc('cdnastart')]):
            print df.at[indice, 'cdnastart']
            break
        startM1 = int(df.iloc[indice, df.columns.get_loc('cdnastart')]) - 1 # position start-1
        startP1 = int(df.iloc[indice, df.columns.get_loc('cdnastart')]) + 1 # position start+1
        #kozak = seq[(startM1 - k):(startP1 + k)]
        contextM1 = (startM1 - k) - j
        contextP1 = (startP1 + k) + j
        contextKozak = seq[((startM1 - k) - j):((startP1 + k) + j)]
        kozak.append(contextKozak)
    return(kozak)

def get_uORF(df):			#avec uORFs chevauchantes
	uORFs_TOT=[]
	for indice in range(len(df)):
		ind_uORF=[]
		seq=df.iloc[indice,df.columns.get_loc('tmp_seq')]
		p5UTR=seq[0:int(df.iloc[indice,df.columns.get_loc("cdnastart")])-1]
		#print p5UTR,len(p5UTR),len(seq),"seq-current\n"
		reg=re.compile('ATG')
		for m_ATG in reg.finditer(p5UTR):
			if m_ATG.start() % 3 == 0:
				ATG=m_ATG.start()
				rex=re.compile('ATG([ATGC]{3}){1,}T(AG|AA|GA)')
				research_zone=p5UTR[ATG:]
				while(len(research_zone)>=9):
					m_rex=rex.search(research_zone)
					if m_rex is not None:
						match=m_rex.group()
						start=m_rex.start()
						end=m_rex.end()
						uORF=(ATG+start,ATG+end)
						ind_uORF.append(uORF)
						#print p5UTR[ATG+start:ATG+end],len(p5UTR[ATG+start:ATG+end])
						research_zone=match[start:end - 3]
					else: break
		uORFs_TOT.append(ind_uORF)
	return(uORFs_TOT)

def get_dORF(df):
    dORFs = []
    for indice in range(len(df)):
        ind_dORF = []
       	seq=df.iloc[indice,df.columns.get_loc('tmp_seq')]
        p3UTR= seq[int(df.iloc[indice,df.columns.get_loc('cdnaend')]):]
        reg = re.compile('ATG')
        for m_ATG in reg.finditer(p3UTR):
            if m_ATG.start() % 3 == 0:
                ATG = m_ATG.start()
                rex = re.compile('ATG([ATGC]{3}){1,}T(AG|AA|GA)')
                research_zone=p3UTR[ATG:]
                while(len(research_zone)>=9):
                    m_rex = rex.search(research_zone)
                    if m_rex is not None:
                        match = m_rex.group()
                        start = m_rex.start()
                        end = m_rex.end()
                        dORF = (ATG + start, ATG + end)
                        ind_dORF.append(dORF)
                        research_zone=match[start:end-3]
                    else:
                        break
        dORFs.append(ind_dORF)
    return(dORFs)

def writeP5utr_fa(df):
    '''
    Write every 5'UTR sequence in a fasta file named "p3utr.fasta".
    '''
    with open("../data/p5utrDEMO.fasta", "w+") as p5utr_fasta:
        for indice in range(len(df)):
        	seq = df.iloc[indice,df.columns.get_loc('tmp_seq')]
        	p5UTR = seq[0:int(df.iloc[indice,df.columns.get_loc("cdnastart")])-1]
        	if len(p5UTR)>0:
				p5utr_fasta.write("{}\n".format(">" + str(indice)))
				p5utr_fasta.write("{}\n".format(p5UTR))

def writeP3utr_fa(df):
    '''
    Write every 3'UTR sequence in a fasta file named "p3utr.fasta".
    '''
    with open("../data/p3utrDEMO.fasta", "w+") as p3utr_fasta:
        for indice in range(len(df)):
        	seq = df.iloc[indice,df.columns.get_loc('tmp_seq')]
        	p3UTR = seq[int(df.iloc[indice,df.columns.get_loc('cdnaend')]):]
        	if len(p3UTR) > 0:
				p3utr_fasta.write("{}\n".format(">" + str(indice)))
				p3utr_fasta.write("{}\n".format(p3UTR))

def getFoldingEnergy(filename, df):
    with open(str(filename), "r") as rnafoldfile:
        tot = rnafoldfile.readlines()
    foldrex = re.compile('(-[0-9]+\.[0-9]+|\s0\.0)')
    foldinf = re.compile('>[0-9]+')
    strtot = ' '.join(tot) # convertion liste en string
    mfe = foldrex.findall(strtot)
    indice = foldinf.findall(strtot)
    real_indice = [None] * len(df) # remplacer 4 par len(cdd)
    for i in range(len(indice)):
        new_id = int(indice[i][1:])
        real_indice[new_id] = float(mfe[i])
    return(real_indice)

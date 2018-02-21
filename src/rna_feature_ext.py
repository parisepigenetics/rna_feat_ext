#! /usr/bin/env python
# -*- coding: utf-8 -*-

# MODULES
import numpy as np
import pandas as pd
import math
import re
import os
import subprocess # doesn't work


# FUCTIONS
# 0.1) Function to get 5' UTR sequences
def writeP5utr_fa(df):
    '''
    Write every 5'UTR sequence in a fasta file named "p3utr.fasta".
    '''
    with open("../data/p5utrT.fasta", "w+") as p5utr_fasta:
        for indice in range(len(df)):
            if math.isnan(df.get_value(indice, 'cdnastart')):
                # 3'UTR start position = length cDNA - length 3'UTR
                startP5 = df.get_value(indice, 'p5start')
                endP5 = df.get_value(indice, 'p5end')
                lenP5 = endP5-startP5+1
                #print(">"+str(indice))
                #print(metatable.get_value(indice, 'tmp_seq')[posP3starTrue:])
                if len(df.get_value(indice, 'tmp_seq')[lenP5 + 1:]) > 0:
                    p5utr_fasta.write("{}\n".format(">" + str(indice)))
                    p5utr_fasta.write("{}\n".format(df.get_value(indice, 'tmp_seq')[lenP5 + 1:]))
            else:
                if len(df.get_value(indice, 'tmp_seq')[0:int(df.get_value(indice, 'cdnastart')) - 1]) > 0:
                    p5utr_fasta.write("{}\n".format(">" + str(indice)))
                    p5utr_fasta.write("{}\n".format(df.get_value(indice, 'tmp_seq')[0:int(df.get_value(indice, 'cdnastart')) - 1]))

# 0.2) Function to get 3' UTR sequences
def writeP3utr_fa(df):
    '''
    Write every 3'UTR sequence in a fasta file named "p3utr.fasta".
    '''
    with open("../data/p3utrT.fasta", "w+") as p3utr_fasta:
        for indice in range(len(df)):
            if math.isnan(df.get_value(indice, 'cdnaend')):
                # 3'UTR start position = length cDNA - length 3'UTR
                startP3 = df.get_value(indice, 'p3start')
                endP3 = df.get_value(indice, 'p3end')
                lenP3 = endP3-startP3+1
                posP3starTrue = int(len(df.get_value(indice, 'tmp_seq')) - lenP3)
                #print(">"+str(indice))
                #print(metatable.get_value(indice, 'tmp_seq')[posP3starTrue:])
                if len(df.get_value(indice, 'tmp_seq')[posP3starTrue:]) > 0:
                    p3utr_fasta.write("{}\n".format(">" + str(indice)))
                    p3utr_fasta.write("{}\n".format(df.get_value(indice, 'tmp_seq')[posP3starTrue:]))
            else:
                # ELSE: 3'UTR starts after the cdna END
                #print(">"+str(indice))
                #print(metatable.get_value(i, 'tmp_seq')[int(metatable.get_value(indice, 'cdnaend')):])
                #q = (metatable.get_value(i, 'tmp_seq')[int(metatable.get_value(indice, 'cdnaend')):])
                if len(df.get_value(indice, 'tmp_seq')[int(df.get_value(indice, 'cdnaend')):]) > 0:
                    p3utr_fasta.write("{}\n".format(">" + str(indice)))
                    p3utr_fasta.write("{}\n".format(df.get_value(indice, 'tmp_seq')[int(df.get_value(indice, 'cdnaend')):]))

# 1) RNA Binding Protein motif (RBP)
# See script named motif_scan_projet_long.py

# 4 and 5) Folding Energy
def getFoldingEnergy(filename, df):
    with open(str(filename), "r") as rnafoldfile:
        tot = rnafoldfile.readlines()
        
    foldrex = re.compile('(-[0-9]+\.[0-9]+|\s0\.0)')
    foldinf = re.compile('>[0-9]+')

    strtot = ' '.join(tot) # convertion liste en string

    mfe = foldrex.findall(strtot)
    indice = foldinf.findall(strtot)

    # real_indice: table à convertir en pd.DataFrame par la suite
    # contiendra les MFE pour chaque séquence
    # s'il n'y avait pas de séquence: la cellule sera 'None'
    real_indice = [None] * len(df) # remplacer 4 par len(metatable)
    for i in range(len(indice)):
        new_id = int(indice[i][1:])
        real_indice[new_id] = float(mfe[i])
    
    return(real_indice)

# 6) Downstrea ORF
def get_dORF(df):
    '''
    This function calculates for each gene the coordinates of the downstream open
    reading frames in the 3' UTR.
    Input: Pandas dataframe, read from a csv file
    (HARD CODED) minimal dORF size = 9 (size found in bibliography).
    Output: return a list of list of tuple (one list of tuple per gene).
    Each tuple (1, more or no) dORF containing the coords of the dORF in the 3'UTR.
    If no dORF are found, the tuple is empty, in concequence, so is the list (of tuple).
    '''
    
    # grande liste à retourner et à transformer en colone 6 dans le tableau final
    # TODO: identification phase de lecture (idem pour uORF)
    dORFs = []
    for indice in range(len(df)):
        ind_dORF = []
        #print("Gene num ",indice)
        q = df.get_value(indice, 'tmp_seq')
        if math.isnan(df.get_value(indice, 'cdnaend')):
            # 3'UTR start position = length cDNA - length 3'UTR
            startP3 = df.get_value(indice, 'p3start')
            endP3 = df.get_value(indice, 'p3end')
            lenP3 = endP3-startP3+1
            posP3starTrue = int(len(q) - lenP3)
            p3utr = q[posP3starTrue:]
        else:
            # ELSE: 3'UTR starts after the cdna END
            p3utr = q[int(df.get_value(indice, 'cdnaend')):]
        # Search the position of the ATGs in the 3'UTR
        # IF the position of the 'A' of ATG is a multiple of 3, keep it
        # IF NOT: the fetched "ATG" might not be a real ATG
        reg = re.compile('ATG')
        for m_ATG in reg.finditer(p3utr):
            # IF Multiple of 3
            if m_ATG.start() % 3 == 0:
                ATG = m_ATG.start()
                #print(ATG)
                rex = re.compile('ATG([ATGC]{3}){1,}T(AG|AA|GA)')
                zone_de_recherche = p3utr[ATG: ]
                # uORF a une taille >= 9 nucleotides
                while(len(zone_de_recherche) >= 9):
                    m_rex = rex.search(zone_de_recherche)
                    if m_rex is not None:
                        match = m_rex.group()
                        start = m_rex.start()
                        end = m_rex.end()
                        #print match, ATG + start, ATG + end
                        dORF = (ATG + start, ATG + end)
                        ind_dORF.append(dORF)
                        zone_de_recherche = match[start:end - 3]
                    else: break
        # append an empty list when no uORF found in 5'UTR
        dORFs.append(ind_dORF)
    return(dORFs)

# 7) Upstream ORF
def get_uORF(df):
    '''
    This function calculates for each gene the coordinates of the upstream open
    reading frames in the 5' UTR.
    Input: Pandas dataframe, read from a csv file
    (HARD CODED) minimal uORF size = 9 (size found in bibliography).
    Output: return a list of list of tuple (one list of tuple per gene).
    Each tuple (1, more or no) uORF containing the coords of the uORF in the 5'UTR.
    If no uORF are found, the tuple is empty, in concequence, so is the list (of tuple).
    '''
    
    # grande liste à retourner et à transformer en colone 6 dans le tableau final
    uORFs = []
    for indice in range(len(df)):
        ind_uORF = []
        #print("Gene num ",indice)
        seq = df.get_value(indice, 'tmp_seq')
        # 5'UTR: starts at the first nucleotide of the tmp_seq sequence (0)
        # ends before the cdnastart (that's why --> -1)
        # >> Only uORF are fetched, not overlapping UTR (oUTR)
        p5utr = seq[0:int(df.get_value(indice, 'cdnastart')) - 1]
        # Search the position of the ATGs in the 5'UTR
        # IF the position of the 'A' of ATG is a multiple of 3, keep it
        # IF NOT: the fetched "ATG" might not be a real ATG
        reg = re.compile('ATG')
        for m_ATG in reg.finditer(p5utr):
            # IF Multiple of 3
            if m_ATG.start() % 3 == 0:
                ATG = m_ATG.start()
                #print(ATG)
                rex = re.compile('ATG([ATGC]{3}){1,}T(AG|AA|GA)')
                zone_de_recherche = p5utr[ATG: ]
                # uORF a une taille >= 9 nucleotides
                while(len(zone_de_recherche) >= 9):
                    m_rex = rex.search(zone_de_recherche)
                    if m_rex is not None:
                        match = m_rex.group()
                        start = m_rex.start()
                        end = m_rex.end()
                        #print match, ATG + start, ATG + end
                        uORF = (ATG + start, ATG + end)
                        ind_uORF.append(uORF)
                        zone_de_recherche = match[start:end - 3]
                    else: break
        # append an empty list when no uORF found in 5'UTR
        uORFs.append(ind_uORF)
    return(uORFs)

# 9) Kozak sequence context
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
        seq = df.get_value(indice, 'tmp_seq')
        if math.isnan(df.get_value(indice, 'cdnastart')):
            next
        else:
            startM1 = int(df.get_value(indice, 'cdnastart')) - 1 # position start-1
            startP1 = int(df.get_value(indice, 'cdnastart')) + 1 # position start+1
            #kozak = seq[(startM1 - k):(startP1 + k)]
            contextM1 = (startM1 - k) - j
            contextP1 = (startP1 + k) + j
            contextKozak = seq[((startM1 - k) - j):((startP1 + k) + j)]
            kozak.append(contextKozak)
    return(kozak)

# MAIN
if __name__ == '__main__':
    # DATA: Prerequis
    metatable = pd.read_csv('02012018/tableMetagene.csv')

    # Preparation of the concatenation of each feature:
    # select only gene id and transcript id of the metatable dataframe
    data = metatable.loc[:,['geneid','trid']]

	# 0) UTRs sequence retrieving (write in a file)
	writeP5utr_fa(metatable)
	writeP3utr_fa(metatable)

    # 1) RNA Binding Protein motif (RBP)
    
    # Either: in the terminal: (if possible use subprocess)
    # python motif_scan/motif_scan_projet_long.py -c 8 ../data/RBP/toread/xab > testxab.csv
    
    # OR in a python environement:
    # path = "data/RBP/toread/"
    # files = []
    # for i in os.listdir(path):
    #     if os.path.isfile(os.path.join(path,i)) and 'x' in i:
    #         files.append(i)
    # for fichier in files[20:30]:
    #     with open('results/RBP_res/RBP_out_'+fichier+'.csv', 'w') as f:
    #         subprocess.call(['python',
    #                           'src/motif_scan/motif_scan_projet_long.py',
    #                           '-c 8',
    #                           'data/RBP/toread/'+fichier],
    #                         stdout = f)
    
    rbp = pd.read_csv('../results/RBP/RBP_outlist.csv')
    data['rbp'] = rbp

	# 2) Length of 3' UTR
	p3len = metatable['p3end']-metatable['p3start']
	pd.DataFrame(p3len, columns = ['utr3_len'])
    data['p3len'] = p3len

	# 4) Folding Energy of 3' UTR
    #TODO: useSubprocess
	os.system("RNAfold < ../data/p3utrT.fasta --noPS > ../results/06012018/p3utr_rnafold.txt")

	p3mfe = getFoldingEnergy("06012018/p3utr_rnafold.txt", metatable)
	p3Fold = pd.DataFrame(p3mfe, columns = ['p3_folding'])
    data['p3mfe'] = p3mfe

	# 5) Folding Energy of 5' UTR
	os.system("RNAfold < ../data/p5utrT.fasta --noPS > ../results/06012018/p5utr_rnafold.txt")

	p5mfe = getFoldingEnergy("06012018/p5utr_rnafold.txt", metatable)
	p5Fold = pd.DataFrame(p5mfe, columns = ['p5_folding'])
	data['p5mfe'] = p5mfe

    # 6) Downstream ORF
    dORFs = get_dORF(metatable)

    empty = [None] * len(metatable)
    # create empty DF of length(metatable)
    dORF = pd.DataFrame(empty, columns=['dORF'])
    # for each cell (and same indice), put uORF value in empty df
    for ind in range(len(uORFs)):
        dORF.set_value(ind, 'dORF', dORFs[ind])

    data['dORF'] = dORF

    # 7) Upstream ORF
    uORFs = get_uORF(metatable)

    empty = [None] * len(metatable)
    uORF = pd.DataFrame(empty, columns=['uORF'])
    for ind in range(len(uORFs)):
        uORF.set_value(ind, 'uORF', uORFs[ind])

    data['uORF'] = uORF

    # 9) Kozak sequence context
    kozak = getKozak(metatable, 10, 20)
    contKozak = pd.DataFrame(kozak, columns = ['contextKozak'])
    data['kozak_context'] = contKozak

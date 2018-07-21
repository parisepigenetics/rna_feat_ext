# -*- coding: utf-8 *-*

#lots de fonction à inclure pour l'execution du script test

import math
import re
import pandas as pd
import subprocess

# PART 1: From IDs to FASTA
def get_utr5MAX(cdna_feat_row):
	#print("getting UtrMax...")
	#take the cdna_feat-row with multiples utrs

	utr5s_start = cdna_feat_row["5' UTR start"].values[0].split(";")
	utr5s_end = cdna_feat_row["5' UTR end"].values[0].split(";")

	size_liste = []
	for i in range(len(utr5s_start)):
		size = int(utr5s_end[i])-int(utr5s_start[i])
		size_liste.append(size)
	indice_max = size_liste.index(max(size_liste))
	max_utr5_start = int(utr5s_start[indice_max])
	max_utr5_end = int(utr5s_end[indice_max])

	return([max_utr5_start,max_utr5_end])

def get_utr3MAX(cdna_feat_row):
    #print("getting UtrMax...")
    #take the cdna_feat-row with multiples utrs

    utr3s_start = cdna_feat_row["3' UTR start"].values[0].split(";")
    utr3s_end = cdna_feat_row["3' UTR end"].values[0].split(";")

    size_liste = []
    for i in range(len(utr3s_start)):
        size = int(utr3s_end[i])-int(utr3s_start[i])
        size_liste.append(size)
    indice_max = size_liste.index(max(size_liste))
    max_utr3_start = int(utr3s_start[indice_max])
    max_utr3_end = int(utr3s_end[indice_max])

    return([max_utr3_start,max_utr3_end])

def get_cDNAstartMIN(cdna_feat_row):
    #print("getting cDNA_Start_Min...")
    #take the cdna_feat-row with multiples utrs

    cDNA_start = cdna_feat_row["cDNA coding start"].values[0].split(";")
    min_cDNA_start = min(map(int, cDNA_start))

    return(min_cDNA_start)

def get_cDNAendMAX(cdna_feat_row):
    #print("getting cDNA_End_Max...")
    #take the cdna_feat-row with multiples utrs

    cDNA_end = cdna_feat_row["cDNA coding end"].values[0].split(";")
    max_cDNA_end = max(map(int, cDNA_end))

    return(max_cDNA_end)

def txt2fasta(cdna_feat_table, fastaOut):
    with open(fastaOut, "w+") as ff:
        for i in range(cdna_feat_table.shape[0]):
            ligne = pd.DataFrame(cdna_feat_table.loc[i,:]).transpose()
            ff.write(">GeneID:{}|TranscriptID:{}|GeneName:{}|5P_UTR_end:{}|5P_UTR_start:{}|3P_UTR_end:{}|3P_UTR_end:{}|cDNAstart:{}|cDNAend:{}\n".format(
                ligne["Gene stable ID"].values[0], ligne["Transcript stable ID"].values[0],
                ligne["Gene name"].values[0], ligne["5' UTR end"].values[0], ligne["5' UTR start"].values[0],
                ligne["3' UTR end"].values[0], ligne["3' UTR start"].values[0],
                ligne["cDNA coding start"].values[0], ligne["cDNA coding end"].values[0]
            ))
            ff.write("{}\n".format(ligne["cDNA sequences"].values[0]))
        print("txt to Fasta conversion done!")



# PART 2: Computation of the feature table
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
    print("extracting Kozak_Context ...")
    kozak = []
    kozakseq = []
    for indice in range(len(df['cDNA sequences'])):
        seq = df.iloc[indice, df.columns.get_loc('cDNA sequences')]
        #Nan supprimes prealablement du dataset df

        startM1 = int(df.iloc[indice, df.columns.get_loc('cDNA coding start')]) - 1 # position start-1
        startP1 = int(df.iloc[indice, df.columns.get_loc('cDNA coding start')]) + 1 # position start+1
        #kozak = seq[(startM1 - k):(startP1 + k)]
        contextM1 = (startM1 - k) - j
        contextP1 = (startP1 + k) + j
        #Cas ou cadre du context_kozak dépasse l'indice de début de sequence
        contextM1=abs(contextM1)
        contextKozak = seq[contextM1:contextP1]
        #print len(contextKozak),contextKozak,contextM1,contextP1, startM1
        kozak.append(contextKozak)
        kozakseq.append(seq[(startM1-k):(startP1+k)])

    print("extracting Kozak_Context ... Done !")
    return([kozak,kozakseq])

def get_uORF(df):			#avec uORFs chevauchantes
    print("getting uORFs_sequence ...")
    uORFs_TOT = []
    for indice in range(len(df)):
    	ind_uORF=[]
    	seq=df.iloc[indice,df.columns.get_loc('cDNA sequences')]
    	p5UTR=seq[0:int(df.iloc[indice,df.columns.get_loc('cDNA coding start')])-1]
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
    print("getting uORFs_sequence ... Done !")
    return(uORFs_TOT)

def get_dORF(df):
    print("getting dORFs_sequence ...")
    dORFs = []
    for indice in range(len(df)):
        ind_dORF = []
       	seq=df.iloc[indice,df.columns.get_loc('cDNA sequences')]
        p3UTR= seq[int(df.iloc[indice,df.columns.get_loc('cDNA coding end')]):]
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
                    else: break
        dORFs.append(ind_dORF)
    print("getting dORFs_sequence ... Done !")
    return(dORFs)

def writeP5utr_fa(df):
    '''
    Write every 5'UTR sequence in a fasta file named "p3utr.fasta".
    '''
    print("Creating p5UTR_seq ...")
    p5len = []
    with open("p5utrDEMO.fasta", "w+") as p5utr_fasta:
        for indice in range(len(df)):
        	seq=df.iloc[indice,df.columns.get_loc('cDNA sequences')]
        	p5UTR=seq[0:int(df.iloc[indice,df.columns.get_loc("cDNA coding start")])-1]
        	if len(p5UTR)>0:
				p5utr_fasta.write("{}\n".format(">" + str(indice)))
				p5utr_fasta.write("{}\n".format(p5UTR))
                p5len.append(len(p5UTR))
    print("Creating p5UTR_seq ... Done !")
    return(p5len)

def writeP3utr_fa(df):
    '''
    Write every 3'UTR sequence in a fasta file named "p3utr.fasta".
    '''
    print("Creating p3UTR_seq ...")
    p3len = []
    with open("p3utrDEMO.fasta", "w+") as p3utr_fasta:
        for indice in range(len(df)):
        	seq=df.iloc[indice,df.columns.get_loc('cDNA sequences')]
        	p3UTR= seq[int(df.iloc[indice,df.columns.get_loc('cDNA coding end')]):]
        	if len(p3UTR)>0:
				p3utr_fasta.write("{}\n".format(">" + str(indice)))
				p3utr_fasta.write("{}\n".format(p3UTR))
                p3len.append(len(p3UTR))
    print("Creating p3UTR_seq ...Done")
    return(p3len)

def RNAfold_calcul(inputfile, output_name_file):
    print("RNAfold_Calcul ...")
    output=open(output_name_file, "w+")
    subprocess.call("RNAfold --noPS --jobs", stdin = inputfile, stdout = output, shell = True)
    output.close()
    print("RNAfold_Calcul ... Done !")
    return(output_name_file)

def getFoldingEnergy(filename, df):
    with open(str(filename), "r") as rnafoldfile:
        tot = rnafoldfile.readlines()
    #TODO use the delimiter of the RNAFold output file to extract the folding energy (the split function). Remove the REGEXPs.
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

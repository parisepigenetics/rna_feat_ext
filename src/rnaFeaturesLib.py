#!/usr/bin/env python
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
    with open(fastaOut + ".fasta", "w+") as ff:
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

def RNAfold_calcul(utr_fasta, out_mfe):
    #print("RNAfold_Calcul ...")
    with open(utr_fasta, "r") as inputfile, open(out_mfe + ".mfe", "w+") as output:
        subprocess.call("RNAfold --noPS --jobs", stdin = inputfile, stdout = output, shell = True)
    print("RNAfold_Calcul ... Done !")
    return(output)

def getFoldingEnergy(input_mfe):
    with open(str(input_mfe), "r") as rnafoldfile:
        tot = rnafoldfile.readlines()
    #TODO use the delimiter of the RNAFold output file to extract the folding energy (the split function). Remove the REGEXPs.
    foldrex = re.compile('(-[0-9]+\.[0-9]+|\s0\.0)')
    foldinf = re.compile('>[0-9]+')
    strtot = ' '.join(tot) # convertion liste en string
    mfe = foldrex.findall(strtot)
    return(float(mfe[0]))

class FastaReader(object):
    # Read Fasta input sequence, store each seq in a generator
    def __init__(self, fastaFile):
        """Constructor"""
        self.fastaFile = fastaFile
    #
    def readSeqs(self):
        mySeq = []
        currentSeq = ''
        with open(self.fastaFile, 'r') as f:
            for line in f:
                if line == '\n':
                    next
                m = re.match('^>.*',line)
                if m:
                    if currentSeq:
                        yield Seq(currentSeq, ''.join(mySeq))
                    #
                    #currentSeq = re.split('[\s\|]+',m.group(1))[0]
                    currentSeq = m.group(0)
                    mySeq = []
                else:
                    mySeq.append(line.replace('\n',''))
            yield Seq(currentSeq, ''.join(mySeq))


class Seq(object):
    # Sequence define by its name and its bases sequence
    def __init__(self, name, bases):
        """Constructor"""
        self.name = name
        self.bases = bases

class Job(object):
    def __init__(self, SeqObj):
        self.SeqObj = SeqObj
        
    def getAttributes(self):
        # TODO get the attibute names dirrectly from the fasta header. Minimise hard coding.
        attributes = self.SeqObj.name.split(">")[1].split("|")
        self.geneID = attributes[0].split(":")[1]
        self.transcriptID = attributes[1].split(":")[1]
        self.geneName = attributes[2].split(":")[1]
        self.UTR5PEnd = attributes[3].split(":")[1]
        self.UTR5PStart = attributes[4].split(":")[1]
        self.UTR3PEnd = attributes[5].split(":")[1]
        self.UTR3PStart = attributes[6].split(":")[1]
        self.cDNAStart = attributes[7].split(":")[1]
        self.cDNAEnd = attributes[8].split(":")[1]
        self.bases = self.SeqObj.bases
    
    def getKozak(self, s, c):
        """ Get both Kozak sequence and context:
        If ATG is located near the 5'UTR, it is most likely probable that
        either seq to fetch will not be retrieved as seld.base[-5:10] returns blank.
        
        In this case, we test whether the value left to ':' is negative or not
        If it is < 0, we simply take the seq from 0 as : self.bases[0:self.cDNAStart) + 2) + s]
        """
        # Kozak sequence
        if int(self.cDNAStart) - 1 - s < 0:
            kozakSeq = self.bases[0:(int(self.cDNAStart) + 2) + s]
        else:
            kozakSeq = self.bases[(int(self.cDNAStart) - 1 - s):(int(self.cDNAStart) + 2)+ s]
        
        # Kozak context
        if int(self.cDNAStart) -1 - s - c < 0:
            kozakContext = self.bases[0:(int(self.cDNAStart) + 2)+ s + c]
        else:
            kozakContext = self.bases[((int(self.cDNAStart) - 1 - s) - c):(int(self.cDNAStart) + 2) + s + c]
        return(kozakSeq, kozakContext)
    
    # NO need to create files rnafold can read/write to stdin and stdout.
    def write3PUTR(self, out_3putr):
        """Fetch 3PUTR sequence, write in external fasta file and return its length"""
        with open(out_3putr + ".fasta", "w") as UTR3P:
            UTR3P.write(">{}_3PUTR\n{}".format(self.geneID,self.bases[int(self.cDNAEnd):]))
        #print("Creating p3UTR_seq ...Done")
        return(len(self.bases[int(self.cDNAEnd):]))
    
    def write5PUTR(self, out_5putr):
        """Fetch 5PUTR sequence, write in external fasta file and return its length"""
        with open(out_5putr + ".fasta", "w") as UTR5P:
            UTR5P.write(">{}_5PUTR\n{}".format(self.geneID,self.bases[0:int(self.cDNAStart)-1]))
        #print("Creating p5UTR_seq ...Done")
        return(len(self.bases[0:int(self.cDNAStart)-1]))
        
class FeatExtract(object):
    
    def __init__(self, fastaFile):
        self.fastaFile = fastaFile
    
    def getFeatures(self):
        ff = FastaReader(self.fastaFile)
        for seq in ff.readSeqs():
            # yield
            job = Job(seq)
            # for job in toto
            job.getAttributes()
            seqKozak, contKozak = job.getKozak(10,20)
            #print "{}\n{}".format(seqKozak, contKozak)
            job_3UTRlen = job.write3PUTR("job3utr")
            #print job_3UTRlen
            job_5UTRlen = job.write5PUTR("job5utr")
            #print job_5UTRlen
            # 3 UTR
            filin = "job3utr.fasta"
            filout = "job3utr"
            RNAfold_calcul(filin, filout)
            p3mfe = getFoldingEnergy(filout + ".mfe")
            engBase_3utr = p3mfe / job_3UTRlen
            #print engBase_3utr
            # 5 UTR
            filin = "job5utr.fasta"
            filout = "job5utr"
            RNAfold_calcul(filin, filout)
            p5mfe = getFoldingEnergy(filout + ".mfe")
            engBase_5utr = p5mfe / job_5UTRlen
            # dico
            features = {
               'ensembl_gene_id':job.geneID,
                'ensembl_transcript_id':job.transcriptID,
                '3PLen':job_3UTRlen,
                '3PMfe':p3mfe,
                '5PLen':job_5UTRlen,
                '5PMfe':p5mfe,
                '5UTR_mfe_Base':engBase_5utr,
                '3UTR_mfe_Base':engBase_3utr,
                'Kozak_Context':contKozak,
                'Kozak_sequence':seqKozak
            }
            yield features
            
    def dicos2table(self):
        """When treating the first yield dictionary, create the final table df_feat
        Each yield dictonaries is transformed into an entry of the table and concatenated to the final table"""
        i = 0
        for yieldFeature in self.getFeatures():
            if i == 0:
                df_feat = pd.DataFrame(columns=yieldFeature.keys())
                i = 1
            df  = pd.DataFrame([yieldFeature],columns=yieldFeature.keys())
            df_feat = pd.concat([df_feat, df], axis=0).reset_index(drop=True)

        return(df_feat)

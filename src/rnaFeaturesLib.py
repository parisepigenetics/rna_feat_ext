#!/usr/bin/env python
# -*- coding: utf-8 *-*

"""Python module for RNA fetures extraction from ENSEMBL derived fasta files.
"""

__version__ = "0.1a01"

import os
import math
import re
import subprocess
import pandas as pd
import Bio
from Bio import SeqIO
from Bio.SeqUtils import GC


class ENSEMBLSeqs(object):
    """Class to represent RNA related sequence features from ENSEMBLself.

    Needs a Bio.Seq.Record generator object to initialise."""

    def __init__(self, bioSeqRecsGen):
        self.gen = bioSeqRecsGen
        self.bioSeqRecs = self.getbioSeqRecs()

    def getbioSeqRecs(self):
        """Expands the Bio.Seq.Rec generator to a list of Bio.Record objects.

        Also puts some SeqIO.record member variables in place, namely the id, the gene name the description and the features."""
        recs = []
        for rec in self.gen:
            # Extract the gene name.
            rec.name = rec.description.split("|")[2].split(":")[1]
            descr = rec.description.split("|")[-1]
            # Extract all the features as a dictionary (apart from the last one which was the description.).
            feat = dict(item.split(":") for item in rec.description.split("|")[1:-1])
            rec.features = feat
            #TODO for the moment all the features are stored as a dictionary and not as proper SeqFeature objects. (perhaps we can stick with that and there is no need to change it.)
            rec.description = descr
            recs.append(rec)
        return recs


class FeaturesExtract(object):
    """Claas to extract features."""

    def __init__(self, bioSeqRecs):
        """Initialise with a list of SeqIO records."""
        self.bioSeqRecs = bioSeqRecs

    def collectFeatures(self):
        """Collect the features that do not need external computations."""
        # Check if files exist and delete them
        try:
            os.remove("3pUTRs.fasta")
        except OSError:
            pass
        utr3p_f = "3pUTRs.fasta"
        try:
            os.remove("5pUTRs.fasta")
        except OSError:
            pass
        utr5p_f = "5pUTRs.fasta"
        for rec in self.bioSeqRecs:
            # Get Koxaks
            seqKozak, contKozak = self.getKozak(rec)
            # Fetch UTRs, lenghts and GCs
            utr3len, utr3gc = self.get3PUTR(rec, utr3p_f)
            utr5len, utr5gc = self.get5PUTR(rec, utr5p_f)
            # Length of coding region.
            codeLen = int(rec.features["cDNAend"]) - int(rec.features["cDNAstart"])
            # dico
            features = {'ensembl_gene_id':rec.features["GeneID"],
                        'ensembl_transcript_id':rec.id,
                        'gene_name':rec.name,
                        'coding_len':codeLen,
                        '3pUTR_len':utr3len,
                        '5pUTR_len':utr5len,
                        '5pUTR_GC':"{0:.2f}".format(utr5gc),
                        '3pUTR_GC':"{0:.2f}".format(utr3gc),
                        'Kozak_Context':contKozak,
                        'Kozak_sequence':seqKozak}
            yield features

    def getKozak(self, rec, s = 10, c = 20):
        """Extract both Kozak sequence and context from a SeqIO record.
        (s and c define the extremeties of a Koxak sequence and are chosen by convention.)
        return a tuple of Kozak sequence and Kozak context.

        If ATG is located near the 5'UTR start, it is most likely that either seqs will not be retrieved as seld.base[-5:10] returns blank.

        In this case, we test whether the value left to ':' is negative or not.
        If it is < 0, we simply take the seq from 0 as : self.bases[0:self.cDNAStart) + 2) + s]"""
        feat = rec.features
        seq = rec.seq
        # Kozak sequence
        if int(feat["cDNAstart"]) - 1 - s < 0:
            kozakSeq = seq[0:(int(feat["cDNAstart"]) + 2) + s]
        else:
            kozakSeq = seq[(int(feat["cDNAstart"]) - 1 - s):(int(feat["cDNAstart"]) + 2)+ s]
        # Kozak context
        if int(feat["cDNAstart"]) -1 - s - c < 0:
            kozakContext = seq[0:(int(feat["cDNAstart"]) + 2)+ s + c]
        else:
            kozakContext = seq[((int(feat["cDNAstart"]) - 1 - s) - c):(int(feat["cDNAstart"]) + 2) + s + c]
        return(str(kozakSeq), str(kozakContext))

    def get3PUTR(self, rec, file):
        """Fetch 3PUTR sequence, append it to fasta file and return its length and GC content."""
        utr3p = rec.seq[int(rec.features["cDNAend"]):]
        with open(file, "a") as UTR3P:
            UTR3P.write(">{}_3PUTR\n{}\n".format(rec.id, utr3p))
        return(len(utr3p), GC(utr3p))

    def get5PUTR(self, rec, file):
        """Fetch 5PUTR sequence, append it to fasta file and return its length and GC content."""
        utr5p = rec.seq[0:int(rec.features["cDNAstart"])-1]
        with open(file, "a") as UTR5P:
            UTR5P.write(">{}_5PUTR\n{}\n".format(rec.id, utr5p))
        return(len(utr5p), GC(utr5p))

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


def get_utr5MAX(cdna_feat_row):
    #take the cdna_feat-row with multiples utrs
    #print("getting UtrMax...")
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
    #take the cdna_feat-row with multiples utrs
    #print("getting UtrMax...")
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
    #take the cdna_feat-row with multiples utrs
    #print("getting cDNA_Start_Min...")
    cDNA_start = cdna_feat_row["cDNA coding start"].values[0].split(";")
    min_cDNA_start = min(map(int, cDNA_start))
    return(min_cDNA_start)


def get_cDNAendMAX(cdna_feat_row):
    #take the cdna_feat-row with multiples utrs
    #print("getting cDNA_End_Max...")
    cDNA_end = cdna_feat_row["cDNA coding end"].values[0].split(";")
    max_cDNA_end = max(map(int, cDNA_end))
    return(max_cDNA_end)


def txt2fasta(cdna_feat_table, fastaOut):
    with open(fastaOut + ".fasta", "w+") as ff:
        for i in range(cdna_feat_table.shape[0]):
            ligne = pd.DataFrame(cdna_feat_table.loc[i,:]).transpose()
            ff.write(">{} |GeneID:{}|GeneName:{}|5P_UTR_end:{}|5P_UTR_start:{}|3P_UTR_end:{}|3P_UTR_end:{}|cDNAstart:{}|cDNAend:{}|{}\n".format(ligne["Transcript stable ID"].values[0], ligne["Gene stable ID"].values[0], ligne["Gene name"].values[0], ligne["5' UTR end"].values[0], ligne["5' UTR start"].values[0], ligne["3' UTR end"].values[0], ligne["3' UTR start"].values[0], ligne["cDNA coding start"].values[0], ligne["cDNA coding end"].values[0], ligne["Gene description"].values[0]))
            ff.write("{}\n".format(ligne["cDNA sequences"].values[0]))
    print("txt to Fasta conversion done!")


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

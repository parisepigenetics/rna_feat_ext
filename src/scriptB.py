#!/usr/bin/env python
# coding: utf-8

"""Calculate feature's Table from a fasta file with the appropriate header.

Authors: AKE Franz-Arnold - LU Antoine
mail: aerod7710@gmail.com lu.zhao.antoine@gmail.com
Juillet 2018
UMR7216 Paris Diderot"""

__version__ = "0.1a01"

# Loading Packages
import pandas as pd
import subprocess
import argparse
import sys
from Bio import SeqIO

import library
import rnaFeaturesLib


# Options/Arguments parser
parser = argparse.ArgumentParser(prog='fasta2table', description='Calculate features from a fasta file (ENSEML header) and output a featurs table for clustering', epilog="Authors: Arnold-Franz AKE, Antoine LU, 2018")
parser.add_argument("infile", help = "input FASTA file", type=str)
parser.add_argument("outfile", help = "output CSV filename; add .csv", type = str)
parser.add_argument("motifs_file", help = "MEME motifs file", default = "", type = str)
parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__))


args = parser.parse_args()

Gene_list = []
print("Collecting data...")
for record in SeqIO.parse(args.infile, "fasta"):
    new_seq = library.Seq(record)
    Gene_list.append(new_seq)
print("done")

tab = pd.DataFrame()
for feat in Gene_list:
    res = tab.append(feat.getFeatures())
    tab = res
tab.index = tab["ensembl_transcript_id"]

# Minimum Folding Energy and Features)
print("Calculate features...")

# Calcul RNAFold and MFE_per Base
print("RNAfolding calcul... P3UTR")
tab["3PMFE"] = library.RNAfold_calcul("P3UTR").astype(float)
print("RNAfolding calcul... P5UTR")
tab["5PMFE"] = library.RNAfold_calcul("P5UTR").astype(float)
tab["3PMFE_per_Base"] = tab["3PMFE"].astype(float).divide(tab["3PLen"])
tab["5PMFE_per_Base"] = tab["5PMFE"].astype(float).divide(tab["5PLen"])
print("Done")

# RBPs
subprocess.call('fimo --verbosity 1 ' + args.mtfile.name + ' ' + args.infile.name, shell=True)
fimo_tab = pd.read_csv("fimo_out/fimo.tsv", sep="\t")
fimo_tab = fimo_tab.reset_index(drop = True)

tab2 = library.RBPs_Motif(fimo_tab)
tab = tab.join(tab2, lsuffix='_tab', rsuffix='_tab2')

# Delete these files
# subprocess.call("rm -f P3UTR", shell=True)
# subprocess.call("rm -f P5UTR", shell=True)

# Sauvegarde
tab = tab.reset_index(drop=True)
print("\nSauvegarde Feature table...\n")
tab.to_csv(args.outfile, sep="\t", index=False)
print("\nDone")

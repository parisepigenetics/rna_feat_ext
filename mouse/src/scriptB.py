#!/usr/bin/env python
# coding: utf-8

# From a fasta file ..getting Feature's Table

# Author: AKE Franz-Arnold - Antoine LU
# mail: aerod7710@gmail.com
# 17 Juillet 2018
# UMR7216 Paris Diderot

# Loading Packages
import library
import pandas as pd
import numpy as np
import subprocess
import argparse
import sys
import re
from Bio import SeqIO

# Options/Arguments parser
parser = argparse.ArgumentParser(description='Fetch features from an Fasta file and output a Feature Tab for Clustering', epilog="Author: Franz-Arnold Ake and Antoine Lu, 2018")

# Positional Arguments
parser.add_argument('infile', metavar="input_file", type=argparse.FileType('r'), help='Path for Fasta file')
parser.add_argument('mtfile', nargs="?", metavar="motif_file", type=argparse.FileType('r'), help='Motif file for determining RBPs')
parser.add_argument('outfile', nargs="?", default=sys.stdout, metavar="output_file", type=argparse.FileType('w'), help='Path for saving Output Feature table file (Default: output in Stdout)')

args = parser.parse_args()

# Delete these files
subprocess.call("rm -f P3UTR", shell=True)
subprocess.call("rm -f P5UTR", shell=True)
#FIXME WHY we have these files in the first place I think you delete them later (code duplication.)

# First Loop
Gene_list = []
print("collecting data...")
for record in SeqIO.parse(args.infile, "fasta"):
    new_seq = library.Seq(record)
    Gene_list.append(new_seq)
print("done")


# 2nd Loop (Minimum Folding Energy and Features)
print("Calcul features...")
tab = pd.DataFrame()
for feat in Gene_list:
    res = tab.append(feat.getFeatures())
    tab = res
tab.index = tab["ensembl_transcript_id"]

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
subprocess.call("rm -f P3UTR", shell=True)
subprocess.call("rm -f P5UTR", shell=True)


# Sauvegarde
tab = tab.reset_index(drop=True)
print("\nSauvegarde Feature table...\n")
tab.to_csv(args.outfile, sep="\t", index=False)
print("\nDone")

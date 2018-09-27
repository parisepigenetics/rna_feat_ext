#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""Calculate feature's Table from a fasta file with the appropriate header.

Authors: AKE Franz-Arnold - LU Antoine
mail: aerod7710@gmail.com lu.zhao.antoine@gmail.com
Juillet 2018
UMR7216 Paris Diderot"""

__version__ = "0.1a01"

#Loading packages
import pandas as pd
import argparse
import subprocess
import sys

import rnaFeaturesLib


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='fasta2table', description='Calculate features from a fasta file (ENSEML header) and output a featurs table for clustering', epilog="Authors: Arnold-Franz AKE, Antoine LU, 2018")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument("infile", help = "input FASTA file", type=str)
    parser.add_argument("outfile", help = "output CSV filename; add .csv", type = str)
    parser.add_argument("motifs_file", help = "MEME motifs file", default = "", type = str)

    args = parser.parse_args()

    featExt = rnaFeaturesLib.FeatExtract(args.fasta_file)
    feaTable = featExt.dicos2table()
    columns = ['ensembl_gene_id', 'ensembl_transcript_id', '3PLen', '3PMfe', '5PLen', '5PMfe', '5UTR_mfe_Base', '3UTR_mfe_Base', 'Kozak_Context', 'Kozak_sequence']
    feaTable.reindex_axis(columns, axis=1)
    feaTable.to_csv(args.out_csv + ".csv", sep="\t", index = False)

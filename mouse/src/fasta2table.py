#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#Loading packages
print("import packages")
import pandas as pd
import rnaFeaturesLib
import argparse
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='fasta2table')
    parser = argparse.ArgumentParser(description='Fasta to table conversion', epilog="Author: Antoine Lu, 2018")
    parser.add_argument('--version', action='version', version='{} {}'.format('fasta2table', 0.1))
    parser.add_argument("fasta_file", help="input FASTA file", type=str)
    parser.add_argument("out_csv", help="output CSV filename; add .csv", type=str)
    
    args = parser.parse_args()

    start = time.time()

    featExt = rnaFeaturesLib.FeatExtract(args.fasta_file)
    feaTable = featExt.dicos2table()
    columns = ['ensembl_gene_id', 'ensembl_transcript_id', '3PLen', '3PMfe', '5PLen', '5PMfe', '5UTR_mfe_Base', '3UTR_mfe_Base', 'Kozak_Context', 'Kozak_sequence']
    feaTable.reindex_axis(columns, axis=1)
    feaTable.to_csv(args.out_csv + ".csv", sep="\t", index = False)

    end = time.time()
    print("Run time: {}".format(end - start))
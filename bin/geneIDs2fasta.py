#!/usr/bin/python3
# -*- coding: UTF-8 -*-

"""From a list of ENSEMBL Gene_IDs get a fasta file with metadata as header.

Author: AKE Franz-Arnold - Antoine LU
mail: aerod7710@gmail.com lu.zhao.antoine@gmail.com
July 2018
UMR7216 Paris Diderot"""

__version__ = "0.3a"

import argparse
import rnaFeaturesLib

parser = argparse.ArgumentParser(prog='geneIDs2Fasta', description='Fetch features for an ENSEMBL gene ID list and output a multi FASTA file.', epilog="Authors: Franz-Arnold Ake and Antoine Lu, 2018")
parser.add_argument('infile', metavar="input_file", type=argparse.FileType('r'), help='Path for Ensembl Gene_ID file')
parser.add_argument("fasta_out", help="output FASTA file name, the suffix '.fasta' is automatically added", type=str)
parser.add_argument('-d', '--dataset', nargs="?", default='hsapiens_gene_ensembl', metavar="ENSEMBL Dataset for collecting information", type=str, help="Choice from Ensembl Datasets (taken from the web API of ENSEMBL) -- Default : hsapiens_gene_ensembl")
parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__))

args = parser.parse_args()

# This is the function wich does all the job.
listID = args.infile.read().splitlines()
rnaFeaturesLib.geneIDs2Fasta(listID, args.fasta_out, args.dataset)

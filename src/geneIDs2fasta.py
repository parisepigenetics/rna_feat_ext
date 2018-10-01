#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""From a list of ENSEMBL Gene_IDs get a fasta file with metadata as header.

Author: AKE Franz-Arnold - Antoine LU
mail: aerod7710@gmail.com lu.zhao.antoine@gmail.com
July 2018
UMR7216 Paris Diderot"""

__version__ = "0.1a01"

#Loading packages
import argparse
from biomart import BiomartServer
import pandas as pd
from pandas.compat import StringIO

import rnaFeaturesLib

# Functions
def geneIDs2Fasta(listID, out_fasta, dataset):
    print("Connection to server....")
    server = BiomartServer( "http://www.ensembl.org/biomart/")
    print("....done!")
    dt = server.datasets[dataset]
    print("Connexion to ENSEMBL dataset.")

    #TODO provide the required attributes in an external parameters file.
    listAttrib = ['ensembl_gene_id', 'ensembl_transcript_id',
                  'external_gene_name', 'transcript_start', 'transcript_end',
                  '5_utr_end', '5_utr_start', '3_utr_end', '3_utr_start',
                  'transcription_start_site', 'transcript_biotype',
                  'cdna_coding_start', 'cdna_coding_end', 'cdna', 'description']
    print("Fetching data.....")
    response = dt.search({'filters': {'ensembl_gene_id': listID}, 'attributes' : listAttrib}, header = 1)
    print("...fetching done!")

    # Convert stringIO to pandas data frame.
    print("...Convert string to PD.data.frame.")
    stream = StringIO(response.text)
    cdna = pd.read_table(stream, sep="\t")

    # remove NA and reset index
    cdna = cdna.dropna()
    cdna = cdna.reset_index()

    #TODO TO TEST all these before publishing anything!
    # select only longest 5'UTR and 3'UTR
    print("...Select the longest UTRs!")
    for i in range(cdna.shape[0]):
        ligne = pd.DataFrame(cdna.loc[i,:]).transpose()
        #For 5_UTR_start
        indices = rnaFeaturesLib.get_utr5MAX(ligne)
        if(len(cdna.iloc[i:i+1,:]["5' UTR start"].values[0].split(";")) > 1):
            cdna.iloc[i:i+1,7:8] = indices[1]
            cdna.iloc[i:i+1,8:9] = indices[0]
        if(len(cdna.iloc[i:i+1,:]["3' UTR start"].values[0].split(";")) > 1):
            cdna.iloc[i:i+1,9:10] = indices[1]
            cdna.iloc[i:i+1,10:11] = indices[0]

    # select longest cDNA
    print("...Select the longesst cDNA")
    for i in range(cdna.shape[0]):
        ligne = pd.DataFrame(cdna.loc[i,:]).transpose()
        # smallest cDNA start value
        min_cdna_start = rnaFeaturesLib.get_cDNAstartMIN(ligne)
        # therefore: biggest cDNA end value
        max_cdna_end = rnaFeaturesLib.get_cDNAendMAX(ligne)
        cdna.iloc[i, cdna.columns.get_loc('cDNA coding start')] = min_cdna_start
        cdna.iloc[i, cdna.columns.get_loc('cDNA coding end')] = max_cdna_end

    # Exportat to FASTA format
    print("Return FASTA file!")
    rnaFeaturesLib.txt2fasta(cdna, out_fasta)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='geneIDs2Fasta', description='Fetch features for an ENSEMBL gene ID list and output a multi FASTA file.', epilog="Authors: Franz-Arnold Ake and Antoine Lu, 2018")
    parser.add_argument('infile', metavar="input_file", type=argparse.FileType('r'), help='Path for Ensembl Gene_ID file')
    parser.add_argument("fasta_out", help="output FASTA file name, the suffix '.fasta' is automatically added", type=str)
    parser.add_argument('-d', '--dataset', nargs="?", default='hsapiens_gene_ensembl', metavar="ENSEMBL Dataset for collecting information", type=str, help="Choice from Ensembl Datasets (taken from the web API of ENSEMBL) -- Default : hsapiens_gene_ensembl")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=__version__))


    args = parser.parse_args()
    #This is the function wich does all the job.
    listID = args.infile.read().splitlines()
    geneIDs2Fasta(listID, args.fasta_out, args.dataset)

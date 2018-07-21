# -*- coding: utf-8 *-*

def geneIDs2Fasta(id_list, out_fasta):
    #Loading packages
    import numpy as np
    import pandas as pd
    import math
    import re
    import sys
    import subprocess
    import rnaFeaturesLib
    import argparse
    from biomart import BiomartServer

    server = BiomartServer( "http://www.ensembl.org/biomart/")
    server.show_datasets()

    dt = server.datasets['mmusculus_gene_ensembl']
    # remplacer par 'hsapiens_gene_ensembl'

    # run a search with custom filters and attributes (no header)
    # 3Plen, 3PEng, 5Plen, 5PEng, 3PNorm, 5PNorm, dORF, uORF, KozakSeq, KozakCont

    with open(id_list, "r") as filin:
        listID = filin.read().splitlines()

    filout = open("tmp.tsv", "w")

    response = dt.search({
    'filters': {
        'ensembl_gene_id': listID
    },
    'attributes': [
        'ensembl_gene_id', 'ensembl_transcript_id','external_gene_name','transcript_start', 
        'transcript_end', '5_utr_end', '5_utr_start', '3_utr_end', '3_utr_start',
        'transcription_start_site', 'transcript_biotype', 'cdna_coding_start',
        'cdna_coding_end','cdna'
    ]
    }, header = 1 )

    filout.write(response.text)
    filout.close()

    # reading data
    cdna = pd.read_csv("tmp.tsv", sep="\t")

    # remove NA and reset index
    cdna = cdna.dropna()
    cdna = cdna.reset_index()

    # select only longest 5'UTR and 3'UTR
    for i in range(cdna.shape[0]):
        ligne = pd.DataFrame(cdna.loc[i,:]).transpose()
        #For 5_UTR_start
        indices = get_utr5MAX(ligne)
        if(len(cdna.iloc[i:i+1,:]["5' UTR start"].values[0].split(";")) > 1):
            cdna.iloc[i:i+1,7:8] = indices[1]
            cdna.iloc[i:i+1,8:9] = indices[0]
        if(len(cdna.iloc[i:i+1,:]["3' UTR start"].values[0].split(";")) > 1):
            cdna.iloc[i:i+1,9:10] = indices[1]
            cdna.iloc[i:i+1,10:11] = indices[0]

    # select longest cDNA
    for i in range(cdna.shape[0]):
        ligne = pd.DataFrame(cdna.loc[i,:]).transpose()
        # smallest cDNA start value
        min_cdna_start = rnaFeaturesLib.get_cDNAstartMIN(ligne)
        # therefore: biggest cDNA end value
        max_cdna_end = rnaFeaturesLib.get_cDNAendMAX(ligne)
        cdna.iloc[i, cdna.columns.get_loc('cDNA coding start')] = min_cdna_start
        cdna.iloc[i, cdna.columns.get_loc('cDNA coding end')] = max_cdna_end

    # Exportation to FASTA format
    rnaFeaturesLib.txt2fasta(cdna, out_fasta)

    return(0)

    if __name__ == '__main__':

        program = 'geneIDs2Fasta'
        version = 0.1
        description = 'Fetch features for an ensembl gene ID list and output a multi FASTA file'
        parser = argparse.ArgumentParser(prog=program)
        parser = argparse.ArgumentParser(description=description, epilog="Author: Franz-Arnold Ake and Antoine Lu, 2018")
        parser.add_argument('--version', action='version', version='{} {}'.format(program, version))
        parser.add_argument("i", help="input id file", type=str)
        parser.add_argument("o", help="output FASTA file name, the suffix '.fasta' is automatically added", type=int)

        args = parser.parse_args()
        geneIDs2Fasta(args.i, args.o)

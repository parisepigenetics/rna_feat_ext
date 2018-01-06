#! /usr/bin/env python
# -*- coding: utf-8 -*-

# MODULES
import sys
import time
import argparse
from BCBio import GFF

# PARSER

parser = argparse.ArgumentParser()
parser.add_argument("in_file", help="INPUT filename (with or without path), gff3 format.", type=str)
parser.add_argument("out_file", help="OUTPUT filename (with or without path).", type=str)
parser.add_argument("-t", "--id_type", help="define what type of sequence to retrive id from.\
    \nExpected value are either 'mRNA' (default) or 'CDS'.", type=str, default="mRNA", choices=["mRNA", "CDS"])
args = parser.parse_args()
if args.id_type:
    print "id_type value is {0}\n".format(args.id_type)


# FUNCTION


def extract_ids(in_file, out_file, id_type):
    '''
    This function aims to extract coding sequences transcript ids using
    biopython's GFF parser.
    Both input file and output file have to be given by the user.
    Takes a GFF file as input, return a txt file containing one id per line.
    ---
    Default value of the id_type parameter is "mRNA" : if we don't exclude
    introns and UTRs the third agument, id_type can take these values : mRNA,
    CDS
    '''
    limit_info = dict(gff_type=[str(id_type)])
    with open(in_file, "r") as in_handle, open(out_file, "w+") as filout:
        tmp_rec = "none"
        for rec in GFF.parse(in_handle, limit_info=limit_info, target_lines=1):
            if tmp_rec != rec.features[0].id:
                filout.write(rec.features[0].id.split(":")[1].strip('"')+"\n")
                tmp_rec = rec.features[0].id
    print "Extraction done.\n"
    return


if __name__ == '__main__':
    start = time.time()
    extract_ids(sys.argv[1], sys.argv[2], args.id_type)

    end = time.time()
    print(end - start)

# USAGE
'''
If ids to extract are from CDS sequences :
python extract_ids.py ../data/Homo_sapiens.GRCh38.90.gff3
                        ../results/Homo_sapiens.GRCh38.90_CDS_id.txt -t CDS

Default : extraction of mRNA sequences :
python extract_ids.py ../data/Homo_sapiens.GRCh38.90.gff3
                        ../results/Homo_sapiens.GRCh38.90_mRNA_id.txt
'''

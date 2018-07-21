# README

authors: Arnold AKE and Antoine LU

# geneIDs2fasta.py

This script takes a text list of ENSEMBL gene IDs and returns a FASTA formated file of the corresponding cDNA sequences. The header is formatted and contains various metadata ordered as:

>Gene stable ID, Transcript stable ID, Gene name, 5' UTR end, 5' UTR start, "3' UTR end, "3' UTR start, cDNA coding start, cDNA coding end

INSIDE A PYTHON TERMINAL
```python
geneIDs2fasta.genesIDs2Fasta(input, output)
```
Works fine, but doesn't as bellow... WHY??
#TODO
```python
python geneIDs2fasta.py input_IDs_file output_FASTA_filename
```   

# rnaFeaturesLib.py

Python file containing all the called functions


## ISSUES:
When lauching the python command line (see previous) nothing happens: it seems like no connection is estabished to biomart ensembl server? 


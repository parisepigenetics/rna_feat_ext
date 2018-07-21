# README

authors: Arnold AKE and Antoine LU

# geneIDs2fasta.py

This script takes a text list of ENSEMBL gene IDs and returns a FASTA formated file of the corresponding cDNA sequences. The header is formatted and contains various metadata ordered as:

>Gene stable ID, Transcript stable ID, Gene name, 5' UTR end, 5' UTR start, "3' UTR end, "3' UTR start, cDNA coding start, cDNA coding end

```python
python geneIDs2Fasta input_IDs_file output_FASTA_filename
```   

# rnaFeaturesLib.py

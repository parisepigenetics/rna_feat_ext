# README

authors: Antoine LU

# geneIDs2fasta.py

This script takes a text list of ENSEMBL gene IDs and returns a FASTA formated file of the corresponding cDNA sequences. The header is formatted and contains various metadata ordered as:

>Gene stable ID, Transcript stable ID, Gene name, 5' UTR end, 5' UTR start, "3' UTR end, "3' UTR start, cDNA coding start, cDNA coding end

```python
python geneIDs2fasta.py input_IDs_file output_FASTA_filename
```

For example:

```python
python geneIDs2fasta.py bics10Ids.txt bics10Ids
```
> returns: bics10Ids.fasta

# fasta2table.py

This script takes the previously returned FASTA file and transforms it into a tabulation separated dataframe (CSV format).

```python
pyhton fasta2table.py input_FASTA_file output_CSV_filename
```

For example:

```python
python fasta2table.py bics10Ids.fasta bics10Ids
```
> returns: bics10Ids.csv

# rnaFeaturesLib.py

Python file containing all the called functions for the two main scripts


## ISSUES:
NONE for now..

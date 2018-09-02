# README

These scripts are part of a project supervised by the PTER Team (@parisepigenetics, @cbouyio) at Paris Diderot. This collaboration now involves two students: **Antoine LU @antoinezl** who started the project as part of a project during his second year of Master's degree and **Franz-Arnold AKE @franzx5** a first year Master's degree student who mainly worked on the clustering aspect of the project (see clustering branch).

## Development of a computational tool for RNA feature extraction

The aim of the project is to develop a tool for RNA feature extraction. These features might have a role in identifying commonly regulated mRNA.

## Data
+ /bics/best_bics_id.txt: Best BICS IDs determined by a post-doc student (Mohamed MACHAT).

## Modules / Tools
(Python 2.7)
+ [`Biopython`](http://biopython.org)
+ [`Pandas`](http://pandas.pydata.org)
+ [`biomart`API](https://pypi.org/project/biomart/)

(Tools and script)
+ [`RNAfold, Vienna RNA Package`](https://www.tbi.univie.ac.at/RNA/index.html#download)(version 2.4.3)
+ [`motif_scan.py`](https://github.com/miha-skalic/motif_scan)

## Installation
+ [`motif_scan.py`](https://github.com/miha-skalic/motif_scan): Please follow the link, the procedure is well detailed.
+ [`RNAfold, Vienna RNA Package`](https://www.tbi.univie.ac.at/RNA/index.html#download): Please follow the link, the procedure is well detailed.

## Usage and descriptions
# geneIDs2fasta.py or scriptA.py

	authors: Antoine LU and Arnold AKE

This script takes a text list of ENSEMBL gene IDs and returns a FASTA formated file of the corresponding cDNA sequences. The header is formatted and contains various metadata ordered as:

>Gene stable ID, Transcript stable ID, Gene name, 5' UTR end, 5' UTR start, "3' UTR end, "3' UTR start, cDNA coding start, cDNA coding end

```python
python geneIDs2fasta.py input_IDs_file output_FASTA_filename
```
or:
```Python
python scriptA.py input_IDs_file output_FASTA_filename
```

For example:

```python
python geneIDs2fasta.py bics10Ids.txt bics10Ids
```
> returns: bics10Ids.fasta

# fasta2table.py or scriptB.py

	author: Antoine LU and Franz-Arnold AKE

This script takes the previously returned FASTA file, MOTIF file, and transforms it into a tabulation separated dataframe (CSV format).

```python
pyhton fasta2table.py input_FASTA_file output_CSV_filename
```

For example:

```python
python fasta2table.py bics10Ids.fasta bics10Ids
```
or:
```python
python scriptB.py bics10Ids.txt  Motif_TXT_file bics10Ids
```
> returns: bics10Ids.tab

# rnaFeaturesLib.py or library.py

Python file containing all the called functions for the two main scripts

# Searching motif

we use for searching motif_id per for each transcript the tool FIMO(find individual motif occurences) of the MEME Suite.
we have to install this tool before to execute the scriptB.py which use it.

Check this link -> + ['fimo'](http://meme-suite.org/doc/fimo.html)


## ISSUES:
NONE for now..

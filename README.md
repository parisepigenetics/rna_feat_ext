## Development of a computational tool for RNA feature extraction
# Projet Long M2BI - 2017-2018
The aim of the project is to develop a tool for RNA feature extraction. These features might have a role in identifying commonly regulated mRNA.

## Modules / Tools
(Python 2.7)
+ [`Biopython`](http://biopython.org)
+ [`Pandas`](http://pandas.pydata.org)

(R 3.3.1)
+ [`biomaRt`](http://bioconductor.org/packages/release/bioc/html/biomaRt.html)

(Tools and script)
+ [`RNAfold, Vienna RNA Package`](https://www.tbi.univie.ac.at/RNA/index.html#download)(version 2.4.3)
+ [`motif_scan.py`](https://github.com/miha-skalic/motif_scan)

## Installation
+ [`motif_scan.py`](https://github.com/miha-skalic/motif_scan): Please follow the link, the procedure is well detailed.
+ [`RNAfold, Vienna RNA Package`](https://www.tbi.univie.ac.at/RNA/index.html#download): Please follow the link, the procedure is well detailed.

## Usage
Because the integration of the scripts couldn't be done in time the scripts have to be executed separately:
### Transcript ID extraction:

Default: the `extract_ids.py` script returns the IDs of mRNA
```
python extract_ids.py ../data/Homo_sapiens.GRCh38.90.gff3 ../results/Homo_sapiens.GRCh38.90_mRNA_id.txt
```
If the CDS' IDs is required, please mention `-t CDS`.
```
python extract_ids.py ../data/Homo_sapiens.GRCh38.90.gff3 ../results/Homo_sapiens.GRCh38.90_CDS_id.txt -t CDS
```

### Dataset Preparation:
Two executable Rscripts have been written for two types of datasets: 
+ [cDNA_write.R] (rna_feat_ext/src/cDNA_write.R): for cDNA sequences (many transcripts sequences for each gene ID)
+ [metagene_script.R] (rna_feat_ext/src/metagene_script.R): for 'metagene' sequences (one transcript sequence for each gene ID) We defined a metagene to be the union of every possible exons of a gene.

In a terminal, to create both a TEXT file containing the fetched sequences and a CSV table (that will be used next) do:
For cDNA sequences: 
Please make sure not to put the extensions to your input and output files
```
Rscript cDNA_write.R <input finename> <output name>
Rscript cDNA_write.R ../results/10112017/Homo_sapiens.GRCh38.90_mARN_id SeqCDNAOut
```
The same syntax has to be applied for the following:
```
Rscript metagene_script.R <input finename> <output name>
Rscript metagene_script.R ../results/10112017/Homo_sapiens.GRCh38.90_mARN_id SeqMetaOut
```
### RNA feature extraction:
Extraction of the features can be done using:
Jupyter notebook: (tested)
```
jupyter notebook rna_feature_ext.ipynb
```
In a terminal (not tested yet, but the script is similar, might have to change the paths)
Optimization to be brought: management of the arguments using Python's argparse module
```
rna_feature_ext.py
```






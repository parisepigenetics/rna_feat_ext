# mRNA feature extraction tools

    authors: Costas BOUYIOUKOS, Antoine LU and Arnold Franz AKE

A set of computational tools to extract user defined mRNA features from a list of ENSEMBL gene IDs by using the web API of ENSEMBL BioMart and custom computations.
Conceived and developed by [**Costas Bouyioukos**](https://github.com/cbouyio) @cbouyio at [Paris Epigenetics](https://github.com/parisepigenetics) @parisepigenetics and Universite Paris Diderot.
Development involved two bioinformatics master students: [**Antoine LU**](https://github.com/antoinezl) @antoinezl who started as part of a coding project during his second year in the degree and [**Franz-Arnold AKE**](https://github.com/franzx5) @franzx5 a second year Master's degree student who mainly worked on the clustering part of the project.


## Installation.
To install the tools in your local python environment (user $HOME directory) type:

```shell
./setup.py install --user
```

(the --user flag installs the software on your personal account (no root privileges required).

### Requirements.
#### Python.
  - [Biopython](http://biopython.org/)
  - [Pandas](http://pandas.pydata.org/)
  - [Biomart](https://pypi.org/project/biomart/)

All are available for installation via `pip install <package_name>`

#### External.
  - [RNA Vienna package](https://www.tbi.univie.ac.at/RNA/)
  - [MEME Suite](http://meme-suite.org/)

For external tools please follow the installation guidelines in the provided links.


## Main Usage.
    geneIDs2fasta.py ENSEMBL_geneIDs_file fasta_output_file

and

    fasta2table.py ENSEMBL_fasta_output_file features_table_file

### geneIDs2fasta.py
This program takes a text file with a list of ENSEMBL gene IDs and returns a FASTA formatted file of the corresponding cDNA sequences. The header is formatted and contains various metadata ordered as:

`>ENSEMBL_transcript_ID |Gene stable ID | Gene name | cDNA start | cDNA end | TSL | APRIS | HAVANA_ENSEMBL | gene description | Source:|`

### fasta2table.py
This program takes the fasta formatted file returned by the previous script geneIDs2fasta in input, and return a semicolon separated table with the following header:

`ensembl_gene_id;gene_name;coding_len;5pUTR_len;5pUTR_GC;5pUTR_MFE;5pUTR_MfeBP;3pUTR_len;3pUTR_GC;3pUTR_MFE;3pUTR_MfeBP;TOP_localScore;CAI;Kozak_Sequence;Kozak_Context`


## Testing
Test directory contains two test files to test and demonstrate the functionality of the tools.

+ test/testENSEMBLids.txt Contains 6 genes with their ENSEMBL IDs.

+ test/testTransExpr.csv Contains the expression levels of each individual transcript of the above genes from a case study.

##### TODO add section for MEME suite integration.

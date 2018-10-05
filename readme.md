# mRNA feature extraction tools

    authors: Costas BOUYIOUKOS, Antoine LU and Arnold Franz AKE

A set of computational tools to extract uder defined mRNA features from a list of ENSEMBL gene IDs by using the web API of ENSEMBL BioMart and custom computations. Conceived and developed by [**Costas Bouyioukos**](https://github.com/cbouyio) @cbouyio at [Paris Epigenetics](https://github.com/parisepigenetics) @parisepigenetics and Universite Paris Diderot. Development involved two bioinformatics master students: [**Antoine LU**](https://github.com/antoinezl) @antoinezl who started as part of a coding project during his second year in the degree and [**Franz-Arnold AKE**](https://github.com/franzx5) @franzx5 a second year Master's degree student who mainly worked on the clustering part of the project.


## Installation.

    setup.py install --user=$HOME
Add the --user flag to install it on your personal account (no root priviledges required).

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
    geneIDs2fasta.py ENSEMBL_geneIDs.txt fasta_file_prefix

and

    fasta2table.py fasta.fasta features_table.csv RNA_MEME_motifs

### geneIDs2fasta.py
This script takes a file text of ENSEMBL gene IDs and returns a FASTA formated file of the corresponding cDNA sequences. The header is formatted and contains various metadata ordered as:

`> Gene stable ID, Transcript stable ID, Gene name, 5' UTR end, 5' UTR start, "3' UTR end, "3' UTR start, cDNA coding start, cDNA coding end, gene description`

### fasta2table.py
This script takes the fasta formated file returned by the previous script geneIDs2fasta in input, and return a semicolon separated table with the following header:

`ensembl_gene_id;gene_name;coding_len;5pUTR_len;5pUTR_GC;5pUTR_MFE;5pUTR_MfeBP;3pUTR_len;3pUTR_GC;3pUTR_MFE;3pUTR_MfeBP;Kozak_Sequence;Kozak_Context`


#TODO Clean up data directory and rewrite this section.
## Data
+ data/bics/best_bics_id.txt: Best BICS IDs determined by a post-doc student (Mohamed MACHAT).

+ data/bics/bicsIds.txt:   some bics ID with relevant information on Ensembl DATABASE.

+ data/motif_databases

  (**Motif_file Database for MEME Suite**) can be also downloaded [`here`](http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.18.tgz)

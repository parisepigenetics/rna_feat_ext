# Development ~ Computational tool for RNA feature extraction

    authors: Antoine LU and Arnold Franz AKE

These scripts are part of a project supervised by the PTER Team (@parisepigenetics, @cbouyio) at Paris Diderot. This collaboration now involves two students: **Antoine LU @antoinezl** who started the project as part of a project during his second year of Master's degree and **Franz-Arnold AKE @franzx5** a second year Master's degree student who mainly worked on the clustering aspect of the project (see clustering branch).


We develop two scripts for achieve this purpose :

**geneIDs2fasta.py** and **fasta2table.py**

### geneIDs2fasta.py
This script takes a file text of ENSEMBL gene IDs and returns a FASTA formated file of the corresponding cDNA sequences. The header is formatted and contains various metadata ordered as:

> Gene stable ID, Transcript stable ID, Gene name, 5' UTR end, 5' UTR start, "3' UTR end, "3' UTR start, cDNA coding start, cDNA coding end

### fasta2table.py
This script takes the fasta formated file returned by the previous script geneIDs2fasta in input, and return an output_table with following variables:

> ensembl\_gene\_id, ensembl\_transcript\_id, uORF, 5PLen, dORF, 3PLen, Kozak\_Context, Kozak\_Sequence, 3PMFE, 5PMFE, 3PMFE\_per\_base, 5\_PMFE\_per\_base, motif_ID

## Data
+ data/bics/best_bics_id.txt: Best BICS IDs determined by a post-doc student (Mohamed MACHAT).

+ data/bics/bicsIds.txt:   some bics ID with relevant information on Ensembl DATABASE.

+ data/motif_databases

  (**Motif_file Database for MEME Suite**) can be also downloaded [`here`](http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.18.tgz)


## Modules / Tools
(Python 2.7)

+ [`Biopython`](http://biopython.org)
+ [`Pandas`](http://pandas.pydata.org)
+ [`biomart`API](https://pypi.org/project/biomart/)
+ [`MEME SUITE`](http://meme-suite.org/meme-software/5.0.2/meme-5.0.2.tar.gz)

(Tools and script)

+ [`RNAfold, Vienna RNA Package`](https://www.tbi.univie.ac.at/RNA/)(version 2.4.9)


## Installation

Theses package are required for the correct execution of the scripts:

+ [`MEME SUITE`](http://meme-suite.org/doc/install.html?man_type=web): Please follow the link, the procedure is well detailed.

+ [`RNAfold, Vienna RNA Package`](https://www.tbi.univie.ac.at/RNA/index.html#download): Please follow the link, the procedure is well detailed.



## Usage and descriptions (How to ...)

### geneIDs2fasta.py or geneIDsfasta_2.py (Franz)

`Unix Environment`

    ./geneIDs2fasta.py input_geneIDs_text output_fasta_filename


For example (*in Unix Terminal*):

    ./geneIDs2fasta.py ../data/example/bicsIds.txt bicsIds.fasta
> returns: bicsIds.fasta


### fasta2table.py or fasta2table_2.py (Franz)

`Unix Environment`

    ./fasta2table.py input_fasta_file motif_file output_tab_file


For example (*in Unix Terminal*):

    ./geneIDs2fasta.py bicsIds.fasta ../data/example/jolma2013.meme bicsIds.tab
> returns: bicsIds.tab

###########################
# Antoine LU              #
# M2BI - Projet Long 2017 #
###########################

# library
library(biomaRt)

# Mart (show available marts):
listMarts()

# we want homo sapiens genes from ensembl database:
ensembl = useMart("ensembl")
listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# split the big file into smaller ones (4,5k lines per file):
system("split -l 4500 Homo_sapiens.GRCh38.90_mARN_id.txt")

# Parsing the database to retrieve sequences :
split_files = system("ls ../results/x*", intern = TRUE)
tmp = NULL
final = NULL

# transcript_exon_intron : full transcript unspliced
# coding_transcript_flank : flanking region of the transcript including the UTRs

for(f in seq(length(split_files))){
  transcript_id = read.table(tmp[f])
  tmp = getSequence(id = transcript_id[,1], 
                     type="ensembl_transcript_id",
                     seqType="transcript_exon_intron", 
                     mart=ensembl)
  final = rbind(final, tmp)
  print(paste("done",f,"/",length(tmp)))
  tmp = NULL
}

# export the sequences into new fasta file
exportFASTA(final, "../results/BioMart/withR/unspliced_transcript_mRNA.txt")

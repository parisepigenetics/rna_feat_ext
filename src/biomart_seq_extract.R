###########################
# Antoine LU              #
# M2BI - Projet Long 2017 #
###########################
# References
# Folding energy calculator:
#http://rna.urmc.rochester.edu/RNAstructureWeb/Servers/AllSub/AllSub.html
#https://www.tbi.univie.ac.at/RNA/tutorial/
  
# library
library("biomaRt")

# Mart (show available marts): ## NEEDS INTERNET
listMarts()
listEnsembl(version=90) # version actuelle: 91/ 90: version du genome utilisé
# we want homo sapiens genes from ensembl database:
ensembl90 = useEnsembl(biomart="ensembl",version=90)
mart = useMart("ensembl")

ld = listDatasets(ensembl90)
View(ld)
ld[which(ld[,]=="hsapiens_gene_ensembl"),1]
mart = useDataset("hsapiens_gene_ensembl",mart=ensembl90)

# Header information: Gene stable ID; Gene Description; Gene name
# + Transcript information: CDS start (within cDNA); CDS end (within cDNA)
# + 5' UTR start/end; 3' UTR start/end; transcript stable ID;
# + Transcript start/end (bp); Transcript length (including UTRs and CDS)
# + Exon stable ID; Exon region start (bp)

filters = listFilters(mart)
attributes = listAttributes(mart)
attribs = attributes[c(1,3,8,14,15,16,17,21,200,201,202,203, 211,1694,1695,1702,1707),1]
#attrib = c("ensembl_gene_id","ensembl_transcript_id", "description", "transcript_start", "transcript_end", "transcription_start_site", "transcript_length", "external_gene_name", "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end", "cdna_coding_start", "cdna_coding_end", "ensembl_exon_id", "exon_chrom_start")
#getBM(attributes = c("refseq_dna", "interpro", "interpro_description"), filt + values = refseqids, mart = ensembl)

# split the big file into smaller ones (4,5k lines per file):
system("split -l 4500 Homo_sapiens.GRCh38.90_mARN_id.txt")

# Parsing the database to retrieve sequences :
split_files = system("ls ../results/split_mRNA_file/x*", intern = TRUE)
tmp = NULL
final = NULL

# transcript_exon_intron: full transcript unspliced
# coding_transcript_flank: flanking region of the transcript including the UTRs

for(f in seq(length(split_files))){
  transcript_id = read.table(split_files[f])
  tmp = getBM(
    attributes = attribs,
    filters = 'ensembl_transcript_id', # order by (as in sql)
    values = transcript_id[,1],
    mart = mart
    )
  final = rbind(final, tmp)
  print(paste("done",f,"/",length(split_files)))
  #print(tmp)
  tmp = NULL
}

# idée : récupérer les informations du header avec getBM()
# récupérer les séquences pour les identifiants soit de transcript_id[,1]
# soit de la bonne colonne 'ensembl_transcript_id' de la matrice de getBM()

for(f in seq(length(split_files))){
  transcript_id = read.table(split_files[f])
  tmp = getSequence(id = transcript_id[,1], 
                     type="ensembl_transcript_id",
                     seqType="transcript_exon_intron", #cdna (que exon et utrs) ou # gene exon
                     mart=mart)
  final = rbind(final, tmp)
  print(paste("done",f,"/",length(split_files)))
  tmp = NULL
}

# export the sequences into new fasta file
exportFASTA(final, "../results/BioMart/withR/unspliced_transcript_mRNA.txt")



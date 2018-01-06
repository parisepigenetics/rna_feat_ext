###########################
# Antoine LU              #
# M2BI - Projet Long 2017 #
###########################
# References
# Folding energy calculator:
#http://rna.urmc.rochester.edu/RNAstructureWeb/Servers/AllSub/AllSub.html
#https://www.tbi.univie.ac.at/RNA/tutorial/

#' @description  
#' This script takes as input the text file containing the transcript ids retrieved using python
#' ../results/10112017/Homo_sapiens.GRCh38.90_mARN_id.txt"
#' the text file is split into smaller files so they can be read and parse with R
#'
#'  The biomaRt library is used. Since the data used in the first place is the Human genome version 90,
#'  this mart will be used here.
#'  
#'  attribs: every metadata wanted in the HEADER of the FASTA sequences
#'  
#'  METADATA table (same table as in cDNA_write.R): the getBM function is used and takes the variable
#'  "attribs" as an input to create a dataframe with every metadata for each of the transcripts.
#'  Since exon features (rank, 5 prime and 3 prime UTRs) are retrieved, the table has this dimension:
#'  843,569 lines and 15 columns
#'  
#'  DEFINITON: A metagene is the union of all the exon of one gene
#'  
#'  METAGENE table: For each transcript id, its length is retrieved using the getBM() function: trlength table.
#'  Then, the listTranscriptId() function takes the previous table as an input (as well as a file name to contain
#'  the results). This function fetches the longest sequences for a given gene_id. At the end: 1 transcript for 1
#'  gene.
#'  The returned txt file, is split in smaller files (4500 lines per file). Then parsed using the getBM sequence
#'  in the attribToMerge() function: retrieves the gene_exon sequence, the gene_id, transcript_id and rank.
#'  A new table is returned, named "ToWrite".
#'      MERGED DF METAGENE: (variable named "ToWrite") : Combination of both previous data frames on the
#'      "ensembl_transcript_id","ensembl_gene_id" and "rank" columns. 
#'      Creates a new dataframe with 843,569 lines and 16 columns
#' 
#'  Finally, the previous table is use as input in the writeMetagene() function.
#'  
#'  OUTPUT:
#'  Write a txt file (can be transformed as a .fasta or .fa file) containing each transcript sequence
#'  with their header
#'  
#'  @usage 
#'  write_cdna(ToWrite, "SequenceCdna")

# library
library("biomaRt")

# Mart (show available marts): ## NEEDS INTERNET
listMarts()
listEnsembl(version=90) # 90: Genome GRCh38.90.gff3 used
# we want homo sapiens genes from ensembl database:
ensembl90 = useEnsembl(biomart="ensembl",version=90)

ld = listDatasets(ensembl90)
ld[which(ld[,]=="hsapiens_gene_ensembl"),1]
mart = useDataset("hsapiens_gene_ensembl",mart=ensembl90)

attributes = listAttributes(mart)
attribs = attributes[c(1,3,8,14,15,16,17,21,200,201,202,203, 211,1694,1695,1702,1707),1]

# split the big file into smaller ones (4,5k lines per file):
#' system("split -l 4500 Homo_sapiens.GRCh38.90_mARN_id.txt")

# Parsing the database to retrieve sequences :
split_files = system("ls ../results/split_mRNA_file/x*", intern = TRUE)

# METADATA Table:
# for each transcript id, retrieve gene id and other attributes mentioned in 'attribs'
tmp = NULL
metadata = NULL
for(f in seq(length(split_files))){
  transcript_id = read.table(split_files[f])
  tmp = getBM(
    attributes = attribs,
    filters = 'ensembl_transcript_id', # order by (as in sql)
    values = transcript_id[,1],
    mart = mart
  )
  metadata = rbind(metadata, tmp)
  print(paste("done",f,"/",length(split_files)))
  #print(tmp)
  tmp = NULL
}

#' 1) get all transcrits and their length ## 95 274 transcrits
#' For each ID: retrieve gene_id, transcript_id, sequence length
print(paste0("START get 'trlength' table: gene transcript length"))

tmp = NULL
trlength = NULL
for(f in seq(length(split_files))){
  transcript_id = read.table(split_files[f])
  tmp = getBM(
    attributes = attributes[c(1,3,17),1],
    filters = 'ensembl_transcript_id', # order by (as in sql)
    values = transcript_id[,1],
    mart = mart
  )
  trlength = rbind(trlength, tmp)
  print(paste("done",f,"/",length(split_files)))
  tmp = NULL
}

#' 2) get transcripts that are the longest ## 19 921 transcrits
#' Only keep the longest transcripts (which contains the smaller ones)
#' --> the biggest overlaps the smaller spliced transcripts
#' --> at the end 1 transcript representing 1 gene

lisTranscriptId <- function(df, filename){
  list_id_max = NULL
  id_max = NULL
  for(gid in unique(df[,1])){
    tmp_max = 0
    #print("change gid")
    #print(gid)
    transcrt = unlist(unique(df[which(df[,1]==gid),2]))
    for(t in transcrt){
      if(length(unlist(strsplit(t, ";")))==1){
        #print(t)
        #print(paste0("TMP MAX ",tmp_max))
        tmp = as.numeric(unique(df[which(df[,2]==t),3]))
        #print(tmp)
        if(tmp > tmp_max){
          tmp_max = tmp
          id_max = t
        }
      }
    }
    #list_id_max = append(list_id_max, id_max)
    write(id_max, paste0(filename,".txt"), append = TRUE)
    print(paste0("wrote ", id_max))
    #print(paste0("Hors de for t: tmp_max = ", tmp_max, " id_max = ", id_max))
  }
  #l = as.data.frame(list_id_max)
  return("writting done")
}

# Write the selected IDs in a txt file
lisTranscriptId(trlength, "../results/listeTrIDs")

# split the the selected IDs file into smaller ones (4,5k lines per file):
#' system("split -l 4500 listeTrIDs.txt")
split_listeTrIDs = system("ls x*", intern = TRUE)

# 3) get attributes: 1676 == gene_exon ## exons sequences
print(paste0("START get 'attribToMerge' table: attributes to merge"))

tmp = NULL
attribToMerge = NULL
for(f in seq(length(split_listeTrIDs))){
  transcript_id = read.table(split_listeTrIDs[f])
  tmp = getBM(
    attributes = attributes[c(1676,1,3,1715),1],
    filters = 'ensembl_transcript_id', # order by (as in sql)
    values = transcript_id[,1],
    mart = mart
  )
  attribToMerge = rbind(attribToMerge, tmp)
  print(paste("done",f,"/",length(split_listeTrIDs)))
  tmp = NULL
}

# 4) Merge both table
print(paste0("START get 'ToWrite' table: Merged table"))

ToWrite =  merge.data.frame(metadata, attribToMerge, by = c("ensembl_transcript_id","ensembl_gene_id", "rank"))

# 5) To write: 

writeMetagene <- function(df, filename){
  #' @description 
  #' Write sequence and header into txt or FASTA file from merged table 'ToWrite'
  #' + construct a dataframe with the concatenated sequence and its metadata
  #' 
  #' @usage 
  #' METAGENE = writeMetagene(ToWrite, "../results/SeqMetagene")
  
  metagene = NULL
  for(gid in unique(df[,2])){
    tim = df[which(df[,2]==gid),]
    tam = tim[order(tim[,"rank"]),]
    tmp_head = NULL
    tmp_seq = NULL
    FLAG = 0
    for(ligne in seq(dim(tam)[1])){
      #print(ligne)
      if(tam[ligne,"rank"]==1){
        if(!is.na(tam[ligne,14])){
          FLAG = 1
          #print(ligne)
          #print(tam[ligne,14])
          geneid = tam[ligne,2]; trid = tam[ligne,1]; gename = tam[ligne,9]
          trstart = tam[ligne,5]; trend = tam[ligne,6]; trss = tam[ligne,7]
          p5start = tam[ligne,10]; p5end = tam[ligne,11]; cdnastart = tam[ligne,14]
          tmp_head = paste0(">",tam[ligne,2],"_",tam[ligne,1],"_",tam[ligne,9],"_",
                            tam[ligne,5],"_",tam[ligne,6],"_",tam[ligne,7],"_",
                            tam[ligne,10],"_",tam[ligne,11],"_",tam[ligne,14])
          tmp_seq = paste0(tam[ligne,18])
          #print(paste("rank =",tam[ligne,"rank"],"seqlen",nchar(tam[ligne,18]), sep = " "))
        }else if(is.na(tam[ligne,14])){
          # si NA: juste récup la seq
          #print(paste("rank =",tam[ligne,"rank"],"seqlen",nchar(tam[ligne,18]), sep = " "))
          tmp_seq = paste0(tam[ligne,18])
        }
      }else if(tam[ligne,"rank"]!=1 & !is.na(tam[ligne,14]) & FLAG == 0){
        # changer le flag, prendre les info et la sequence 
        FLAG = 1
        #print(ligne)
        #print(tam[ligne,14])
        geneid = tam[ligne,2]; trid = tam[ligne,1]; gename = tam[ligne,9]
        trstart = tam[ligne,5]; trend = tam[ligne,6]; trss = tam[ligne,7]
        p5start = tam[ligne,10]; p5end = tam[ligne,11]; cdnastart = tam[ligne,14]
        tmp_head = paste0(">",tam[ligne,2],"_",tam[ligne,1],"_",tam[ligne,9],"_",
                          tam[ligne,5],"_",tam[ligne,6],"_",tam[ligne,7],"_",
                          tam[ligne,10],"_",tam[ligne,11],"_",tam[ligne,14])
        tmp_seq = paste0(tmp_seq, tam[ligne,18])
        #print(paste("rank =",tam[ligne,"rank"],"seqlen",nchar(tam[ligne,18]), sep = " "))
      }
      else if(tam[ligne,"rank"]==dim(tam)[1]){
        # prendre la seq et les attributs de fin de seq
        p3start = tam[ligne,12]; p3end = tam[ligne,13]; cdnaend = tam[ligne,15]
        tmp_head = paste(tmp_head, tam[ligne,12],tam[ligne,13],tam[ligne,15], sep = "_")
        tmp_seq = paste0(tmp_seq, tam[ligne,18])
        #print(paste("rank =",tam[ligne,"rank"],"seqlen",nchar(tam[ligne,18]), sep = " "))
      }
      else{
        # pour tous les autres rangs, prendre que la seq
        tmp_seq = paste0(tmp_seq, tam[ligne,18])
        #print(paste("rank =",tam[ligne,"rank"],"seqlen",nchar(tam[ligne,18]), sep = " "))
      }
    }
    write(tmp_head, paste0(filename,".txt"), append = TRUE)
    write(tmp_seq, paste0(filename,".txt"), append = TRUE)
    test = data.frame(tmp_seq, geneid, trid, gename, trstart, trend, trss, p5start,
                      p5end, p3start, p3end, cdnastart, cdnaend)
    metagene = rbind(metagene, test)
  }
  print("done writting")
  return(metagene)
}

#' METAGENE = writeMetagene(ToWrite, "../results/SeqMetagene")

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
#'  METADATA table (same table as in metagene_script.R): the getBM function is used and takes the variable
#'  "attribs" as an input to create a dataframe with every metadata for each of the transcripts.
#'  Since exon features (rank, 5 prime and 3 prime UTRs) are retrieved, the table has this dimension:
#'  843,569 lines and 15 columns
#'  
#'  cDNA table: For each transcript id, the matching cDNA sequence is retrieved using the getSequence
#'  function. cDNA : coding sequence (all exon) + UTRs 
#'  dataframe dimension : 95,274 lines and 2 columns
#'  
#'  MERGED DF CDNA: (variable named "merge.cdna") : Combination of both previous data frames on the
#'  "ensembl_transcript_id" column. 
#'  Creates a new dataframe with 843,569 lines and 16 columns
#'  
#'  OUTPUT:
#'  Write a txt file (can be transformed as a .fasta or .fa file) containing each transcript sequence
#'  with their header
#'  
#'  @usage 
#'  write_cdna(merge.cdna, "SequenceCdna")


# library
library("biomaRt")

# Mart (show available marts): ## NEEDS INTERNET
listMarts()
listEnsembl(version=90) # version actuelle: 91/ 90: version du genome utilisé
# we want homo sapiens genes from ensembl database:
ensembl90 = useEnsembl(biomart="ensembl",version=90)

ld = listDatasets(ensembl90)
ld[which(ld[,]=="hsapiens_gene_ensembl"),1]
mart = useDataset("hsapiens_gene_ensembl",mart=ensembl90)

attributes = listAttributes(mart)
attribs = attributes[c(1,3,8,14,15,16,17,21,200,201,202,203, 211,1694,1695,1702),1]

# split the big file into smaller ones (4,5k lines per file):
#' system("split -l 4500 ../results/10112017/Homo_sapiens.GRCh38.90_mARN_id.txt")

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

# Retrieve cDNA sequences : 1 cDNA per transcript, many transcript per gene
tmp3 = NULL
cdna_table = NULL

for(f in seq(length(split_files))){
  transcript_id = read.table(split_files[f])
  tmp3 = getSequence(id = transcript_id[,1], 
                     type="ensembl_transcript_id",
                     seqType="cdna",
                     mart=mart)
  cdna_table = rbind(cdna_table, tmp3)
  print(paste("done",f,"/",length(split_files)))
  tmp3 = NULL
}


## MERGING TABLES
# Merge by ensembl_transcript_id
merge.cdna = merge.data.frame(metadata, cdna_table, by="ensembl_transcript_id")


#############################################################################
# FUNCTION:
write_cdna <- function(df, filename){
  #' @description 
  #' Write a txt file (can be transformed as a .fasta or .fa file) containing each 
  #' transcript sequence with their header
  #' @usage
  #' write_cdna(merge.cdna, "SequenceCdna")
  
  for(gid in unique(df[,2])){
    #print(gid)
    transcrt = unlist(unique(df[which(df[,2]==gid),1]))
    for(t in transcrt){
      FLAG = 1
      nbtranscrits = length(which(df[,1]==t)) #pour un transcrit du gene
      #print(t)
      for(ligne in which(df[,1]==t)){
        #print(ligne)
        if(df[ligne,2] == gid){
          if(nbtranscrits == 1){
            tmp_md = paste0(">", df[ligne,2], "_", df[ligne,1], "_", 
                            df[ligne,8], "_", df[ligne,4], "_", df[ligne,5], "_",
                            df[ligne,6], "_", df[ligne,7], "_", df[ligne,9], "_",
                            df[ligne,10], "_", df[ligne,11], "_", df[ligne,12], "_",
                            df[ligne,14], "_", df[ligne,15], "_", df[ligne,16])
            #print(tmp_md)
            write(tmp_md, paste0(filename,".txt"), append = TRUE)
            #print(df[ligne, 17])
            write(df[ligne, 17], paste0(filename,".txt"), append = TRUE)
            break
          }
          if(nbtranscrits > 1){
            if(FLAG == 1){
              tmp_md = paste0(">", df[ligne,2], "_", df[ligne,1], "_", 
                              df[ligne,8], "_", df[ligne,4], "_", df[ligne,5], "_",
                              df[ligne,6], "_", df[ligne,7], "_", df[ligne,9], "_",
                              df[ligne,10], "_", df[ligne,11], "_", df[ligne,12], "_",
                              df[ligne,14], "_", df[ligne,15], "_", df[ligne,16])
              FLAG = FLAG + 1
              #print(c("apres if(FLAG == 1): FLAG value = ", FLAG))
              next
            }else if(FLAG < nbtranscrits){
              tmp_md = paste(tmp_md, df[ligne,9], df[ligne,10], 
                             df[ligne,11], df[ligne,12], df[ligne,14],
                             df[ligne,15], df[ligne,16], sep = "_")
              FLAG = FLAG + 1
              #print(c("apres if(FLAG < nbtranscrits): FLAG value = ", FLAG))
              next
            }else if(FLAG == nbtranscrits){
              tmp_md = paste(tmp_md, df[ligne,9], df[ligne,10], 
                             df[ligne,11], df[ligne,12], df[ligne,14],
                             df[ligne,15], df[ligne,16], sep = "_")
              FLAG = FLAG + 1
              #print(c("apres if(FLAG == nbtranscrits): FLAG value = ", FLAG))
              #print(tmp_md)
              #print(df[ligne, 17])
              write(tmp_md, paste0(filename,".txt"), append = TRUE)
              write(df[ligne, 17], paste0(filename,".txt"), append = TRUE)
              next
            }
          }
        }
      }
    }
  }
  return("writting done")
}
#############################################################################
#' Create a dataframe as well as a txt/fasta file VERY LONG (6 hours)
write_cdna <- function(df, filename){
  #' @description 
  #' Write a txt file (can be transformed as a .fasta or .fa file) containing each 
  #' transcript sequence with their header
  #' @usage
  #' write_cdna(merge.cdna, "SequenceCdna")
  
  final = NULL
  for(gid in unique(df[,2])){
    #print(gid)
    transcrt = unlist(unique(df[which(df[,2]==gid),1]))
    for(t in transcrt){
      FLAG = 1
      nbtranscrits = length(which(df[,1]==t)) #pour un transcrit du gene
      #print(t)
      for(ligne in which(df[,1]==t)){
        #print(ligne)
        if(df[ligne,2] == gid){
          if(nbtranscrits == 1){
            geneid = df[ligne,2]; trid = df[ligne,1]; gename = df[ligne,8]
            trstart = df[ligne,4]; trend = df[ligne,5]; trss = df[ligne,6]
            p5start = df[ligne,9]; p5end = df[ligne,10]
            p3start = df[ligne,11]; p3end = df[ligne,12]
            cdnastart = df[ligne,14]; cdnaend = df[ligne,15]
            trtype = df[ligne,16]
            #
            tmp_md = paste0(">", df[ligne,2], "_", df[ligne,1], "_", 
                            df[ligne,8], "_", df[ligne,4], "_", df[ligne,5], "_",
                            df[ligne,6], "_", df[ligne,7], "_", df[ligne,9], "_",
                            df[ligne,10], "_", df[ligne,11], "_", df[ligne,12], "_",
                            df[ligne,14], "_", df[ligne,15], "_", df[ligne,16])
            #print(tmp_md)
            write(tmp_md, paste0(filename,".txt"), append = TRUE)
            #print(df[ligne, 18])
            write(df[ligne, 18], paste0(filename,".txt"), append = TRUE)
            #
            tmp_seq = df[ligne,18]
            test = data.frame(tmp_seq, geneid, trid, gename, trstart, trend, trss, trtype,
                              p5start, p5end, p3start, p3end, cdnastart, cdnaend)
            final = rbind(final, test)
            break
          }
          if(nbtranscrits > 1){
            if(FLAG == 1){
              geneid = df[ligne,2]; trid = df[ligne,1]; gename = df[ligne,8]
              trstart = df[ligne,4]; trend = df[ligne,5]; trss = df[ligne,6]
              p5start = df[ligne,9]; p5end = df[ligne,10]; cdnastart = df[ligne,14]
              trtype = df[ligne,16]
              #
              tmp_md = paste0(">", df[ligne,2], "_", df[ligne,1], "_", 
                              df[ligne,8], "_", df[ligne,4], "_", df[ligne,5], "_",
                              df[ligne,6], "_", df[ligne,9], "_", df[ligne,10], "_", 
                              df[ligne,14], "_", df[ligne,16])
              FLAG = FLAG + 1
              #print(c("apres if(FLAG == 1): FLAG value = ", FLAG))
              next
            }else if(FLAG < nbtranscrits){
              FLAG = FLAG + 1
              #print(c("apres if(FLAG < nbtranscrits): FLAG value = ", FLAG))
              next
            }else if(FLAG == nbtranscrits){
              p3start = df[ligne,11]; p3end = df[ligne,12]; cdnaend = df[ligne,15]
              #
              tmp_md = paste(tmp_md, df[ligne,11], df[ligne,12], df[ligne,15], sep = "_")
              FLAG = FLAG + 1
              #print(c("apres if(FLAG == nbtranscrits): FLAG value = ", FLAG))
              #print(tmp_md)
              #print(df[ligne, 18])
              write(tmp_md, paste0(filename,".txt"), append = TRUE)
              write(df[ligne, 18], paste0(filename,".txt"), append = TRUE)
              #
              tmp_seq = df[ligne,18]
              #
              test = data.frame(tmp_seq, geneid, trid, gename, trstart, trend, trss, trtype,
                                p5start, p5end, p3start, p3end, cdnastart, cdnaend)
              final = rbind(final, test)
              next
            }
          }
        }
      }
    }
  }
  print("writting done")
  return(final)
}

#############################################################################

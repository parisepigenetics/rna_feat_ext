###########################
# Antoine LU              #
# M2BI - Projet Long 2017 #
###########################
# References
# Folding energy calculator:
#http://rna.urmc.rochester.edu/RNAstructureWeb/Servers/AllSub/AllSub.html
#https://www.tbi.univie.ac.at/RNA/tutorial/

#' @description  
#' This script aims to write a FASTA file of the "metagenes" of each gene
#' (total of genes: 19 921)
#' A metagene is the union of all the exon of one gene
#' 

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

# 1) get all transcrits and their length ## 95 274 transcrits
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

# 2) get transcripts that are the longest ## 19 921 transcrits
### Ecriture dans un fichier txt (Rapide)
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

# pas besoin de stocker dans une variable, sera écrit dans le fichier
lisTranscriptId(trlength, "../results/listeTrIDs")

# split the big file into smaller ones (4,5k lines per file):
#' system("split -l 4500 listeTrIDs.txt")
split_listeTrIDs = system("ls x*", intern = TRUE)


# 3) get attributes dont 1676 == gene_exon ## exons
print(paste0("START get 'attribToMerge' table: attributes to merge"))

tmp = NULL
attribToMerge = NULL
for(id in seq(dim(listeTrID)[1])){
  tmp = getBM(
    attributes = attributes[c(1676,1,3,1715),1],
    filters = 'ensembl_transcript_id', # order by (as in sql)
    values = listeTrID[id,1],
    mart = mart
  )
  attribToMerge = rbind(attribToMerge, tmp)
  print(paste("done",id,"/",dim(listeTrID)[1]))
  tmp = NULL
}

### IF READ FROM SPLIT FILE:
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
  
  final = NULL
  for(gid in unique(df[,2])){
    tmp_head = NULL
    tmp_seq = NULL
    #print(gid)
    tim = df[which(df[,2]==gid),]
    tam = tim[order(tim[,"rank"]),]
    for(ligne in seq(dim(tam)[1])){
      #print(ligne)
      if(tam[ligne,"rank"]==1){
        geneid = tam[ligne,2]; trid = tam[ligne,1]; gename = tam[ligne,9]
        trstart = tam[ligne,5]; trend = tam[ligne,6]; trss = tam[ligne,7]
        p5start = tam[ligne,10]; p5end = tam[ligne,11]; cdnastart = tam[ligne,14]
        tmp_head = paste0(">",tam[ligne,2],"_",tam[ligne,1],"_",tam[ligne,9],"_",
                          tam[ligne,5],"_",tam[ligne,6],"_",tam[ligne,7],"_",
                          tam[ligne,10],"_",tam[ligne,11],"_",tam[ligne,14])
        tmp_seq = paste0(tam[ligne,18])
      }
      else if(tam[ligne,"rank"]==dim(tam)[1]){
        p3start = tam[ligne,12]; p3end = tam[ligne,13]; cdnaend = tam[ligne,15]
        tmp_head = paste(tmp_head, tam[ligne,12],tam[ligne,13],tam[ligne,15], sep = "_")
        tmp_seq = paste0(tmp_seq, tam[ligne,18])
      }
      else{
        tmp_seq = paste0(tmp_seq, tam[ligne,18])
      }
    }
    #print(tmp_head)
    #print(tmp_seq)
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

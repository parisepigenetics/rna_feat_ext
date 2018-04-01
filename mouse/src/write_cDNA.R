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
library("methods")
myArgs <- commandArgs(trailingOnly = TRUE)

ENSEMBL_IDs_file = myArgs[1]

# Mart (show available marts): ## NEEDS INTERNET
mmusculus = useMart("ENSEMBL_MART_ENSEMBL",
                    dataset="mmusculus_gene_ensembl")

#ld[which(ld[,]=="mc57bl6nj_gene_ensembl"),1]
#mart = useDataset("mc57bl6nj_gene_ensembl",mart=ensembl)

filters = listFilters(mart)
attributes = listAttributes(mart)
attribs = attributes[c(1,3,13,14,15,161,165,166,167,168,179,180,22,17),1]

# if file = idfile & length(file)<=4500 

metadata = NULL
gene_id = read.table(ENSEMBL_IDs_file)
metadata = getBM(
  attributes = attribs,
  filters = 'ensembl_gene_id', # order by (as in sql)
  values = gene_id[,1],
  mart = mmusculus
)

# idée : récupérer les informations du header avec getBM()
# récupérer les séquences pour les identifiants soit de transcript_id[,1]
# soit de la bonne colonne 'ensembl_transcript_id' de la matrice de getBM()

cdna_table = getSequence(id = metadata[,2], 
                  type="ensembl_transcript_id",
                  seqType="cdna", #(que exon et utrs) ou # gene exon
                  mart=mmusculus)

# merging tables: mRNA
merge.cdna = merge.data.frame(metadata, cdna_table, by="ensembl_transcript_id")

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
            geneid = df[ligne,2]; trid = df[ligne,1]; gename = df[ligne,14]
            trstart = df[ligne,3]; trend = df[ligne,4]; trss = df[ligne,5]
            p5start = df[ligne,7]; p5end = df[ligne,8]
            p3start = df[ligne,9]; p3end = df[ligne,10]
            cdnastart = df[ligne,11]; cdnaend = df[ligne,12]
            trtype = df[ligne,13]
            #
            tmp_md = paste0(">", df[ligne,2], "_", df[ligne,1], "_",
                            df[ligne,14], "_", df[ligne,3], "_", df[ligne,4], "_",
                            df[ligne,5], "_", df[ligne,6], "_", df[ligne,7], "_",
                            df[ligne,8], "_", df[ligne,9], "_", df[ligne,10], "_",
                            df[ligne,11], "_", df[ligne,12], "_", df[ligne,13])
            #print(tmp_md)
            write(tmp_md, paste0(filename,".txt"), append = TRUE)
            #print(df[ligne, 15]) # cdna, print la sequence
            write(df[ligne, 15], paste0(filename,".txt"), append = TRUE)
            #
            tmp_seq = df[ligne,15]
            test = data.frame(tmp_seq, geneid, trid, gename, trstart, trend, trss, trtype,
                              p5start, p5end, p3start, p3end, cdnastart, cdnaend)
            final = rbind(final, test)
            break
          }
          if(nbtranscrits > 1){
            if(FLAG == 1){
              geneid = df[ligne,2]; trid = df[ligne,1]; gename = df[ligne,14]
              trstart = df[ligne,3]; trend = df[ligne,4]; trss = df[ligne,5]
              p5start = df[ligne,7]; p5end = df[ligne,8]; cdnastart = df[ligne,11]
              trtype = df[ligne,13]
              #
              tmp_md = paste0(">", df[ligne,2], "_", df[ligne,1], "_",
                              df[ligne,14], "_", df[ligne,3], "_", df[ligne,4], "_",
                              df[ligne,5], "_", df[ligne,6], "_", df[ligne,7], "_",
                              df[ligne,8], "_", df[ligne,9], "_", df[ligne,10], "_",
                              df[ligne,11], "_", df[ligne,12], "_", df[ligne,13])
              FLAG = FLAG + 1
              #print(c("apres if(FLAG == 1): FLAG value = ", FLAG))
              next
            }else if(FLAG < nbtranscrits){
              tmp_md = paste(tmp_md, df[ligne,7], df[ligne,8], 
                             df[ligne,9], df[ligne,10], df[ligne,11],
                             df[ligne,12], df[ligne,13], sep = "_")
              FLAG = FLAG + 1
              #print(c("apres if(FLAG < nbtranscrits): FLAG value = ", FLAG))
              next
            }else if(FLAG == nbtranscrits){
              p3start = df[ligne,9]; p3end = df[ligne,10]; cdnaend = df[ligne,12]
              #
              tmp_md = paste(tmp_md, df[ligne,7], df[ligne,8], 
                             df[ligne,9], df[ligne,10], df[ligne,11],
                             df[ligne,12], df[ligne,13], sep = "_")
              FLAG = FLAG + 1
              #print(c("apres if(FLAG == nbtranscrits): FLAG value = ", FLAG))
              #print(tmp_md)
              #print(df[ligne, 15])
              write(tmp_md, paste0(filename,".txt"), append = TRUE)
              write(df[ligne, 15], paste0(filename,".txt"), append = TRUE)
              #
              tmp_seq = df[ligne,15]
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
  print("Writting done!")
  return(final)
}

writtenCDNA = write_cdna(merge.cdna, myArgs[2])
write.csv(writtenCDNA, paste0(myArgs[2],".csv"))

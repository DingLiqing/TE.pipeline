#### te.pipeline 2025-10-19
#### TE-Gene list
#### Based on the TE_Gene_position.tsv to te_gene_m.txt(position<40kb)
library(stringr)
library(tidyr)
library(ggplot2)
library(progress)
library(readr)
##TE_gene
te_gene <- read_tsv("TE_Gene_position.tsv")
te_gene <- te_gene[which(te_gene$Position!="None"),]
colnames(te_gene)
gene_id <- as.data.frame(te_gene$Gene_name)
colnames(gene_id) <- 'gene'
a <- paste0("a",c(1:11))
gene_id2 <- separate(gene_id,gene,into =a ,
                     ,sep = '"')
gene_id2 <- gene_id2[,c("a2","a4","a6")]
colnames(gene_id2) <- c('gene_id',"gene_type","gene_name")
te_id <- as.data.frame(te_gene$TE_name)
colnames(te_id ) <- 'te'
a <- paste0("a",c(1:11))
te_id2 <- separate(te_id,te,into =a ,
                   ,sep = '"')
te_id2 <- te_id2[,c("a2","a4","a6","a8","a10")]
colnames(te_id2) <- c('gene_id',"transcript_id","family_id",
                      "class_id","gene_name")
#combind
te_gene2 <- cbind(gene_id2,te_id2,te_gene)
te_gene2 <- te_gene2[,-c(9,13)]
colnames(te_gene2) <- c("Gene_id","Gene_type","Gene_name",    
                        "TE_id","TE_transcript_id","TE_family_id",    
                        "TE_class_id","TE_name",
                        "Gene_start" ,"Gene_end","Gene_strand",
                        "TE_start","TE_end","Chromosome","Position")
write.table(te_gene2,"te_gene_m.txt")


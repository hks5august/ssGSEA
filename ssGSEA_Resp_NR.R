#install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSVAdata")
library(GSVA)
library(GSVAdata)

library(Biobase)

set.seed(7)
setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Responder_vs_Non_Responder/ssGSEA/")

#complete data files
#training <- read.table("Training_23270_Complete_log_data.txt", sep="\t", header=T, row.names = 1, check.names = F)
#Ext_test <- read.table("External_test_23270_Complete_log_data.txt", sep="\t", header=T, row.names = 1, check.names = F)

######## LTS and STS treated samples files ######
training <- read.table("train_res_NR_log_data.txt", sep="\t", header=T, row.names = 1, check.names = F)
Ext_test <- read.table("Ext_test_NR_log_data.txt", sep="\t", header=T, row.names = 1, check.names = F)



mat <- as.matrix(training)
mat_test <- as.matrix(Ext_test)

msigdb_GMTs <- "/Users/kaurh8/Documents/mSigDB_ssGSEA_pathways"
msigdb_kegg <- "c2.cp.kegg.v2022.1.Hs.symbols.gmt"
msigdb_reactome <- "c2.cp.reactome.v2022.1.Hs.symbols.gmt"
msigdb_biocarta <- "c2.cp.biocarta.v2022.1.Hs.symbols.gmt"
msigdb_GO <- "c5.all.v2022.1.Hs.symbols.gmt"
msigdb_oncogenic <- "c6.all.v2022.1.Hs.symbols.gmt"
msigdb_hallmark<- "h.all.v2022.1.Hs.symbols.gmt"
msigdb_immune_sig<- "c7.all.v2022.1.Hs.symbols.gmt"

geneset_reactome <- getGmt(file.path(msigdb_GMTs, msigdb_reactome),  geneIdType=SymbolIdentifier())
geneset_kegg <- getGmt(file.path(msigdb_GMTs, msigdb_kegg),  geneIdType=SymbolIdentifier())
geneset_biocarta <- getGmt(file.path(msigdb_GMTs, msigdb_biocarta),  geneIdType=SymbolIdentifier())
geneset_GO <- getGmt(file.path(msigdb_GMTs, msigdb_GO),  geneIdType=SymbolIdentifier())
geneset_oncogenic <- getGmt(file.path(msigdb_GMTs, msigdb_oncogenic),  geneIdType=SymbolIdentifier())
geneset_hallmark <- getGmt(file.path(msigdb_GMTs, msigdb_hallmark),  geneIdType=SymbolIdentifier())
geneset_immune_sig <- getGmt(file.path(msigdb_GMTs, msigdb_immune_sig),  geneIdType=SymbolIdentifier())

#training mats
ssgsea_train_reactome <- gsva(mat, geneset_reactome, method="ssgsea",  min.sz=5, max.sz=500) #tau=0.25, 
ssgsea_train_kegg <- gsva(mat, geneset_kegg, method="ssgsea", min.sz=5, max.sz=500) #tau=0.25, 
ssgsea_train_GO <- gsva(mat, geneset_GO, method="ssgsea",  min.sz=5, max.sz=500) #tau=0.25, 
ssgsea_train_oncogenic <- gsva(mat, geneset_oncogenic, method="ssgsea",  min.sz=5, max.sz=500) #tau=0.25, 
ssgsea_train_hallmark <- gsva(mat, geneset_hallmark, method="ssgsea",  min.sz=5, max.sz=500) #tau=0.25, 



write.table(cbind("ID"=rownames(ssgsea_train_hallmark), ssgsea_train_hallmark),file="ssgsea_train_hallmark.txt",sep="\t",quote=F, row.names=F)


# save files 
write.table(cbind("ID"=rownames(ssgsea_train_reactome), ssgsea_train_reactome),file="ssgsea_train_reactome.txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(ssgsea_train_kegg), ssgsea_train_kegg),file="ssgsea_train_kegg.txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(ssgsea_train_GO), ssgsea_train_GO),file="ssgsea_train_GO.txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(ssgsea_train_oncogenic ), ssgsea_train_oncogenic ),file="ssgsea_train_oncogenic .txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(ssgsea_train_hallmark), ssgsea_train_hallmark),file="ssgsea_train_hallmark.txt",sep="\t",quote=F, row.names=F)



#test mats
ssgsea_test_reactome <- gsva(mat_test, geneset_reactome, method="ssgsea",   min.sz=5, max.sz=500) #tau=0.25,
ssgsea_test_kegg <- gsva(mat_test, geneset_kegg, method="ssgsea",  min.sz=5, max.sz=500) #tau=0.25,
ssgsea_test_GO <- gsva(mat_test, geneset_GO, method="ssgsea", min.sz=5, max.sz=500) # tau=0.25, 
ssgsea_test_oncogenic <- gsva(mat_test, geneset_oncogenic, method="ssgsea",  min.sz=5, max.sz=500) #tau=0.25,
ssgsea_test_hallmark <- gsva(mat_test, geneset_hallmark, method="ssgsea",  min.sz=5, max.sz=500) #tau=0.25,


# save files 
write.table(cbind("ID"=rownames(ssgsea_test_reactome), ssgsea_test_reactome),file="ssgsea_test_reactome.txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(ssgsea_test_kegg), ssgsea_test_kegg),file="ssgsea_test_kegg.txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(ssgsea_test_GO), ssgsea_test_GO),file="ssgsea_test_GO.txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(ssgsea_test_oncogenic ), ssgsea_test_oncogenic ),file="ssgsea_test_oncogenic .txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(ssgsea_test_hallmark), ssgsea_test_hallmark),file="ssgsea_test_hallmark.txt",sep="\t",quote=F, row.names=F)

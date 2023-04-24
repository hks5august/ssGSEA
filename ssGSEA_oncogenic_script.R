#https://bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html
#install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSVAdata")
library(GSVA)
library(GSVAdata)

library(Biobase)
args <- commandArgs(TRUE)
set.seed(7)


#path wher have input data all files and folders
path <- paste0("/Users/kaurh8/Documents/GDC_TCGA_Biolinks/GDC_All_samples_Data/updated_folder2/ssGSEA_oncogenic/")
#setwd("/Users/kaurh8/Documents/GDC_TCGA_Biolinks/GDC_All_samples_Data/updated_folder2/ssGSEA_oncogenic/")

#set working directory
setwd(path)
getwd()


#provide list of cancer type
cancer_list <- read.table(paste0(path, "Cancer_list2"),header =T, sep = "\t",  row.names=1,  check.names = FALSE)
cancer_list <- read.table(paste0(path, args[1]),header =T, sep = "\t",  row.names=1,  check.names = FALSE)



folder <- as.character(row.names(cancer_list))
folder
#print("folder:" folder)


for (c in 1:length(folder))
{
  setwd(paste0(path,folder[c], "/"))
  getwd()
  
  path_c <- paste0(path,folder[c], "/")
  #print("cancer folder", path_c)
  path_c
  setwd(path_c)
  getwd()
  
#complete data files
######## input log files ###### Note: Samples in column and genes in rows
#data <- read.table("ACC/ACC_Exp_log_data.tsv", sep="\t", header=T, row.names = 1, check.names = F)
#data <- read.table("../../TCGA-ACC/TCGA-ACC_Matrix_FPKM_log_data_with_unique_sample_ID_without_SD_tpose_with_gene_annotation_only_prot_lncrna_mirna_only_symbol_with_unique_gene_ID_with_SD.txt", sep="\t", header=T, row.names = 1, check.names = F)



data<-read.table(paste0(path_c,folder[c], "_Exp_log_data.tsv"),header =TRUE, sep = "\t", row.names=1,  check.names = FALSE)
dim(data)



mat <- as.matrix(data)
mat <- data.matrix(data)

msigdb_GMTs <- "/Users/kaurh8/Documents/mSigDB_ssGSEA_pathways"
msigdb_oncogenic <- "c6.all.v2022.1.Hs.symbols.gmt"
msigdb_hallmark<- "h.all.v2022.1.Hs.symbols.gmt"

geneset_oncogenic <- getGmt(file.path(msigdb_GMTs, msigdb_oncogenic),  geneIdType=SymbolIdentifier())

geneset_hallmark <- getGmt(file.path(msigdb_GMTs, msigdb_hallmark),  geneIdType=SymbolIdentifier())

geneset_hallmark 

#training mats
#ssgsea_data_oncogenic <- gsva(mat, geneset_oncogenic, method="ssgsea",  min.sz=5, max.sz=500) #tau=0.25, 
#ssgsea_data_oncogenic_1 <- gsva(mat, geneset_oncogenic, method="gsva",  min.sz=5, max.sz=500, kcdf="Gaussian", mx.diff=TRUE, abs.ranking=TRUE, tau=0.25) #tau=0.25, 

#head(ssgsea_data_oncogenic,5)

sgsea_data_hallmark <- gsva(mat, geneset_hallmark, method="ssgsea",  min.sz=5, max.sz=500) #tau=1 

sgsea_data_hallmark <- round(sgsea_data_hallmark , 3)

head(sgsea_data_hallmark,5)
#ssgsea_data_hallmark_1 <- gsva(mat, geneset_hallmark, method="gsva",  min.sz=5, max.sz=500, kcdf="Gaussian", mx.diff=TRUE, abs.ranking=TRUE, tau=0.25) #tau=0.25, 




# save files 
#write.table(cbind("ID"=rownames(ssgsea_data_oncogenic), ssgsea_data_oncogenic),file="ssgsea_oncogenic.tsv",sep="\t",quote=F, row.names=F)

#write.table(cbind("ID"=rownames(sgsea_data_hallmark), sgsea_data_hallmark),file="sgsea_hallmark.tsv",sep="\t",quote=F, row.names=F)


#write.table(cbind("ID"=rownames(ssgsea_data_oncogenic), ssgsea_data_oncogenic),file=paste0(path,folder[c], "/", folder[c], "_", "ssgsea_oncogenic.tsv"),sep="\t",quote=F, row.names=F)


write.table(cbind("ID"=rownames(sgsea_data_hallmark), sgsea_data_hallmark),file=paste0(path,folder[c], "/", folder[c], "_", "ssgsea_hallmark.tsv"),sep="\t",quote=F, row.names=F)


}

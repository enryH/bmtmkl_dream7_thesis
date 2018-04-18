###################################################################################################
# Author: Henry Webel 
# Content: Extract pathways (set of genes) from databases and select information of 
# corresponding genes in data-sets (Kind of subsampling of datasets)
#
# Thanks to Manuela Benary  (manuela.benary@charite.de) on using databases!
###################################################################################################
# DZA:
  setwd("C:/Users/Webel/Desktop/Versuche/statisticsmaster/cancerprototyp/")
# Win-Laptop:
  setwd("C:/Users/Henry/Documents/1_Statistik/statisticsmaster/cancerprototyp/")
# Linux-Laptop:
  setwd("/home/enryh/Desktop/statisticsmaster/cancerprototyp/")
###################################################################################################
load("traindata.RData")  ### Gene Names as column names
###################################################################################################
## Install:
# source("https://bioconductor.org/biocLite.R")
# biocLite("reactome.db")
library(reactome.db)  ## for Version of R 3.4.2 (changes in database along versions of R!!!)
###################################################################################################
## Install: 
# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
###################################################################################################
###################################################################################################
# ### to test
# file <- "DREAM7_DrugSensitivity1_Exomeseq.txt"
# data <- traindata$ExomeSeqRealVal
# df.reactome <- as.list(reactomePATHID2EXTID)
# pathways.reactome <- names(df.reactome)
# pathways.reactome <- pathways.reactome[grep("R-HSA-", pathways.reactome)]
# pathways.reactome <- pathways.reactome[1:2119]
###############################################################################
loadDream7DataPathwaysReactome <- function(file=NULL, data= NULL, fct=max)
{  
  # assign new filename:
  file <- paste0(strsplit(file, ".txt"),".pathway.Reactome")
  cat("newfilename: ", file)
  # get genes in data
  data <- t(data)
  genes <- rownames(data)
  ##################################################### 
  # Extract pathways from "data":
  extracted.pathway.data <- NULL
  j=0
  for (i in pathways.reactome){
    j= j+ 1
    cat(j,". pathway name: ", i,"\n")
    match.gene   <- select(org.Hs.eg.db, keys = df.reactome[i][[i]] , columns =c("SYMBOL", "GENENAME"), keytype = "ENTREZID")
    genes.path   <- match.gene$SYMBOL[(duplicated(match.gene)==FALSE)]
    select       <- genes[c(which(genes %in% genes.path))]  # select genes of path
    pathway_data <- data[select, ]
    if (is.matrix(pathway_data)==TRUE) {
      pathway_data <- apply((pathway_data), MARGIN= 2, FUN=fct)  # max, min, median, mean
    }  # max pooling over genes
    extracted.pathway.data <- rbind(extracted.pathway.data, pathway_data)
    rownames(extracted.pathway.data)[j] <- i
  }
  # save and return results:
  file <- paste0(file,".RData")
  save(extracted.pathway.data, file = file)
  return(extracted.pathway.data)
}
###################################################################################################
# ### to test
# file <- "DREAM7_DrugSensitivity1_Exomeseq.txt"
# data <- traindata$ExomeSeqRealVal
# k <- keys(org.Hs.eg.db,keytype="PATH")  ## 229 Paths
### 
loadDream7DataPathwaysOrgHsEg <- function (file= NULL, data = NULL, fct=max) {
  # assign new filename:
  file <- paste0(strsplit(file, ".txt"),".pathway.OrgHsEg")
  cat("newfilename: ", file)
  # get genes in data
  data <- t(data)
  genes <- rownames(data)
  ##################################################### 
  # Extract pathways from "data":
  extracted.pathway.data <- NULL
  for (i in (1:length(k))) {
    cat(i,". path which identifier: ", k[i])
    path_temp <- select(org.Hs.eg.db, keys=k[i], columns=c("SYMBOL","GENENAME"), keytype="PATH")
    genes.path <- path_temp$SYMBOL
    select <- genes[c(which(genes %in% genes.path))]  # select genes of path
    pathway_data <-  data[select, ]
    if (is.matrix(pathway_data)==TRUE) {
      pathway_data <- apply((pathway_data), MARGIN= 2, FUN=fct)  # max, min, median, mean
    }  # max pooling over genes
    extracted.pathway.data <- rbind(extracted.pathway.data, pathway_data)
    rownames(extracted.pathway.data)[i] <- k[i]
  }
  # save and return results:
  file <- paste0(file,".RData")
  save(extracted.pathway.data, file = file)
  return(extracted.pathway.data)
}
###################################################################################################
###################################################################################################
# load Org.Hs.eg.db Pathways
k <- keys(org.Hs.eg.db,keytype="PATH")  ## 229 Paths
# load Reactome Dataset
df.reactome <- as.list(reactomePATHID2EXTID)
pathways.reactome <- names(df.reactome)
pathways.reactome <- pathways.reactome[grep("R-HSA-", pathways.reactome)]
pathways.reactome <- pathways.reactome[1:2119]
###################################################################################################
###################################################################################################
file <- "DREAM7_DrugSensitivity1_Exomeseq.txt"
dream7.ExomeSeq.pathway.Reactome <- loadDream7DataPathwaysReactome(file, traindata$ExomeSeqRealVal)
dream7.ExomeSeq.pathway.OrgHsEg <-loadDream7DataPathwaysOrgHsEg(file, traindata$ExomeSeqRealVal)
###################################################################################################
file <- "DREAM7_DrugSensitivity1_GeneExpression.txt"
dream7.GeneExpression.pathway.Reactome <- loadDream7DataPathwaysReactome(file, traindata$GeneExpression, mean)
dream7.GeneExpression.pathway.OrgHsEg <-loadDream7DataPathwaysOrgHsEg(file, traindata$GeneExpression, mean)
###################################################################################################
###################################################################################################
file <- "DREAM7_DrugSensitivity1_SNP6_gene_level.txt"
dream7.CNV.pathway.Reactome <- loadDream7DataPathwaysReactome(file, traindata$CNV)
dream7.CNV.pathway.OrgHsEg <-loadDream7DataPathwaysOrgHsEg(file, traindata$CNV)
###################################################################################################
###################################################################################################
file <- "DREAM7_DrugSensitivity1_Methylation.txt"
dream7.Methylation.pathway.Reactome <- loadDream7DataPathwaysReactome(file, traindata$Methylation)
dream7.Methylation.pathway.Reactome  <-loadDream7DataPathwaysOrgHsEg(file, traindata$Methylation)
###################################################################################################
###################################################################################################
file <- "DREAM7_DrugSensitivity1_RNAseq.txt"
dream7.RNAseq.pathway.Reactome <- loadDream7DataPathwaysReactome(file, traindata$RNASeq_Exp)
dream7.RNAseq.pathway.OrgHsEg <-loadDream7DataPathwaysOrgHsEg(file, traindata$RNASeq_Exp)
###################################################################################################
save.image(file="data/pathways_raw.Rdata")
###################################################################################################
# currentwd <- getwd()
# setwd("C:/Users/Webel/Desktop/Versuche/statisticsmaster/cancerprototyp/data/pathways/")
# filelist <- dir()
# pathways = vector(mode='list',length=10)
# for (f in 1:length(filelist)) {
#   print(filelist[f])
#   pathways[[f]] <- load(filelist[f])
# }
# names(pathways)=c('GeneExpression_Reactome','CNV_Reactome','Methylation_Reactome', 
#                   'RNASeq_Reactome','ExomeSeqRealVal_Reactome',
#                   'GeneExpression_OrgHsEg','CNV_OrgHsEg','Methylation_OrgHsEg', 'RNASeq_OrgHsEg','ExomeSeqRealVal_OrgHsEg')
# 
pathways = vector(mode='list',length=10)

pathways[[1]] <- t(dream7.GeneExpression.pathway.Reactome)
pathways[[2]] <- t(dream7.GeneExpression.pathway.OrgHsEg)
pathways[[3]] <- t(dream7.Methylation.pathway.Reactome)
pathways[[4]] <- t(dream7.Methylation.pathway.OrgHsEg)
pathways[[5]] <- t(dream7.CNV.pathway.Reactome)
pathways[[6]] <- t(dream7.CNV.pathway.OrgHsEg)
pathways[[7]] <- t(dream7.ExomeSeq.pathway.Reactome)
pathways[[8]]<- t(dream7.ExomeSeq.pathway.OrgHsEg)
pathways[[9]] <- t(dream7.RNAseq.pathway.Reactome)
pathways[[10]] <- t(dream7.RNAseq.pathway.OrgHsEg)
names(pathways)=c('GeneExpression_Reactome','GeneExpression_OrgHsEg',
                  'Methylation_Reactome', 'Methylation_OrgHsEg',
                  'CNV_Reactome','CNV_OrgHsEg',
                  'ExomeSeqRealVal_Reactome','ExomeSeqRealVal_OrgHsEg',
                  'RNASeq_Reactome','RNASeq_OrgHsEg')
save(pathways, file="pathways.RData")
###################################################################################################
# #################### automatisches Einlesen der Dateien... (spielerei)
# DREAM7_DrugSensitivity1_SNP6_README
# test <- strsplit("aDREAM7_DrugSensitivity1_SNP6_README","DREAM7_DrugSensitivity1_")[1]
# sapply(dir(), function(a) (strsplit(a,"DREAM7_DrugSensitivity1_"))[[1]][2])
# ###
# ##################### Dokumentation
# keytypes(org.Hs.eg.db)  ## available Key-Types 
# head(keys(org.Hs.eg.db, keytype="SYMBOL"))  # Gene-Name (HGNEID)

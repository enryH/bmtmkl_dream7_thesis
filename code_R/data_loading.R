# Author:  	Ammad-ud-din Muhammad, https://users.ics.aalto.fi/ammad/
#
#' Load DREAM7 challenge data sources
#'
#' \code{loadDream7Data} returns the list object containing DREAM7 challenge data sources.
#' The data sources contains genomic and molecular profiles (i.e, .Gene Expression, 
#' Copy Number Variations, Methylations, Reverse Protein Lysate Array (RPPA), 
#' RNA-Seq Expression and Exome-Seq Mutations) of 53 breast cancer cell lines
#' and their sensitivity measurments (i.e., GI50 values) on 31 anti-cancer drugs.
#' @param dataLocation a character string containing the name and location of the 
#' folder where data sources have been stored. Note , the directory must contain 
#' the following files (downlaodable from DREAM7 challenge website)
#' 1. DREAM7_DrugSensitivity1_Predictions.csv
#' 2. DREAM7_DrugSensitivity1_Drug_Response_Training.txt
#' 3. DREAM7_DrugSensitivity1_GeneExpression.txt
#' 4. DREAM7_DrugSensitivity1_SNP6_gene_level.txt 
#' 5. DREAM7_DrugSensitivity1_Methylation.txt
#' 6. DREAM7_DrugSensitivity1_RNAseq_expressed_calls.txt
#' 7. DREAM7_DrugSensitivity1_RNAseq_quantification.txt
#' 8. DREAM7_DrugSensitivity1_RPPA.txt
#' 9. DREAM7_DrugSensitivity1_Exomeseq.txt.
#' @return Y A list containing DREAM7 challenge data sources.
#' @examples
#' dataLocation='DrugSensitivity1/'
#' Y=loadDream7Data(dataLocation)
loadDream7Data <- function(dataLocation=NULL)
{
  # list all the files in the dataLocation folder
  datafiles=list.files(dataLocation)
  
  #Read DREAM7_DrugSensitivity1_Predictions.csv file to get 
  # the names of all cell lines (training + test)
  predictionFile=datafiles[grep('_Predictions.csv',datafiles)]  #select-file-name
  predictions = read.table(paste0(dataLocation,predictionFile),header=T, sep=',', row.names=1)
  allCelllineNames <- rownames(predictions) #Names of all cell lines
  
  #Load the training drug sensitivity/response data
  drugResponseFile=datafiles[grep('_Drug_Response_Training.txt',datafiles)]
  drugResponse=as.matrix(read.table(paste0(dataLocation,drugResponseFile), header=T, row.names=1))  
  trainingCelllineNames=rownames(drugResponse)
  
  #Order the cell lines such that test cell lines are added 
  #at the end of the training cell lines. The length of training cell lines
  #is 1:35 and 36:53 are the test cell lines. 
  allCelllineNames <- c(allCelllineNames[c(which(allCelllineNames %in% trainingCelllineNames))], 
                        allCelllineNames[-c(which(allCelllineNames %in% trainingCelllineNames))])
  
  #Load Gene Expression data 
  dream7.celllineGeneExpression <- t(as.matrix(read.table(paste0(dataLocation,
                                  datafiles[grep('_GeneExpression.txt',datafiles)]), header=T, row.names=1)))
  #Preprocess gene expression data such that the order of the cell lines 
  #matches the order in the "allCelllineNames" variable and for missing cell lines
  #add NAs.
  dream7.celllineGeneExpression <- preprocessData(dream7.celllineGeneExpression, allCelllineNames)
  
  
  #Load Copy Number Variations data 
  dream7.celllineCNV <- read.table(paste0(dataLocation,datafiles[grep('_SNP6_gene_level.txt',datafiles)]))
  entrezIds_cnv <- dream7.celllineCNV[,1]
  hgncIds_cnv   <- dream7.celllineCNV[,2]
  CNV.hgnc2entrezIds_map <- as.data.frame(cbind(as.character(hgncIds_cnv), as.character(entrezIds_cnv)))
  colnames(CNV.hgnc2entrezIds_map) <- c("HGNC_ID", "EntrezID")
  dream7.celllineCNV <- dream7.celllineCNV[,c(-1,-2)]  # without IDs
  cellNames.cnv=as.vector(as.matrix(dream7.celllineCNV[1,]))
  colnames(dream7.celllineCNV)=cellNames.cnv
  dream7.celllineCNV =  dream7.celllineCNV[-1,]
  dream7.celllineCNV <- t(as.matrix(dream7.celllineCNV))
  colnames(dream7.celllineCNV) <- hgncIds_cnv[2:length(hgncIds_cnv)]
  class(dream7.celllineCNV) = 'double'
  #Preprocess copy number variations data such that the order of the cell lines 
  #matches the order in the "allCelllineNames" variable and for missing cell lines
  #add NAs.
  dream7.celllineCNV <- preprocessData(dream7.celllineCNV, allCelllineNames)
  
  #Load Methylation data
  dream7.celllineMethylation <- read.table(paste0(dataLocation,datafiles[grep('_Methylation.txt',datafiles)]),header=T, row.names=1)
  dream7.celllineMethylation <- t(as.matrix(dream7.celllineMethylation))
  dream7.celllineMethylation_annot <- dream7.celllineMethylation[1:3,]
  dream7.celllineMethylation <- dream7.celllineMethylation[-c(1:3),]
  class(dream7.celllineMethylation) <- "numeric"
  colnames(dream7.celllineMethylation) <- dream7.celllineMethylation_annot[1,]
  #Preprocess methylation data such that the order of the cell lines 
  #matches the order in the "allCelllineNames" variable and for missing cell lines
  #add NAs.
  dream7.celllineMethylation <- preprocessData(dream7.celllineMethylation, allCelllineNames)
  
  
  #Load RNA-Seq Expression (real-valued) data
  dream7.celllineRNASeq_exp <- read.table(paste0(dataLocation,datafiles[grep('_RNAseq_expressed_calls.txt',datafiles)]),header=T)
  ensembleIds_exp <- dream7.celllineRNASeq_exp[,2]
  dream7.celllineRNASeq_exp <- dream7.celllineRNASeq_exp[,-2]
  hgncIds_exp <- dream7.celllineRNASeq_exp[,1]
  dream7.celllineRNASeq_exp <- dream7.celllineRNASeq_exp[,-1]
  RNAexp.hgnc2ensembleIds_map <- as.data.frame(cbind(as.character(hgncIds_exp), as.character(ensembleIds_exp)))
  colnames(RNAexp.hgnc2ensembleIds_map) <- c("HGNC_ID", "EnsembleID")
  dream7.celllineRNASeq_exp <- t(as.matrix(dream7.celllineRNASeq_exp))
  colnames(dream7.celllineRNASeq_exp) <- hgncIds_exp
  #Preprocess RNA-seq data such that the order of the cell lines 
  #matches the order in the "allCelllineNames" variable and for missing cell lines
  #add NAs.
  dream7.celllineRNASeq_exp <- preprocessData(dream7.celllineRNASeq_exp, allCelllineNames)
  
  
  #Load RNA-Seq Expression (binary-valued) data
  dream7.celllineRNASeq_quant <- read.table(paste0(dataLocation,datafiles[grep('_RNAseq_quantification.txt',datafiles)]),header=T)
  ensembleIds_quant <- dream7.celllineRNASeq_quant[,2]
  dream7.celllineRNASeq_quant <- dream7.celllineRNASeq_quant[,-2]
  hgncIds_quant <- dream7.celllineRNASeq_quant[,1]
  dream7.celllineRNASeq_quant <- dream7.celllineRNASeq_quant[,-1]
  RNAexp.hgnc2ensembleIds_map <- as.data.frame(cbind(as.character(hgncIds_exp), as.character(ensembleIds_exp)))
  colnames(RNAexp.hgnc2ensembleIds_map) <- c("HGNC_ID", "EnsembleID")
  dream7.celllineRNASeq_quant <- t(as.matrix(dream7.celllineRNASeq_quant))
  colnames(dream7.celllineRNASeq_quant) <- hgncIds_exp
  #Preprocess RNA-seq data such that the order of the cell lines 
  #matches the order in the "allCelllineNames" variable and for missing cell lines
  #add NAs.
  dream7.celllineRNASeq_quant <- preprocessData(dream7.celllineRNASeq_quant, allCelllineNames)
  
  
  #Load RPPA data
  dream7.celllineRPPA <- read.table(paste0(dataLocation,datafiles[grep('_RPPA.txt',datafiles)]))
  colnames(dream7.celllineRPPA)=as.vector(as.matrix(dream7.celllineRPPA[1,]))
  dream7.celllineRPPA=dream7.celllineRPPA[-1,]
  dream7.celllineRPPA <- dream7.celllineRPPA[which(dream7.celllineRPPA$FullyValidated=='Yes'),]
  colNames <- apply(dream7.celllineRPPA,1,function(x) paste(x[1],x[2], sep=":"))
  dream7.celllineRPPA <- dream7.celllineRPPA[,c(-1,-2)]
  dream7.celllineRPPA <- t(as.matrix(dream7.celllineRPPA))
  colnames(dream7.celllineRPPA) <- colNames
  class(dream7.celllineRPPA) = 'double'
  #Preprocess RPPA such that the order of the cell lines 
  #matches the order in the "allCelllineNames" variable and for missing cell lines
  #add NAs.
  dream7.celllineRPPA <- preprocessData(dream7.celllineRPPA, allCelllineNames)
  
  #Preprocess training drug sensitivity/response data such that the 
  #order of the cell lines matches the order in the "allCelllineNames" 
  #variable and for missing cell lines add NAs.
  dream7.drugResponse <- preprocessData(drugResponse, allCelllineNames)
  dream7.drugResponse=as.matrix(dream7.drugResponse)
  
  #Load and pre-process Exome-Seq data
  # dream7.celllineExomeSeq = read.delim('DrugSensitivity1/DREAM7_DrugSensitivity1_Exomeseq.txt')
  dream7.celllineExomeSeq = read.delim('H:/WS1617/bioinformatics_topic/data/NCI_dream_challenge/task1/training/DrugSensitivity1/DREAM7_DrugSensitivity1_Exomeseq.txt')
  dream7.celllineExomeSeq.filtered <- dream7.celllineExomeSeq[(dream7.celllineExomeSeq$Avg.Mismatch.alt. <= 0.5 & dream7.celllineExomeSeq$MismatchQualitySum.alt. <= 7.5), ]
  dream7.celllineExomeSeq <- processExomeSeqDataset(dream7.celllineExomeSeq.filtered)
  dream7.celllineExomeSeq[dream7.celllineExomeSeq < 0.5 ] = 0 
  
  #Convert Exome-Seq real-valued to binary valued data (heuristic)
  dream7.celllineExomeSeq.binary <- dream7.celllineExomeSeq
  dream7.celllineExomeSeq.binary[dream7.celllineExomeSeq.binary >= 0.5 ] = 1
  
  #Preprocess Exome-seq data such that the order of the cell lines 
  #matches the order in the "allCelllineNames" variable and for missing cell lines
  #add NAs.
  dream7.celllineExomeSeq <- preprocessData(dream7.celllineExomeSeq, allCelllineNames)
  dream7.celllineExomeSeq.binary <- preprocessData(dream7.celllineExomeSeq.binary, allCelllineNames)  
  
  
  
  #dream7.celllineExomeSeq.hethom.binary <- processExomeSeqDatasetIntoBinaryMat(dream7.celllineExomeSeq.filtered)
  #dream7.celllineExomeSeq.hethom.binary[1:35,1:20]
  #class(dream7.celllineExomeSeq.matrix.hethom.binary) <- 'numeric'
  #dream7.celllineExomeSeq.hethom.binary <- preprocessData(dream7.celllineExomeSeq.matrix.hethom.binary, allCelllineNames)
  #dream7.celllineExomeSeq[1:5,1:50]
  #dream7.celllineExomeSeq.binary[1:5,1:50]
  #dream7.celllineExomeSeq.hethom.binary[1:5,1:20]
  
  Y = vector(mode='list',length=9)
  
  Y[[1]] <- dream7.drugResponse
  Y[[2]] <- dream7.celllineGeneExpression
  Y[[3]] <- dream7.celllineCNV
  Y[[4]] <- dream7.celllineMethylation
  Y[[5]] <- dream7.celllineRNASeq_exp
  Y[[6]] <- dream7.celllineRNASeq_quant
  Y[[7]] <- dream7.celllineExomeSeq
  Y[[8]] <- dream7.celllineExomeSeq.binary
  Y[[9]] <- dream7.celllineRPPA
  names(Y)=c('DrugResponse','GeneExpression','CNV','Methylation',
             'RNASeq_Quant','RNASeq_Exp','ExomeSeqRealVal','ExomeSeqBinary',
             'RPPA')
  return(Y=Y)
}

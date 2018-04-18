#' @title A function for processing the DREAM7 Exome-Seq data source
#'
#' \code{processExomeSeqDataset} is a utility function and is only called
#' from \code{\link{loadDream7Data}} to preprocess the exome-seq data.
#'
#' @param exomeSeqDataset the exome-seq data source.
#' @return exomeSeqDataset.processed exome-seq data processed based on 
#' heuristic choices (not recommended to use)
processExomeSeqDataset <- function(exomeSeqDataset)
{
  uniqueCelllineNames_ExomeSeq <- unique(exomeSeqDataset[,1])
  count <- 0
  exomeStatList <- list()
  geneName <- (exomeSeqDataset[which(uniqueCelllineNames_ExomeSeq[1] == exomeSeqDataset[,1]),7])
  for(i in 2:length(uniqueCelllineNames_ExomeSeq))
  {
    geneName <- union(geneName, (exomeSeqDataset[which(uniqueCelllineNames_ExomeSeq[i] == exomeSeqDataset[,1]),7]))
  }
  
  uniqueGeneNames <- unique(geneName)
  exomeSeqDataset.processed <- matrix(0, nrow=length(uniqueCelllineNames_ExomeSeq),ncol=length(uniqueGeneNames))
  rownames(exomeSeqDataset.processed) <- uniqueCelllineNames_ExomeSeq
  colnames(exomeSeqDataset.processed) <- uniqueGeneNames
  class(exomeSeqDataset.processed) <- "numeric"
  
  for(k in 1:length(uniqueCelllineNames_ExomeSeq))
  {
    dataSubSet <- (exomeSeqDataset[which(uniqueCelllineNames_ExomeSeq[k] == exomeSeqDataset[,1]),])
    length(dataSubSet[,7] %in% uniqueGeneNames)
    alt.ref.count <- vector('numeric',length=length(dataSubSet[,7]))
    names(alt.ref.count) <- dataSubSet[,7]
    
    for(i in 1:dim(dataSubSet)[1])
    {
      alt.ref.count[i] <- (dataSubSet[i,20]/(dataSubSet[i,20]+dataSubSet[i,19]))
    }
    
    duplicateGenes <- unique(names(alt.ref.count)[duplicated(names(alt.ref.count))])
    for(d in 1:length(duplicateGenes))
    {
      index <- which(duplicateGenes[d] == names(alt.ref.count))
      maxValue <- alt.ref.count[index[1]]
      for(a in 2:length(index))
      {
        if(alt.ref.count[index[a]] >= maxValue)
        {
          maxValue <- alt.ref.count[index[a]]
          alt.ref.count[index[a-1]] <- 0
        }
        else
        {
          alt.ref.count[index[a]] <- 0
        }
      }
    }
    maxValue2 <- max(alt.ref.count)
    if(maxValue2 < 0.5)
    {
      alt.ref.count[which(maxValue2==alt.ref.count)] <- 0
    }
    
    exomeSeqDataset.processed[which(as.character(uniqueCelllineNames_ExomeSeq[k]) == rownames(exomeSeqDataset.processed)),names(alt.ref.count)] <- alt.ref.count
    
  }
  
  return(exomeSeqDataset.processed=exomeSeqDataset.processed)
  
}
# #' @title A function for processing the DREAM7 data sources
# #'
# #' \code{processExomeSeqDatasetIntoBinaryMat} is a utility function and is only called
# #' from \code{\link{loadDream7Data}}.
# #'
# #' @param exomeSeqDataset the exome-seq data source.
# #' @return exomeSeqDataset.binary exome-seq real-valued data converted to binary data based on 
# #' heuristic choices (not recommended to use)
# processExomeSeqDatasetIntoBinaryMat <- function(exomeSeqDataset=NULL)
# {
#   uniqueCelllineNames_ExomeSeq <- unique(exomeSeqDataset[,1])
#   count <- 0
#   exomeStatList <- list()
#   geneName <- (exomeSeqDataset[which(uniqueCelllineNames_ExomeSeq[1] == exomeSeqDataset[,1]),7])
#   for(i in 2:length(uniqueCelllineNames_ExomeSeq))
#   {
#     geneName <- union(geneName, (exomeSeqDataset[which(uniqueCelllineNames_ExomeSeq[i] == exomeSeqDataset[,1]),7]))
#     length((geneName)) 
#   }
#   uniqueGeneNames <- unique(geneName)
#   exomeSeqDataset.binary <- matrix(0, nrow=length(uniqueCelllineNames_ExomeSeq),ncol=length(uniqueGeneNames))
#   rownames(exomeSeqDataset.binary) <- uniqueCelllineNames_ExomeSeq
#   colnames(exomeSeqDataset.binary) <- uniqueGeneNames
#   class(exomeSeqDataset.binary) <- "numeric"
#   
#   for(k in 1:length(uniqueCelllineNames_ExomeSeq))
#   {
#     dataSubSet <- (exomeSeqDataset[which(uniqueCelllineNames_ExomeSeq[k] == exomeSeqDataset[,1]),])
#     het.hom.binary <- vector('character',length=length(dataSubSet[,2]))
#     names(het.hom.binary) <- as.character(dataSubSet$HGNC_ID)
#     het.hom.binary <- as.character(dataSubSet$Zygosity.tumor.)
#     het.hom.binary[which("het"== het.hom.binary)] <- 0
#     het.hom.binary[which("hom"== het.hom.binary)] <- 1
#     names(het.hom.binary) <- as.character(dataSubSet$HGNC_ID)
#     exomeSeqDataset.binary[which(as.character(uniqueCelllineNames_ExomeSeq[k]) == rownames(exomeSeqDataset.binary)),names(het.hom.binary)] <- het.hom.binary
#   }
#   
#   
#   return(exomeSeqDataset.binary=exomeSeqDataset.binary)
#   
# }
###################################################################################################
###################################################################################################
###################################################################################################
#' @title A function for processing the DREAM7 data sources
#'
#' \code{preprocessData} is a utility function and is only called
#' from \code{\link{loadDream7Data}}.
#'
#' @param inputData the input data source.
#' @param allCelllineNames cell line names such that 
#' 1:35 are training cell lines and 36:53 are test cell lines.
#' @return processedInputData processed data such that the order of the cell lines 
#' matches the order in the "allCelllineNames" variable and 
#' for missing cell lines, NA's are added
preprocessData <- function(inputData, allCelllineNames)
{
  
  processedInputData <- matrix(nrow=length(allCelllineNames), ncol=dim(inputData)[2])
  dim(processedInputData)
  rownames(processedInputData) <- allCelllineNames
  colnames(processedInputData) <- colnames(inputData)
  for(i in 1:length(allCelllineNames))
  {
    if(allCelllineNames[i] %in% rownames(inputData))
    {
      processedInputData[i,1:dim(inputData)[2]] <- inputData[which(allCelllineNames[i]==rownames(inputData)),1:dim(inputData)[2]]
    }
    else
    {
      processedInputData[i,] <- NA
    }
  }
  return(processedInputData=processedInputData)
  
}
###################################################################################################
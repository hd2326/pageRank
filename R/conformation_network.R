#' Build Network from Conformation Peaks.
#'
#' Build network from conformation, e.g. HiChIP records.
#'
#' @param table (data.frame) Records, with "Chr1", "Position1", "Strand1", "Chr2", "Position2" and "Strand2" in column name, and record ID in row names.
#'
#' @param promoter (GRanges) Promoter regions.
#' 
#' @param pfm (PFMatrixList) Positon Frequency Matrices (PFMs) of regulators.
#' 
#' @param genome (BSgenome or character) Genome build in which regulator motifs will be searched.
#'
#' @param range (numeric) Search radius from "Position1" and "Position2" for promoters.
#'
#' @param p.cutoff (numeric) P-value cutoff for motifs searching within peaks for TF identificaton.
#'
#' @param w (numeric) Window size for motifs searching within peaks for TF identificaton.
#'
#' @return (data.frame) Network, with "reg" and "target" in column name.
#'
#' @examples
#' table <- data.frame(Chr1=c("chr1", "chr1"), Position1=c(569265, 713603), Strand1=c("+", "+"),
#'                     Chr2=c("chr4", "chr1"), Position2=c(206628, 715110), Strand2=c("+", "-"),
#'                     row.names=c("A", "B"), stringsAsFactors=FALSE)
#' regulators=c("FOXF2", "MZF1")
#' #peaks and regulators to be analyzed
#' 
#' library(GenomicRanges)
#' library(GenomicFeatures)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(org.Hs.eg.db)
#' library(annotate)
#' promoter <- promoters(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
#' names(promoter) <- getSYMBOL(names(promoter), data="org.Hs.eg")
#' promoter <- promoter[!is.na(names(promoter))]
#' #get promoter regions
#' 
#' library(JASPAR2018)
#' library(TFBSTools)
#' library(motifmatchr)
#' pfm <- getMatrixSet(JASPAR2018, list(species="Homo sapiens"))
#' pfm <- pfm[unlist(lapply(pfm, function(x) x@name)) %in% regulators]
#' #get regulator position frequency matrix (PFM) list
#' 
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' conformation_network(table, promoter, pfm, "BSgenome.Hsapiens.UCSC.hg19")
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @importFrom motifmatchr matchMotifs
#' @importFrom motifmatchr motifMatches
#'
#' @export

conformation_network <- function(table, promoter, pfm, genome, range=500, p.cutoff=5e-05, w=7){
  table1 <- makeGRangesFromDataFrame(data.frame(Chr=table$Chr1, Start=table$Position1-range, End=table$Position1+range, Strand=table$Strand1,
                                                row.names=rownames(table), stringsAsFactors=FALSE))
  table2 <- makeGRangesFromDataFrame(data.frame(Chr=table$Chr2, Start=table$Position2-range, End=table$Position2+range, Strand=table$Strand2,
                                                row.names=rownames(table), stringsAsFactors=FALSE))
  #convert record table to GRanges, "1" and "2" indicate the"read1" and "read2" can both be in TF and target regions.
  target1 <- as.data.frame(findOverlaps(table1, promoter))
  target1$subjectHits <- names(promoter)[target1$subjectHits]
  target1$queryHits <- names(table1)[target1$queryHits]
  target2 <- as.data.frame(findOverlaps(table2, promoter))
  target2$subjectHits <- names(promoter)[target2$subjectHits]
  target2$queryHits <- names(table2)[target2$queryHits]
  #map between ranges and promoters by GRanges overlapping searching
  reg1 <- matchMotifs(pfm, table1[unique(target2$queryHits)], genome=genome, p.cutoff=p.cutoff, w=w)
  reg1 <- as.matrix(motifMatches(reg1))
  reg2 <- matchMotifs(pfm, table2[unique(target1$queryHits)], genome=genome, p.cutoff=p.cutoff, w=w)
  reg2 <- as.matrix(motifMatches(reg2))
  colnames(reg1) <- colnames(reg2) <- structure(unlist(lapply(pfm, function(x) x@name)), names=NULL)
  #map between ranges and regulators by motif searching
  regul1 <- lapply(unique(target1$queryHits), function(x, target1, reg2){
    expand.grid(target=target1$subjectHits[target1$queryHits == x], reg=colnames(reg2)[reg2[x, ]], stringsAsFactors=FALSE)}, target1=target1, reg2=reg2)
  regul2 <- lapply(unique(target2$queryHits), function(x, target2, reg1){
    expand.grid(target=target2$subjectHits[target2$queryHits == x], reg=colnames(reg1)[reg1[x, ]], stringsAsFactors=FALSE)}, target2=target2, reg1=reg1)
  regul <- rbind(do.call(rbind, regul1), do.call(rbind, regul2))
  regul <- regul[!duplicated(paste(regul$reg, regul$target, sep="-")), ]
  #regulator-target network
  return(regul)}

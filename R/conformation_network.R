#' Build Network from Conformation Peaks.
#'
#' Build network from conformation, e.g. HiChIP records.
#'
#' @param table (data.frame) Records, with "Chr1", "Position1", "Strand1", "Chr2", "Position2" and "Strand2" in column name, and record ID in row names.
#'
#' @param regulators (character) Regulators to be analyzed.
#'
#' @param species (character) Could be "hg19", "hg38", "mm9" or "mm10".
#'
#' @param range (numeric) Search region defined from "Position1" and "Position2".
#'
#' @param upstream @param downstream (numeric) Promoter region defined by distance from TSS.
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
#'                     row.names=c("A", "B"), stringsAsFactors=F)
#' regulators=c("FOXF2", "MZF1")
#' conformation_network(table, regulators, species = "hg19")
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import TxDb.Mmusculus.UCSC.mm9.knownGene
#' @import BSgenome.Mmusculus.UCSC.mm9
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import TxDb.Mmusculus.UCSC.mm10.knownGene
#' @import BSgenome.Mmusculus.UCSC.mm10
#' @import JASPAR2018
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicRanges promoters
#' @importFrom GenomicRanges findOverlaps
#' @importFrom annotate getSYMBOL
#' @importFrom TFBSTools getMatrixSet
#' @importFrom motifmatchr matchMotifs
#' @importFrom motifmatchr motifMatches
#'
#' @export

conformation_network <- function(table, regulators, species=c("hg19", "mm9", "hg38", "mm10"), range=500, upstream=2000, downstream=200, p.cutoff=5e-05, w=7){
  table1 <- makeGRangesFromDataFrame(data.frame(Chr=table$Chr1, Start=table$Position1-range, End=table$Position1+range, Strand=table$Strand1,
                                                row.names=rownames(table), stringsAsFactors=F))
  table2 <- makeGRangesFromDataFrame(data.frame(Chr=table$Chr2, Start=table$Position2-range, End=table$Position2+range, Strand=table$Strand2,
                                                row.names=rownames(table), stringsAsFactors=F))
  #convert record table to GRanges, "1" and "2" indicate the"read1" and "read2" can both be in TF and target regions.
  if (species == "hg19"){
    library(org.Hs.eg.db)
    promoter <- promoters(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), upstream=upstream, downstream=downstream)
    names(promoter) <- getSYMBOL(names(promoter), data="org.Hs.eg")
  }else if (species == "mm9"){
    library(org.Mm.eg.db)
    promoter <- promoters(genes(TxDb.Mmusculus.UCSC.mm9.knownGene), upstream=upstream, downstream=downstream)
    names(promoter) <- getSYMBOL(names(promoter), data="org.Mm.eg")
  }else if (species == "hg38"){
    library(org.Hs.eg.db)
    promoter <- promoters(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), upstream=upstream, downstream=downstream)
    names(promoter) <- getSYMBOL(names(promoter), data="org.Hs.eg.db")
  }else if (species == "mm10"){
    library(org.Mm.eg.db)
    promoter <- promoters(genes(TxDb.Mmusculus.UCSC.mm10.knownGene), upstream=upstream, downstream=downstream)
    names(promoter) <- getSYMBOL(names(promoter), data="org.Mm.eg.db")
  }else message("only hg38 and mm10 are supported")
  promoter <- promoter[!is.na(names(promoter))]
  #get promoter region of genes
  target1 <- as.data.frame(findOverlaps(table1, promoter))
  target1$subjectHits <- names(promoter)[target1$subjectHits]
  target1$queryHits <- names(table1)[target1$queryHits]
  target2 <- as.data.frame(findOverlaps(table2, promoter))
  target2$subjectHits <- names(promoter)[target2$subjectHits]
  target2$queryHits <- names(table2)[target2$queryHits]
  #map between ranges and promoters by GRanges overlapping searching
  if (species == "hg19"){
    PFMatrixList <- getMatrixSet(JASPAR2018, list(species="Homo sapiens"))
    PFMatrixList <- PFMatrixList[unlist(lapply(PFMatrixList, function(x) x@name)) %in% regulators]
    reg1 <- matchMotifs(PFMatrixList, table1[unique(target2$queryHits)], genome="BSgenome.Hsapiens.UCSC.hg19", p.cutoff=p.cutoff, w=w)
    reg2 <- matchMotifs(PFMatrixList, table2[unique(target1$queryHits)], genome="BSgenome.Hsapiens.UCSC.hg19", p.cutoff=p.cutoff, w=w)
  }else if (species == "mm9"){
    PFMatrixList <- getMatrixSet(JASPAR2018, list(species="Mus musculus"))
    PFMatrixList <- PFMatrixList[unlist(lapply(PFMatrixList, function(x) x@name)) %in% regulators]
    reg1 <- matchMotifs(PFMatrixList, table1[unique(target2$queryHits)], genome="BSgenome.Mmusculus.UCSC.mm9", p.cutoff=p.cutoff, w=w)
    reg2 <- matchMotifs(PFMatrixList, table2[unique(target1$queryHits)], genome="BSgenome.Mmusculus.UCSC.mm9", p.cutoff=p.cutoff, w=w)
  }else if (species == "hg38"){
    PFMatrixList <- getMatrixSet(JASPAR2018, list(species="Homo sapiens"))
    PFMatrixList <- PFMatrixList[unlist(lapply(PFMatrixList, function(x) x@name)) %in% regulators]
    reg1 <- matchMotifs(PFMatrixList, table1[unique(target2$queryHits)], genome="BSgenome.Hsapiens.UCSC.hg38", p.cutoff=p.cutoff, w=w)
    reg2 <- matchMotifs(PFMatrixList, table2[unique(target1$queryHits)], genome="BSgenome.Hsapiens.UCSC.hg38", p.cutoff=p.cutoff, w=w)
  }else if (species == "mm10"){
    PFMatrixList <- getMatrixSet(JASPAR2018, list(species="Mus musculus"))
    PFMatrixList <- PFMatrixList[unlist(lapply(PFMatrixList, function(x) x@name)) %in% regulators]
    reg1 <- matchMotifs(PFMatrixList, table1[unique(target2$queryHits)], genome="BSgenome.Mmusculus.UCSC.mm10", p.cutoff=p.cutoff, w=w)
    reg2 <- matchMotifs(PFMatrixList, table2[unique(target1$queryHits)], genome="BSgenome.Mmusculus.UCSC.mm10", p.cutoff=p.cutoff, w=w)}
  reg1 <- as.matrix(motifMatches(reg1))
  reg2 <- as.matrix(motifMatches(reg2))
  colnames(reg1) <- colnames(reg2) <- structure(unlist(lapply(PFMatrixList, function(x) x@name)), names=NULL)
  #map between ranges and regulators by motif searching
  regul1 <- lapply(unique(target1$queryHits), function(x, target1, reg2){
    expand.grid(target=target1$subjectHits[target1$queryHits == x], reg=colnames(reg2)[reg2[x, ]], stringsAsFactors=F)}, target1=target1, reg2=reg2)
  regul2 <- lapply(unique(target2$queryHits), function(x, target2, reg1){
    expand.grid(target=target2$subjectHits[target2$queryHits == x], reg=colnames(reg1)[reg1[x, ]], stringsAsFactors=F)}, target2=target2, reg1=reg1)
  regul <- rbind(do.call(rbind, regul1), do.call(rbind, regul2))
  regul <- regul[!duplicated(paste(regul$reg, regul$target, sep="-")), ]
  #regulator-target network
  return(regul)}

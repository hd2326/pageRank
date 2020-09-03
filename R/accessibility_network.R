#' Build Network from Accessibility Peaks.
#'
#' Build network from accessibility, e.g. ATAC-Seq peaks.
#'
#' @param table (data.frame) Peaks, with "Chr", "Start" and "End" in column name, and peak ID in row names.
#'
#' @param regulators (character) Regulators to be analyzed.
#'
#' @param species (character) Could be "hg19", "hg38", "mm9" or "mm10".
#'
#' @param upstream (numeric) Promoter region from TSS for target identificaton.
#' 
#' @param downstream (numeric) Promoter region from TSS for target identificaton.
#'
#' @param p.cutoff (numeric) P-value cutoff for motifs searching within peaks for TF identificaton.
#'
#' @param w (numeric) Window size for motifs searching within peaks for TF identificaton.
#'
#' @return (data.frame) Network, with "reg" and "target" in column name.
#'
#' @examples
#' table <- data.frame(Chr=c("chr1", "chr1"), Start=c(713689, 856337), End=c(714685, 862152),
#'                     row.names=c("A", "B"), stringsAsFactors=F)
#' regulators=c("FOXF2", "MZF1")
#' accessibility_network(table, regulators, species = "hg19")
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

accessibility_network <- function(table, regulators, species=c("hg19", "mm9", "hg38", "mm10"),
                                  upstream=2000, downstream=200, p.cutoff=5e-05, w=7){
  table <- makeGRangesFromDataFrame(table)
  #convert peak table to GRanges
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
    names(promoter) <- getSYMBOL(names(promoter), data="org.Hs.eg")
  }else if (species == "mm10"){
    library(org.Mm.eg.db)
    promoter <- promoters(genes(TxDb.Mmusculus.UCSC.mm10.knownGene), upstream=upstream, downstream=downstream)
    names(promoter) <- getSYMBOL(names(promoter), data="org.Mm.eg")
  }else message("only hg19, hg38, mm9 and mm10 are supported")
  promoter <- promoter[!is.na(names(promoter))]
  #get promoter region of genes
  target <- as.data.frame(findOverlaps(table, promoter))
  target$subjectHits <- names(promoter)[target$subjectHits]
  target$queryHits <- names(table)[target$queryHits]
  #map between peaks and promoters by GRanges overlapping searching
  if (species == "hg19"){
    PFMatrixList <- getMatrixSet(JASPAR2018, list(species="Homo sapiens"))
    PFMatrixList <- PFMatrixList[unlist(lapply(PFMatrixList, function(x) x@name)) %in% regulators]
    reg <- matchMotifs(PFMatrixList, table[unique(target$queryHits)], genome="BSgenome.Hsapiens.UCSC.hg19")
  }else if (species == "mm9"){
    PFMatrixList <- getMatrixSet(JASPAR2018, list(species="Mus musculus"))
    PFMatrixList <- PFMatrixList[unlist(lapply(PFMatrixList, function(x) x@name)) %in% regulators]
    reg <- matchMotifs(PFMatrixList, table[unique(target$queryHits)], genome="BSgenome.Mmusculus.UCSC.mm9")
  }else if (species == "hg38"){
    PFMatrixList <- getMatrixSet(JASPAR2018, list(species="Homo sapiens"))
    PFMatrixList <- PFMatrixList[unlist(lapply(PFMatrixList, function(x) x@name)) %in% regulators]
    reg <- matchMotifs(PFMatrixList, table[unique(target$queryHits)], genome="BSgenome.Hsapiens.UCSC.hg38")
  }else if (species == "mm10"){
    PFMatrixList <- getMatrixSet(JASPAR2018, list(species="Mus musculus"))
    PFMatrixList <- PFMatrixList[unlist(lapply(PFMatrixList, function(x) x@name)) %in% regulators]
    reg <- matchMotifs(PFMatrixList, table[unique(target$queryHits)], genome="BSgenome.Mmusculus.UCSC.mm10")}
  reg <- as.matrix(motifMatches(reg))
  colnames(reg) <- structure(unlist(lapply(PFMatrixList, function(x) x@name)), names=NULL)
  #map between peaks and regulators by motif searching
  regul <- do.call(rbind, lapply(unique(target$queryHits), function(x, target, reg){
    expand.grid(target=target$subjectHits[target$queryHits == x], reg=colnames(reg)[reg[x, ]], stringsAsFactors=F)}, target=target, reg=reg))
  #regulator-target network
  return(regul)}

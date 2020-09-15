#' Build Network from Accessibility Peaks.
#'
#' Build network from accessibility, e.g. ATAC-Seq peaks.
#'
#' @param table (data.frame) Peaks, with "Chr", "Start" and "End" in column
#' name, and peak ID in row names.
#'
#' @param promoter (GRanges) Promoter regions.
#' 
#' @param pfm (PFMatrixList) Positon Frequency Matrices (PFMs) of regulators.
#' 
#' @param genome (BSgenome or character) Genome build in which regulator motifs
#' will be searched.
#'
#' @param p.cutoff (numeric) P-value cutoff for motifs searching within peaks
#' for TF identificaton.
#'
#' @param w (numeric) Window size for motifs searching within peaks for TF
#' identificaton.
#'
#' @return (data.frame) Network, with "reg" and "target" in column name.
#'
#' @examples
#' 
#' table <- data.frame(Chr=c("chr1", "chr1"), Start=c(713689, 856337),
#'                     End=c(714685, 862152), row.names=c("A", "B"),
#'                     stringsAsFactors=FALSE)
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
#' pfm <- pfm[unlist(lapply(pfm, function(x) name(x))) %in% regulators]
#' #get regulator position frequency matrix (PFM) list
#' 
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' accessibility_network(table, promoter, pfm, "BSgenome.Hsapiens.UCSC.hg19")
#' #generate network
#'
#' @author DING, HONGXU (hd2326@columbia.edu)
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @importFrom motifmatchr matchMotifs
#' @importFrom motifmatchr motifMatches
#'
#' @export

accessibility_network <- function(table, promoter, pfm, genome,
                                  p.cutoff=5e-05, w=7){
    table <- makeGRangesFromDataFrame(table)
    #convert peak table to GRanges
    target <- as.data.frame(findOverlaps(table, promoter))
    target$subjectHits <- names(promoter)[target$subjectHits]
    target$queryHits <- names(table)[target$queryHits]
    #map between peaks and promoters by GRanges overlapping searching
    reg <- matchMotifs(pfm, table[unique(target$queryHits)], genome=genome,
                       p.cutoff=p.cutoff, w=w)
    reg <- as.matrix(motifMatches(reg))
    colnames(reg) <- structure(unlist(lapply(pfm, function(x){
        name(x)})), names=NULL)
    #map between peaks and pfms by motif searching
    regul <- lapply(unique(target$queryHits), function(x, target, reg){
        expand.grid(target=target$subjectHits[target$queryHits == x],
                    reg=colnames(reg)[reg[x, ]],
                    stringsAsFactors=FALSE)}, target=target, reg=reg)
    regul <- do.call(rbind, regul)
    #regulator-target network
    return(regul)}

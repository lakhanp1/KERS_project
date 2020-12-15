library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(TxDb.Anidulans.FGSCA4.AspGD.GFF)
library(here)


rm(list = ls())

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
analysisName <- "KERS_tRNA"
workDir <- here::here("kdmB_analysis", "09_KERS_tRNA")
outPrefix <- paste(workDir, "/", analysisName, sep = "")

file_plotSamples <- paste(workDir, "/", "samples.txt", sep = "")
orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

matrixType <- "normalizedMatrix_5kb"
up <- 5000
body <- 2000
down <- 1000
binSize <- 10
matrixDim = c(c(up, body, down)/binSize, binSize)

showExpressionHeatmap = FALSE

## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")

anLables <- list()

##################################################################################

sampleList <- suppressMessages(readr::read_tsv(file = file_plotSamples, col_names = T, comment = "#"))

## read the experiment sample details and select only those which are to be plotted
tempSInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = sampleList$sampleId,
                                    dataPath = TF_dataPath,
                                    profileMatrixSuffix = matrixType)

polII_ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "polII")]
tfIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag %in% c("HA", "MYC", "TAP") & tempSInfo$TF != "untagged")]
inputIds <- tempSInfo$sampleId[which(! tempSInfo$IP_tag %in% c("polII", "HIST") & tempSInfo$TF == "untagged")]
histIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST")]

## read the experiment sample details and select only those which are to be plotted
tfData <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = tfIds,
                                 dataPath = TF_dataPath,
                                 profileMatrixSuffix = matrixType)


inputData <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = inputIds,
                                    dataPath = TF_dataPath,
                                    profileMatrixSuffix = matrixType)

polIIData <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = polII_ids,
                                    dataPath = polII_dataPath,
                                    profileMatrixSuffix = matrixType)

histData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = histIds,
                                   dataPath = hist_dataPath,
                                   profileMatrixSuffix = matrixType)

exptData <- dplyr::left_join(
  sampleList, dplyr::bind_rows(tfData, inputData, histData, polIIData), by = "sampleId"
)

exptDataList <- purrr::transpose(exptData) %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))


tfCols <- sapply(
  c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
    "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
    "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
  simplify = F, USE.NAMES = T)

##################################################################################

## prepare tRNA information
txInfo <- suppressMessages(
  AnnotationDbi::select(
    x = txDb, keys = AnnotationDbi::keys(x = txDb, keytype = "TXID"),
    columns = c("GENEID", "TXNAME", "TXTYPE", "TXCHROM"), keytype = "TXID")) %>%
  dplyr::mutate(TXID = as.character(TXID)) %>%
  dplyr::rename(geneId = GENEID, txName = TXNAME, txType = TXTYPE)

trnaInfo <- dplyr::filter(txInfo, TXCHROM != "mito_A_nidulans_FGSC_A4") %>% 
  dplyr::filter(txType %in% c("tRNA"))

trnaGr <- GenomicFeatures::genes(x = txDb, filter = list(gene_id = trnaInfo$geneId))
rtracklayer::export(object = trnaGr, format = "bed",
                    con = paste(workDir, "/AN_tRNA.bed", sep = ""))


trnaGr$name <- stringr::str_replace_all(
  string = trnaGr$gene_id, pattern = "\\(|\\)", replacement = "_"
)
trnaGr$codon <- stringr::str_replace(
  string = trnaGr$gene_id, pattern = "^\\w+.(\\w+).\\d+", replacement = "\\1"
)
trnaGr$aa <- stringr::str_replace(
  string = trnaGr$gene_id, pattern = "^t(\\w)\\W.+", replacement = "\\1"
)


hist(width(trnaGr))

trnaRegions <- GenomicRanges::resize(x = trnaGr, width = 1, fix = "center")

## get strands of the immediate upstream and downstream genes
txGr <- GenomicFeatures::transcripts(
  x = txDb, filter = list(tx_name = setdiff(txInfo$txName, trnaInfo$txName)))

downTx <- GenomicRanges::precede(x = trnaGr, subject = txGr, select = "first", ignore.strand = TRUE)
downHits <- tibble::tibble(
  trna = trnaGr$name, from = 1:length(trnaGr),
  tx = txGr$tx_name[downTx], to = downTx,
  strand_trna = as.vector(strand(trnaGr)),
  strand_tx = as.vector(strand(txGr[downTx])),
  dist = distance(x = trnaGr, y = txGr[downTx], ignore.strand = TRUE),
  orientation = "down"
)

upTx <- GenomicRanges::follow(x = trnaGr, subject = txGr, select = "last", ignore.strand = TRUE)
upHits <- tibble::tibble(
  trna = trnaGr$name, from = 1:length(trnaGr),
  tx = txGr$tx_name[upTx], to = upTx,
  strand_trna = as.vector(strand(trnaGr)),
  strand_tx = as.vector(strand(txGr[upTx])),
  dist = distance(x = trnaGr, y = txGr[upTx], ignore.strand = TRUE),
  orientation = "up"
)

ovlp <- GenomicRanges::findOverlaps(query = trnaGr, subject = txGr)
ovlpHits <- tibble::tibble(
  trna = trnaGr$name[ovlp@from], from = ovlp@from,
  tx = txGr$tx_name[ovlp@to], to = ovlp@to,
  strand_trna = as.vector(strand(trnaGr[ovlp@from])),
  strand_tx = as.vector(strand(txGr[ovlp@to])),
  dist = distance(x = trnaGr[ovlp@from], y = txGr[ovlp@to], ignore.strand = TRUE),
  orientation = "overlap"
)

if(!all(upHits$trna == downHits$trna)){
  stop("upHits and downHits tRNAs are different")
}

trnaNeighbors <- tibble::tibble(
  trna = upHits$trna, strand_trna = upHits$strand_trna,
  gene_up = upHits$tx, strand_up = upHits$strand_tx, dist_up = upHits$dist * -1,
  gene_down = downHits$tx, strand_down = downHits$strand_tx, dist_down = downHits$dist
) %>% 
  dplyr::group_by_at(.vars = vars(starts_with("strand_"))) %>% 
  dplyr::mutate(group = group_indices()) %>% 
  dplyr::ungroup()

##################################################################################

## generate tRNA region profile matrices
i <- 1

exptData$matFile <- paste(workDir, "/profile_matrix/", exptData$sampleId,
                          ".normalizedMatrix.tab.gz", sep = "")

# ## generate profile matrix
# for (i in 1:nrow(exptData)) {
# 
#   mat <- bigwig_profile_matrix(
#     bwFile = exptData$bwFile[i],
#     regions = trnaRegions,
#     signalName = exptData$sampleId[i],
#     extend = c(500, 500),
#     targetName = "center",
#     storeLocal = TRUE,
#     localPath = exptData$matFile[i]
#   )
# 
# }


matList <- import_profiles(exptInfo = exptData, geneList = trnaNeighbors$trna,
                           source = "normalizedmatrix",
                           targetType = "point", targetName = "center",
                           up = 50, target = 0, down = 50)


## tf colors
tfMeanProfile <- NULL
if(length(exptData$sampleId) == 1){
  tfMeanProfile <- matList[[exptData$sampleId[1]]]
} else{
  tfMeanProfile <- getSignalsFromList(lt = matList[exptData$sampleId])
}

quantile(tfMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# tfMeanColor <- colorRamp2(quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T), c("white", "red"))
tfColorList <- sapply(
  X = exptData$sampleId,
  FUN = function(x){
    return(colorRamp2(breaks = quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T),
                      colors = unlist(strsplit(x = exptDataList[[x]]$color, split = ",")))
    )
  }
)


profilePlots <- multi_profile_plots(
  exptInfo = exptData, genesToPlot = trnaNeighbors$trna,
  profileColors = tfColorList,
  targetType = "point",
  targetName = "center",
  matBins = c(50, 0, 50, 10), matSource = "normalizedmatrix",
  column_title_gp = gpar(fontsize = 12)
)

upAn <- rowAnnotation(
  upGenes = anno_points(
    pmax(trnaNeighbors$dist_up, -2000),
    width = unit(2, "cm"),
    ylim = c(-2000, 0),
    # pcs = dplyr::recode(trnaNeighbors$strand_up, `-` = 4, `+` = 3),
    gp = gpar(col = dplyr::recode(trnaNeighbors$strand_up, `-` = "red", `+` = "blue"))
  )
)

downAn <- rowAnnotation(
  upGenes = anno_points(
    pmin(trnaNeighbors$dist_down, 2000),
    width = unit(2, "cm"),
    ylim = c(0, 2000),
    # pcs = 1,
    # pcs = dplyr::recode(trnaNeighbors$strand_down, `-` = 4, `+` = 3),
    gp = gpar(col = dplyr::recode(trnaNeighbors$strand_down, `-` = "red", `+` = "blue"))
  )
)

htList <- upAn + profilePlots$heatmapList + downAn

# pdf(file = paste(outPrefix, "_profiles.pdf", sep = ""), width = 18, height = 13)
png(file = paste(outPrefix, ".profiles.png", sep = ""), width = 5000, height = 2500, res = 250)

ht <- draw(
  htList,
  main_heatmap = exptData$profileName[1],
  column_title = "KERS binding profile centered at tRNA",
  split = trnaNeighbors$group,
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_sub_title_side = "left",
  heatmap_legend_side = "bottom",
  gap = unit(7, "mm"),
  padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()
















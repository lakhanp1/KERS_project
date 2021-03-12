suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(summarytools))


rm(list = ls())

source(file = "https://raw.githubusercontent.com/lakhanp1/omics_utils/master/04_GO_enrichment/s01_enrichment_functions.R")

##################################################################################
## main configuration
comparisonName <- "kdmB_20h_del_vs_WT"
outDir <- here::here("analysis", "06_KERS_del_vs_WT", comparisonName)
outPrefix <- paste(outDir, "/", comparisonName, sep = "")

if(!dir.exists(outDir)) dir.create(path = outDir, recursive = T)

file_polIIRatioConf <- here::here("data", "reference_data", "polII_ratio.config.tab")

mainPolIIPair <- "polII_20h.kdmB_del_vs_untagged"

polIIPairs <- c(mainPolIIPair)

tfSamples <- c("An_kdmB_20h_HA_1")

# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
matrixType <- "normalizedmatrix_5kb"
matrixDim = c(500, 200, 100, 10)

## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/annotation_resources/geneid2go.ANidulans.topGO.map"

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db

showExpressionHeatmap <- FALSE

##################################################################################

## genes to read and annotations
geneSet <- suppressMessages(
  readr::read_tsv(file = file_genes, col_names = c("chr", "start", "end", "geneId", "score", "strand"))
)

geneDesc <- suppressMessages(
  AnnotationDbi::select(x = orgDb, keys = geneSet$geneId,
                        columns = c("DESCRIPTION"), keytype = "GID")
)

smInfo <- suppressMessages(
  AnnotationDbi::select(x = orgDb, keys = keys(orgDb, keytype = "SM_CLUSTER"),
                        columns = c("GID", "SM_ID"), keytype = "SM_CLUSTER")) %>% 
  dplyr::group_by(GID) %>% 
  dplyr::mutate(SM_CLUSTER = paste(SM_CLUSTER, collapse = ";"),
                SM_ID = paste(SM_ID, collapse = ";")) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup()

geneInfo <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID")) %>% 
  dplyr::left_join(y = smInfo, by = c("geneId" = "GID"))

glimpse(geneInfo)

## polII ratio configuration
polIIRatioConf <- suppressMessages(readr::read_tsv(file = file_polIIRatioConf)) %>% 
  dplyr::filter(comparison %in% polIIPairs)

polIIDiffPairs <- purrr::transpose(polIIRatioConf)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))

##################################################################################

## read the experiment sample details and select only those which are to be plotted
tfData <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSamples,
  dataPath = TF_dataPath,
  profileMatrixSuffix = matrixType
)

polIIData <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = unique(purrr::map(polIIDiffPairs, `[`, c("group2", "group1")) %>% unlist() %>% unname()),
  dataPath = polII_dataPath,
  profileMatrixSuffix = "normalized_profile"
)

exptData <- dplyr::bind_rows(tfData, polIIData)
exptDataList <- purrr::transpose(exptData)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

polII_ids <- exptData$sampleId[which(exptData$IP_tag == "polII")]
tfIds <- exptData$sampleId[which(exptData$IP_tag %in% c("HA", "MYC", "TAP") & exptData$TF != "untagged")]


polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)


tfCols <- sapply(
  c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
    "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
    "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
  simplify = F, USE.NAMES = T)


expressionData <- get_polII_expressions(exptInfo = exptData,
                                        genesDf = geneInfo)

## add fold change columns
for (i in names(polIIDiffPairs)) {
  expressionData <- get_fold_change(
    df = expressionData,
    nmt = polIIDiffPairs[[i]]$group1,
    dmt = polIIDiffPairs[[i]]$group2,
    newCol = polIIDiffPairs[[i]]$comparison,
    isExpressedCols = polIICols$is_expressed)
}


expressionData <- get_TF_binding_data(exptInfo = tfData,
                                      genesDf = expressionData)

glimpse(expressionData)
readr::write_tsv(x = expressionData, file = paste(outPrefix, ".data.tab", sep = ""))

##################################################################################
## colors for profile matrix
tfMats <- import_profiles(
  exptInfo = tfData,
  geneList = geneInfo$geneId,
  source = matrixType,
  up = matrixDim[1], target = matrixDim[2], down = matrixDim[3]
)

polIIMats <- import_profiles(
  exptInfo = polIIData,
  geneList = geneInfo$geneId,
  source = "deeptools",
  up = 200, target = 200, down = 100
)

matList <- c(tfMats, polIIMats)
scalledMatList <- purrr::map(.x = matList, .f = chipmine::scale_matrix_columns, add_attr = FALSE)

## tf colors
tfMeanProfile <- NULL
if(length(c(tfIds)) == 1){
  tfMeanProfile <- matList[[tfIds]]
} else{
  tfMeanProfile <- getSignalsFromList(lt = matList[tfData$sampleId])
}

quantile(tfMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
tfMeanColorList <- sapply(
  X = tfData$sampleId,
  FUN = function(x){
    return(
      colorRamp2(
        breaks = quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T),
        colors = unlist(strsplit(x = exptDataList[[x]]$color, split = ","))
      )
    )
  }
)


## polII colors
polIIMeanProfile <- NULL
polIIColorList <- NULL
if(nrow(polIIData) == 1){
  polIIMeanProfile <- matList[[polIIData$sampleId]]
} else{
  polIIMeanProfile <- getSignalsFromList(lt = matList[polIIData$sampleId])
}
quantile(polIIMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
polIIMeanColor <- colorRamp2(quantile(polIIMeanProfile, c(0.01, 0.5, 0.995), na.rm = T), c("blue", "white", "red"))
polIIColorList <- sapply(X = polIIData$sampleId, FUN = function(x){return(polIIMeanColor)})


##################################################################################
















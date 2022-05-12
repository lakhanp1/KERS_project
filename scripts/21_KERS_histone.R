suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))


## This script plots the profile heatmaps for multiple samples. If the expression values are present, it
## also plots a simple heatmap using these expression values

rm(list = ls())

##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
comparisonName <- "R_H3KAc_48h"
outDir <- here::here("analysis", "10_KERS_histone", comparisonName)

file_commonGenes <- file.path(outDir, "KERS_complex_48h.common_genes.tab")
file_plotSamples <- file.path(outDir, "samples.txt")

outPrefix <- file.path(outDir, comparisonName)


# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
regionMatType <- "deeptools"
regionMatSuffix <- "normalized_profile"
regionMatDim = c(200, 200, 100, 10)

tssMatType <- "normalizedmatrix"
tssMatSuffix <- "normalizedmatrix_3kbATG3kb"
tssMatDim <- c(300, 1, 300, 10)

showExpressionHeatmap <- FALSE

## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF


## colors
colList <- list()


##################################################################################
sampleList <- suppressMessages(readr::read_tsv(file = file_plotSamples, comment = "#"))

tempSInfo <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = sampleList$sampleId,
  dataPath = TF_dataPath, profileMatrixSuffix = regionMatSuffix
)

tfIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag %in% c("HA", "MYC", "TAP") & tempSInfo$TF != "untagged")]
inputIds <- tempSInfo$sampleId[which(! tempSInfo$IP_tag %in% c("polII", "HIST") & tempSInfo$TF == "untagged")]
polII_ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "polII")]
histIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST" & tempSInfo$TF != "H3")]
histH3Ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST" & tempSInfo$TF == "H3")]

## read the experiment sample details and select only those which are to be plotted
tfData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = tfIds,
  dataPath = TF_dataPath, profileMatrixSuffix = tssMatSuffix
)

inputData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = inputIds,
  dataPath = TF_dataPath, profileMatrixSuffix = tssMatSuffix
)

polIIData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = polII_ids,
  dataPath = polII_dataPath, profileMatrixSuffix = regionMatSuffix
)

histData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = histIds,
  dataPath = hist_dataPath, profileMatrixSuffix = tssMatSuffix
)

h3Data <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = histH3Ids,
  dataPath = hist_dataPath, profileMatrixSuffix = tssMatSuffix
)


regionProfileData <- dplyr::bind_rows(polIIData)
tssProfileData <- dplyr::bind_rows(tfData, inputData, histData, h3Data)

exptData <- dplyr::bind_rows(tfData, inputData, histData, h3Data, polIIData)

exptDataList <- purrr::transpose(exptData)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)


tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
  simplify = F, USE.NAMES = T)

##################################################################################

## genes to read
geneSet <- data.table::fread(
  file = file_genes, header = F,
  col.names = c("chr", "start", "end", "geneId", "score", "strand")
)

kersGenes <- suppressMessages(readr::read_tsv(file = file_commonGenes)) %>% 
  dplyr::mutate(binding = "KERS")

# kmClust <- dplyr::left_join(
#   x = suppressMessages(readr::read_tsv(file = tfData$clusterFile[1])),
#   y = geneSet, by = c("geneId" = "geneId")
# )

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, keytype = "GID",
                                  columns = c("GENE_NAME", "DESCRIPTION"))

geneInfo <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID")) %>% 
  dplyr::left_join(y = kersGenes, by = "geneId") %>% 
  tidyr::replace_na(replace = list(binding = "not_KERS"))

glimpse(geneInfo)

expressionData <- get_TF_binding_data(
  exptInfo = tfData, genesDf = geneInfo
)

expressionData <- get_polII_expressions(
  exptInfo = polIIData, genesDf = expressionData
) %>% 
  dplyr::mutate_at(
    .vars = vars(unname(tfCols$hasPeak)),
    .funs = ~replace_na(data = ., replace = FALSE)
  )

# glimpse(expressionData)


##################################################################################
## color list
tssMatList <- NULL
regionMatList <- NULL

regionMatList <- import_profiles(
  exptInfo = regionProfileData,
  geneList = expressionData$geneId,
  source = regionMatType,
  up = regionMatDim[1], target = regionMatDim[2], down = regionMatDim[3],
  targetType = "region", targetName = "gene"
)

tssMatList <- import_profiles(
  exptInfo = tssProfileData,
  geneList = expressionData$geneId,
  up = tssMatDim[1], target = tssMatDim[2], down = tssMatDim[3],
  targetType = "point", targetName = "ATG"
)

matList <- c(regionMatList, tssMatList)

##################################################################################

## tf colors
tfMeanProfile <- NULL
tfColorList <- NULL
if(length(c(tfIds)) == 1){
  tfMeanProfile <- matList[[tfIds]]
} else{
  tfMeanProfile <- getSignalsFromList(lt = matList[tfIds])
}

quantile(tfMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# tfMeanColor <- colorRamp2(quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T), c("white", "red"))
tfColorList <- sapply(
  X = c(tfIds, inputIds),
  FUN = function(x){
    return(
      colorRamp2(breaks = quantile(tfMeanProfile, c(0.50, 0.99), na.rm = T),
                 colors = unlist(strsplit(x = exptDataList[[x]]$color, split = ","))))
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
polIIMeanColor <- colorRamp2(quantile(polIIMeanProfile, c(0.01, 0.5, 0.99), na.rm = T), c("blue", "white", "red"))
polIIColorList <- sapply(X = polIIData$sampleId, FUN = function(x){return(polIIMeanColor)})

## histone colors
histMeanProfile <- NULL
histColorList <- NULL
if(nrow(histData) == 1){
  histMeanProfile <- matList[[histIds]]
} else{
  histMeanProfile <- getSignalsFromList(lt = matList[histIds])
}
quantile(histMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# histMeanColor <- colorRamp2(quantile(histMeanProfile, c(0.30, 0.995), na.rm = T), c("black", "yellow"))
histColorList <- sapply(
  X = c(histIds, histH3Ids),
  FUN = function(x){
    return(colorRamp2(breaks = quantile(histMeanProfile, c(0.60, 0.90), na.rm = T),
                      colors =  unlist(strsplit(x = exptDataList[[x]]$color, split = ","))))
  }
)

colorList <- unlist(list(tfColorList, polIIColorList, histColorList))
ylimList <- list()
# ylimList <- sapply(c(histIds, histH3Ids), function(x){return(c(5, 13))}, simplify = FALSE)
# ylimList <- append(x = ylimList,
#                    values = sapply(c(tfIds, inputIds), function(x){return(c(0, 25))}, simplify = FALSE))


##################################################################################
# plot genes which has TF peak in any of the TF samples
commonPeaks <- dplyr::filter(
  .data = expressionData, binding == "KERS"
)

nGenes <- nrow(commonPeaks)

## genes which do not have peak and no polII signal
noPeakPolIIDf <-dplyr::filter_at(
  .tbl = expressionData,
  .vars = unname(c(tfCols$hasPeak, polIICols$is_expressed)),
  .vars_predicate = all_vars(. == FALSE)
) %>% 
  dplyr::filter(binding == "not_KERS") %>% 
  dplyr::slice_min(
    order_by = !!sym(unname(polIICols$exp)), n = nGenes, with_ties = FALSE
  )

plotData <- dplyr::bind_rows(commonPeaks, noPeakPolIIDf) %>% 
  dplyr::mutate(
    binding = forcats::fct_relevel(.f = binding, "KERS", "not_KERS")
  )

regionProfiles_peak <- multi_profile_plots(
  exptInfo = regionProfileData,
  genesToPlot = plotData$geneId,
  matSource = regionMatType,
  matBins = regionMatDim,
  drawClusterAn = FALSE,
  profileColors = colorList,
  expressionColor = NULL,
  plotExpression = showExpressionHeatmap,
  column_title_gp = gpar(fontsize = 12),
  ylimFraction = ylimList
)

tssProfiles_peak <- multi_profile_plots(
  exptInfo = tssProfileData,
  genesToPlot = plotData$geneId,
  targetType = "point",
  targetName = "ATG",
  matSource = tssMatType,
  matBins = tssMatDim,
  clusters = dplyr::select(plotData, geneId, cluster = binding),
  drawClusterAn = TRUE,
  profileColors = colorList,
  column_title_gp = gpar(fontsize = 12),
  ylimFraction = ylimList
)



## gene length annotation
anGl_peaks <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = plotData$geneId,
  # pointSize = unit(4, "mm"),
  axis_param = list(at = c(2000, 4000), labels = c("2kb", "> 4kb"))
)

peaks_htlist <- NULL
peaks_htlist <- peaks_htlist + anGl_peaks$an
peaks_htlist <- peaks_htlist + tssProfiles_peak$heatmapList
peaks_htlist <- peaks_htlist + regionProfiles_peak$heatmapList

## make sure that the order of genes in the heatmap list and in the dataframe is same
if(all(rownames(peaks_htlist@ht_list[[ tfData$profileName[1] ]]@matrix) == plotData$geneId)){
  
  rowOrd_peaks <- order(plotData[[ tfCols$peakDist[[1]] ]], decreasing = TRUE)
  
}

pdfWd <- 2 + 
  (length(peaks_htlist) * 2) +
  (length(polII_ids) * 0.25 * showExpressionHeatmap)

# wd <- 500 + (nrow(exptData) * 700) + (length(polII_ids) * 500 * showExpressionHeatmap)
title_peak= paste(comparisonName, ": KERS target genes (n = ", nGenes, ")", sep = "")


pdf(file = paste(outPrefix, ".common_peaks.profiles.pdf", sep = ""), width = pdfWd, height = 13)

draw(
  peaks_htlist,
  main_heatmap = exptData$profileName[1],
  # annotation_legend_list = list(profile1$legend),
  column_title = title_peak,
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_sub_title_side = "left",
  heatmap_legend_side = "bottom",
  gap = unit(7, "mm"),
  # row_order = rowOrd_peaks,
  # split = rep(1, nrow(plotData)),
  padding = unit(rep(0.5, times = 4), "cm")
)


dev.off()




##################################################################################







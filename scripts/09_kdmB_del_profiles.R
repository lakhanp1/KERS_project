library(chipmine)
library(here)
library(summarytools)

## This script plots the profile heatmaps for kdmB_del and sntB_del samples
## three plots are generate: all genes, macs2 target genes: grouped and ungrouped

rm(list = ls())

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
comparisonName <- "kdmB_del_sntB_del_48h_combined"
outPrefix <- here::here("kdmB_analysis", comparisonName, comparisonName)


file_plotSamples <- here::here("kdmB_analysis", comparisonName, "samples.txt")

# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
matrixType <- "normalizedmatrix_5kb"
matrixDim = c(500, 200, 100, 10)

showExpressionHeatmap <- FALSE


## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"


TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")


## colors
colList <- list()

outPrefix_all <- paste0(outPrefix, "_allGenes", collapse = "")
outPrefix_expressed <- paste0(outPrefix, "_expressedGenes", collapse = "")
outPrefix_sm <- paste0(outPrefix, "_SM_genes", collapse = "")
outPrefix_peaks <- paste0(outPrefix, "_peaksGenes", collapse = "")
outPrefix_pkExp <- paste0(outPrefix, "_pkExpGenes", collapse = "")



##################################################################################
sampleList <- readr::read_tsv(file = file_plotSamples, col_names = T, comment = "#")

tempSInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = sampleList$sampleId,
                                    dataPath = TF_dataPath,
                                    matrixSource = matrixType)

tfIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag %in% c("HA", "MYC", "TAP") & tempSInfo$TF != "untagged")]
inputIds <- tempSInfo$sampleId[which(! tempSInfo$IP_tag %in% c("polII", "HIST") & tempSInfo$TF == "untagged")]
polII_ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "polII")]
histIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST")]


## read the experiment sample details and select only those which are to be plotted
tfData <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = tfIds,
                                 dataPath = TF_dataPath,
                                 matrixSource = matrixType)


inputData <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = inputIds,
                                    dataPath = TF_dataPath,
                                    matrixSource = matrixType)

polIIData <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = polII_ids,
                                    dataPath = polII_dataPath,
                                    matrixSource = matrixType)

histData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = histIds,
                                   dataPath = hist_dataPath,
                                   matrixSource = matrixType)

# histTssData <- get_sample_information(exptInfoFile = file_exptInfo,
#                                       samples = histIds,
#                                       dataPath = hist_dataPath,
#                                       matrixSource = tssMatType, profileType = "TSS_profile")
# 
# histTesData <- get_sample_information(exptInfoFile = file_exptInfo,
#                                       samples = histIds,
#                                       dataPath = hist_dataPath,
#                                       matrixSource = tesMatType, profileType = "TES_profile")


exptData <- dplyr::bind_rows(tfData, inputData, histData, polIIData)

exptDataList <- purrr::transpose(exptData)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)


tfCols <- sapply(c("hasPeak", "pval", "peakType", "tesPeakType", "peakDist", "summitDist", "upstreamExpr", "peakExpr", "relativeDist"),
                 FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
                 simplify = F, USE.NAMES = T)

##################################################################################

## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::mutate(length = end - start)

kmClust <- dplyr::left_join(x = readr::read_tsv(file = tfData$clusterFile[1], col_names = T),
                            y = geneSet,
                            by = c("gene" = "gene"))

## gene information annotations: cluster and TF and polII expression values
geneInfo <- add_gene_info(file = file_geneInfo, clusterDf = kmClust)

head(geneInfo)


expressionData <- get_TF_binding_data(exptInfo = tfData,
                                      genesDf = geneInfo)

# expressionData <- get_polII_expressions(exptInfo = polIIData,
#                                         genesDf = expressionData)

view(dfSummary(expressionData))


anLables <- list()
# anLables[[tssPeakTypeCol]] = gsub("peakType", "TSS peak type\n", tssPeakTypeCol) %>% gsub("\\(|\\)", "", .)
# anLables[[tesPeakTypeCol]] = gsub("tesPeakType", "TES peak type\n", tesPeakTypeCol) %>% gsub("\\(|\\)", "", .)
# anLables[[isExpCol]] = txt = gsub("is_expressed", "is expressed\n", isExpCol) %>% gsub("\\(|\\)", "", .)
anLables[["is_SM_gene"]] = "SM gene"
anLables[["is_TF"]] = "Transcription Factor"

##################################################################################
## color list
matList <- import_profiles(exptInfo = dplyr::bind_rows(tfData, inputData, histData, polIIData),
                               geneList = geneInfo$gene,
                               source = matrixType,
                               up = matrixDim[1], target = matrixDim[2], down = matrixDim[3])


## tf colors
tfMeanProfile <- NULL
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
    return(colorRamp2(breaks = quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T),
                      colors = c("white", exptDataList[[x]]$color)))
  }
)

## polII colors
polIIMeanProfile <- NULL
polIIColorList <- NULL
# if(nrow(polIIData) == 1){
#   polIIMeanProfile <- matList[[polIIData$sampleId]]
# } else{
#   polIIMeanProfile <- getSignalsFromList(lt = matList[polIIData$sampleId])
# }
# quantile(polIIMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# polIIMeanColor <- colorRamp2(quantile(polIIMeanProfile, c(0.01, 0.5, 0.995), na.rm = T), c("blue", "white", "red"))
# polIIColorList <- sapply(X = polIIData$sampleId, FUN = function(x){return(polIIMeanColor)})

## histone colors
histMeanProfile <- NULL
histColorList <- NULL
if(nrow(histData) == 1){
  histMeanProfile <- matList[[histData$sampleId]]
} else{
  histMeanProfile <- getSignalsFromList(lt = matList[histData$sampleId])
}
quantile(histMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# histMeanColor <- colorRamp2(quantile(histMeanProfile, c(0.30, 0.995), na.rm = T), c("black", "yellow"))
histColorList <- sapply(
  X = histData$sampleId,
  FUN = function(x){
    return(colorRamp2(breaks = quantile(histMeanProfile, c(0.30, 0.995), na.rm = T),
                      colors =  c("black", exptDataList[[x]]$color)))
  }
)

colorList <- unlist(list(tfColorList, polIIColorList, histColorList))
ylimList <- list()
# ylimList <- sapply(c(polII_ids, histIds), function(x){return(0.996)}, simplify = FALSE)
ylimList <- append(x = ylimList,
                   values = sapply(c(tfIds, inputIds), function(x){return(c(0, 25))}, simplify = FALSE))



##################################################################################

# expressionData$group <- dplyr::group_by_at(expressionData, .vars = unname(polIICols$is_expressed)) %>% 
#   dplyr::group_indices()
# 
# newClusters <- dplyr::select(expressionData, gene, group) %>% 
#   dplyr::rename(cluster = group)

multiProfiles_all <- multi_profile_plots(exptInfo = exptData,
                                         expressionData = expressionData,
                                         genesToPlot = geneInfo$gene,
                                         matSource = matrixType,
                                         matBins = matrixDim,
                                         clusters = expressionData,
                                         clusterColor = NULL,
                                         profileColors = colorList,
                                         expressionColor = NULL,
                                         plotExpression = showExpressionHeatmap,
                                         column_title_gp = gpar(fontsize = 12),
                                         ylimFraction = ylimList)


## gene length annotation
anGl <- gene_length_heatmap_annotation(bedFile = file_genes, genes = expressionData$gene)


htlist_allGenes <- anGl$an +
  multiProfiles_all$heatmapList
# tssTesPlotList

# ## row order as decreasing polII signal
# if( all(rownames(htlist_allGenes@ht_list[[ polIIData$profileName[1] ]]@matrix) == expressionData$gene) ){
#   
#   rowOrd <- order(expressionData[[polIICols$exp[1]]], decreasing = TRUE)
# }

pdfWd <- 2 + 
  (length(multiProfiles_all$heatmapList@ht_list) * 2) +
  (length(polII_ids) * 0.25 * showExpressionHeatmap)

title_all <- paste(comparisonName, ": all genes", collapse = "")

# draw Heatmap and add the annotation name decoration
# png(filename = paste(outPrefix_all, "_profiles.png", sep = ""), width=wd, height=3500, res = 250)
pdf(file = paste(outPrefix_all, "_profiles.pdf", sep = ""), width = pdfWd, height = 13)

draw(htlist_allGenes,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_all,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "bottom",
     gap = unit(7, "mm"),
     # row_order = rowOrd,
     padding = unit(rep(0.5, times = 4), "cm")
)


row_annotation_axis(an = "gene_length",
                    at = c(0, 2000, 4000),
                    labels = c("0kb", "2kb", ">4kb"),
                    slice = length(unique(geneInfo$cluster))
                    # slice = length(unique(newClusters$cluster))
)

dev.off()

##################################################################################
# plot genes which has TF peak in any of the TF samples
hasPeakDf <- filter_at(.tbl = expressionData,
                       .vars = unname(tfCols$hasPeak[ sampleList$sampleId[sampleList$usePeaks] ]),
                       .vars_predicate = any_vars(. == "TRUE"))

# plot genes which has TF peak in a specific TF sample
# hasPeakDf <- filter_at(.tbl = expressionData, .vars = unname(tfCols$hasPeak[1]), .vars_predicate = any_vars(. == "TRUE"))


multiProfiles_peak <- multi_profile_plots(exptInfo = exptData,
                                          expressionData = hasPeakDf,
                                          genesToPlot = hasPeakDf$gene,
                                          matSource = matrixType,
                                          matBins = matrixDim,
                                          clusters = hasPeakDf,
                                          clusterColor = multiProfiles_all$profileHeatmaps[[1]]$clusterColor,
                                          profileColors = colorList,
                                          expressionColor = multiProfiles_all$expressionColor,
                                          column_title_gp = gpar(fontsize = 12),
                                          plotExpression = showExpressionHeatmap,
                                          ylimFraction = ylimList)

## gene length annotation
anGl_peaks <- gene_length_heatmap_annotation(bedFile = file_genes, genes = hasPeakDf$gene)


peaks_htlist <- anGl_peaks$an +
  multiProfiles_peak$heatmapList


## make sure that the order of genes in the heatmap list and in the dataframe is same
if(all(rownames(peaks_htlist@ht_list[[ tfData$profileName[1] ]]@matrix) == hasPeakDf$gene)){
  
  rowOrd_peaks <- order(hasPeakDf[[ tfCols$peakDist[[1]] ]], decreasing = TRUE)
  
}


# wd <- 500 + (nrow(exptData) * 700) + (length(polII_ids) * 500 * showExpressionHeatmap)
title_peak= paste(comparisonName, ": genes with binding signal (macs2 peak targets)", collapse = "")

# draw Heatmap and add the annotation name decoration
pdf(file = paste(outPrefix_peaks, "_profiles.pdf", sep = ""), width = pdfWd, height = 13)

draw(peaks_htlist,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_peak,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "bottom",
     gap = unit(7, "mm"),
     # row_order = rowOrd_peaks,
     padding = unit(rep(0.5, times = 4), "cm")
)

row_annotation_axis(an = "gene_length",
                    at = c(0, 2000, 4000),
                    labels = c("0kb", "2kb", ">4kb"),
                    slice = length(unique(hasPeakDf$cluster)))
dev.off()


##################################################################################
## ungrouped profile plot
peaks_htlist2 <- anGl_peaks$an

peaks_htlist2 <- anGl_peaks$an
i <- 1
for(i in 1:nrow(exptData)){
  peaks_htlist2 <- peaks_htlist2 + multiProfiles_peak$profileHeatmaps[[exptData$sampleId[i]]]$heatmap
}


# draw Heatmap and add the annotation name decoration
pdf(file = paste(outPrefix_peaks, "_profiles_ungrouped.pdf", sep = ""), width = pdfWd, height = 13)

draw(peaks_htlist2,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_peak,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "bottom",
     gap = unit(7, "mm"),
     split = rep(1, nrow(hasPeakDf)),
     # row_order = rowOrd_peaks,
     padding = unit(rep(0.5, times = 4), "cm")
)

row_annotation_axis(an = "gene_length",
                    at = c(0, 2000, 4000),
                    labels = c("0kb", "2kb", ">4kb"),
                    slice = 1)
dev.off()


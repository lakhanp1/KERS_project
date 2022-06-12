suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(gghalves))
suppressPackageStartupMessages(library(ggpubr))


## compare Histone signal in different KERS binding group combinations

rm(list = ls())

##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
timePoint <- "20h"

analysisName <- sprintf(fmt = "KERS_groups_vs_H3KAc_%s", timePoint)
outDir <- here::here("analysis", "14_KERS_groups_vs_data", analysisName)

file_plotSamples <- file.path(outDir, "samples.txt")
file_kersBinding <- sprintf(
  fmt = here::here("analysis", "03_KERS_complex", "KERS_complex_%s",
                   "KERS_complex_%s.peaks_data.tab"),
  timePoint, timePoint
)

outPrefix <- file.path(outDir, analysisName)


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

kersGenes <- suppressMessages(readr::read_tsv(file = file_kersBinding))

# kmClust <- dplyr::left_join(
#   x = suppressMessages(readr::read_tsv(file = tfData$clusterFile[1])),
#   y = geneSet, by = c("geneId" = "geneId")
# )

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, keytype = "GID",
                                  columns = c("GENE_NAME", "DESCRIPTION"))

geneInfo <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID")) %>% 
  dplyr::left_join(y = kersGenes, by = "geneId") %>% 
  tidyr::replace_na(replace = list(bindingGroup = "0000"))

glimpse(geneInfo)

expressionData <- geneInfo

genesGr <- rtracklayer::import.bed(con = file_genes)
tssUp <- GenomicFeatures::promoters(x = genesGr, upstream = 500, downstream = 0)
tssDown <- GenomicFeatures::promoters(x = genesGr, upstream = 0, downstream = 500)


##################################################################################
## color list
tssMatList <- NULL
regionMatList <- NULL

# regionMatList <- import_profiles(
#   exptInfo = regionProfileData,
#   geneList = expressionData$geneId,
#   source = regionMatType,
#   up = regionMatDim[1], target = regionMatDim[2], down = regionMatDim[3],
#   targetType = "region", targetName = "gene"
# )

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
} else if(length(c(tfIds)) > 1){
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
if(length(polII_ids) == 1){
  polIIMeanProfile <- matList[[polIIData$sampleId]]
} else if(length(polII_ids) > 1){
  polIIMeanProfile <- getSignalsFromList(lt = matList[polIIData$sampleId])
}
quantile(polIIMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
polIIMeanColor <- colorRamp2(quantile(polIIMeanProfile, c(0.01, 0.5, 0.99), na.rm = T), c("blue", "white", "red"))
polIIColorList <- sapply(X = polIIData$sampleId, FUN = function(x){return(polIIMeanColor)})

## histone colors
histMeanProfile <- NULL
histColorList <- NULL
if(length(histIds) == 1){
  histMeanProfile <- matList[[histIds]]
} else if(length(histIds) > 1){
  histMeanProfile <- getSignalsFromList(lt = matList[histIds])
}
quantile(histMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# histMeanColor <- colorRamp2(quantile(histMeanProfile, c(0.30, 0.995), na.rm = T), c("black", "yellow"))
histColorList <- sapply(
  X = c(histIds, histH3Ids),
  FUN = function(x){
    return(
      colorRamp2(
        breaks = quantile(histMeanProfile, c(0.60, 0.90), na.rm = T),
        colors =  unlist(strsplit(x = exptDataList[[x]]$color, split = ","))
      )
    )
  }
)

colorList <- unlist(list(tfColorList, polIIColorList, histColorList))
ylimList <- list()
# ylimList <- sapply(c(histIds, histH3Ids), function(x){return(c(5, 13))}, simplify = FALSE)
# ylimList <- append(x = ylimList,
#                    values = sapply(c(tfIds, inputIds), function(x){return(c(0, 25))}, simplify = FALSE))


##################################################################################
## histone coverage comparison plot
coverageDf <- chipmine::region_coverage_matrix(
  regions = tssDown, exptInfo = exptData
) %>% 
  dplyr::select(-score)


expressionData <- dplyr::left_join(
  x = expressionData, y = coverageDf, by = c("geneId" = "name")
)


plotData <- dplyr::filter(
  expressionData, bindingGroup %in% c("KERS", "00RS")
) %>% 
  dplyr::mutate(
    bindingGroup = forcats::fct_relevel(.f = bindingGroup, "KERS")
  )

ggplotDf <- dplyr::select(
  plotData, geneId, bindingGroup, exptData$sampleId
) %>% 
  tidyr::pivot_longer(
    cols = !c(geneId, bindingGroup),
    names_to = "sampleId",
    values_to = "coverage"
  ) %>% 
  dplyr::mutate(
    sampleId = forcats::fct_relevel(.f = sampleId, sampleList$sampleId),
    bindingGroup = forcats::fct_relevel(.f = bindingGroup, "KERS")
  )


ggpubr::compare_means(
  formula = coverage ~ bindingGroup, data = ggplotDf, group.by = "sampleId",
  method = "wilcox.test", paired = FALSE
)

ks.test(
  x = dplyr::filter(ggplotDf, sampleId == "An_H3K9Ac_20h_HIST_1", bindingGroup == "KERS") %>% 
    dplyr::pull(coverage),
  y = dplyr::filter(ggplotDf, sampleId == "An_H3K9Ac_20h_HIST_1", bindingGroup == "00RS") %>% 
    dplyr::pull(coverage)
)


pt_violin <- ggplot2::ggplot(
  data = ggplotDf,
  mapping = aes(x = sampleId, y = coverage,  group = bindingGroup)
) +
  geom_quasirandom(mapping = aes(color = bindingGroup), dodge.width = 1) +
  geom_boxplot(
    position = position_dodge(1), width=0.3, 
    color = "black", fill = NA, outlier.shape = NA
  ) +
  stat_compare_means(paired = FALSE, label = "p.format", size = 6) +
  guides(color = guide_legend(override.aes = list(size = 5) ) ) +
  facet_grid(cols = vars(sampleId), scales = "free") +
  labs(
    title = "Histone acetylation vs KERS binding"
  ) +
  scale_color_manual(
    name = "binding",
    values = c("KERS" = "#E69F00", "00RS" = "#56B4E9")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom"
  )


ggsave(
  filename = paste(outPrefix, ".violine.pdf", sep = ""),
  plot = pt_violin, width = 10, height = 8
)

##################################################################################
## profile heatmap

tssProfiles_peak <- multi_profile_plots(
  exptInfo = tssProfileData,
  genesToPlot = plotData$geneId,
  targetType = "point",
  targetName = "ATG",
  matSource = tssMatType,
  matBins = tssMatDim,
  clusters = dplyr::select(plotData, geneId, cluster = bindingGroup),
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
# peaks_htlist <- peaks_htlist + regionProfiles_peak$heatmapList



pdfWd <- 2 + 
  (length(peaks_htlist) * 2) +
  (length(polII_ids) * 0.25 * showExpressionHeatmap)

# wd <- 500 + (nrow(exptData) * 700) + (length(polII_ids) * 500 * showExpressionHeatmap)
title_peak= paste(analysisName, ": Histone modification vs KERS binding groups", sep = "")


pdf(file = paste(outPrefix, ".profiles.pdf", sep = ""), width = pdfWd, height = 13)

draw(
  peaks_htlist,
  main_heatmap = exptData$profileName[1],
  # annotation_legend_list = list(profile1$legend),
  column_title = title_peak,
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_sub_title_side = "left",
  heatmap_legend_side = "bottom",
  gap = unit(7, "mm"),
  padding = unit(rep(0.5, times = 4), "cm")
)


dev.off()




##################################################################################

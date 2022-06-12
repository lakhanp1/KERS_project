suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(summarytools))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(gghalves))
suppressPackageStartupMessages(library(ggpubr))

## KERS binding groups vs polII ChIPseq signal

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")
##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
analysisName <- "KERS_group_vs_polII_20h"
outDir <- here::here("analysis", "03_KERS_complex", "KERS_complex_20h")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_plotSamples <- paste(outDir, "/", "samples.txt", sep = "")
file_kersBinding <- paste(outDir, "/", "KERS_complex_20h.peaks_data.tab", sep = "")


regionMatType <- "deeptools"
regionMatSuffix <- "normalized_profile"
regionMatDim = c(200, 200, 100, 10)

tssMatType <- "normalizedmatrix"
tssMatSuffix <- "normalizedmatrix_3kbATG3kb"
tssMatDim <- c(300, 1, 300, 10)

showExpressionHeatmap = TRUE

## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db

anLables <- list()

##################################################################################

sampleList <- suppressMessages(readr::read_tsv(file = file_plotSamples, comment = "#"))

## genes to read
geneSet <- suppressMessages(
  readr::read_tsv(file = file_genes, col_names = c("chr", "start", "end", "geneId", "score", "strand"))
) %>% 
  dplyr::select(geneId)

geneDesc <- suppressMessages(
  AnnotationDbi::select(x = orgDb, keys = geneSet$geneId,
                        columns = c("DESCRIPTION"), keytype = "GID")
)

geneInfo <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID"))

glimpse(geneInfo)


##################################################################################

## read the experiment sample details and select only those which are to be plotted
tempSInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = sampleList$sampleId,
  dataPath = TF_dataPath
)


polII_ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "polII")]
tfIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag %in% c("HA", "MYC", "TAP") & tempSInfo$TF != "untagged")]
inputIds <- tempSInfo$sampleId[which(! tempSInfo$IP_tag %in% c("polII", "HIST") & tempSInfo$TF == "untagged")]
histIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST")]
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

kersGenes <- suppressMessages(readr::read_tsv(file = file_kersBinding)) %>% 
  dplyr::select(geneId, bindingGroup)

geneInfo <- dplyr::left_join(x = geneInfo, y = kersGenes, by = "geneId") %>% 
  tidyr::replace_na(replace = list(bindingGroup = "!!!!"))

expressionData <- get_polII_expressions(exptInfo = exptData,
                                        genesDf = geneInfo)

plotData <- dplyr::filter(
  expressionData, bindingGroup %in% c("KERS", "!!RS")
) %>% 
  dplyr::mutate(
    bindingGroup = forcats::fct_relevel(.f = bindingGroup, "KERS")
  )


ggplotDf <- dplyr::select(
  plotData, geneId, bindingGroup, polIICols$exp
) %>% 
  tidyr::pivot_longer(
    cols = !c(geneId, bindingGroup),
    names_to = "sampleId",
    values_to = "polII_signal"
  ) %>% 
  dplyr::mutate(
    sampleId = forcats::fct_relevel(.f = sampleId, polII_ids),
    bindingGroup = forcats::fct_relevel(.f = bindingGroup, "KERS")
  )

ggpubr::compare_means(
  formula = polII_signal ~ bindingGroup, data = ggplotDf,
  method = "wilcox.test", paired = FALSE
)

ks.test(
  x = dplyr::filter(ggplotDf, bindingGroup == "KERS") %>% 
    dplyr::pull(polII_signal),
  y = dplyr::filter(ggplotDf, bindingGroup == "!!RS") %>% 
    dplyr::pull(polII_signal)
)

pt_violin <- ggplot2::ggplot(
  data = ggplotDf,
  mapping = aes(x = sampleId, y = log2(polII_signal),  group = bindingGroup)
) +
  geom_quasirandom(mapping = aes(color = bindingGroup), dodge.width = 1) +
  geom_boxplot(position = position_dodge(1), width=0.3, color = "black", fill = NA) +
  stat_compare_means(paired = FALSE, label = "p.format", size = 6) +
  guides(color = guide_legend(override.aes = list(size = 5) ) ) +
  labs(
    title = paste("Histone acetylation vs polII ChIPseq", polII_ids, sep = "\n", collapse = ""),
    y = "log2(polII-ChIPseq-FPKM)"
  ) +
  scale_color_manual(
    name = "binding",
    values = c("KERS" = "#E69F00", "!!RS" = "#56B4E9")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom"
  )



ggsave(
  filename = paste(outPrefix, ".violine.pdf", sep = ""),
  plot = pt_violin, width = 6, height = 8
)









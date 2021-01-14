suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(tidyverse))

## This script:
## variation of 48h vs 20h binding intensity for KERS ChIPseq w.r.t. polII 48h/20h fold change


rm(list = ls())


##################################################################################
## main configuration
comparisonName <- "TF_polII.48h_vs_20h"

file_polIIRatioConf <- here::here("data", "reference_data", "polII_ratio.config.tab")

outDir <- here::here("analysis", "05_KERS_48h_vs_20h", "combined_analysis", comparisonName)
outPrefix <- paste(outDir, "/", comparisonName, sep = "")

if(!dir.exists(outDir)) dir.create(path = outDir, recursive = TRUE)

polIIComparison <- "polII_untagged.48h_vs_20h"


tfDiffPairs <- list(
  p1 = list(
    name = "kdmB_48h_vs_20h",
    samples = c("An_kdmB_20h_HA_1", "An_kdmB_48h_HA_1")
  ),
  p2 = list(
    name = "rpdA_48h_vs_20h",
    samples = c("An_rpdA_20h_HA_1", "An_rpdA_48h_HA_1")
  ),
  p3 = list(
    name = "sntB_48h_vs_20h",
    samples = c("An_sntB_20h_HA_1", "An_sntB_48h_HA_1")
  ),
  p4 = list(
    name = "ecoA_48h_vs_20h",
    samples = c("An_ecoA_20h_HA_1", "An_ecoA_48h_HA_1")
  )
)

mainTfPair <- "p1"

tf1 <- tfDiffPairs[[mainTfPair]]$samples[1]
tf2 <- tfDiffPairs[[mainTfPair]]$samples[2]

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

## genes to read
geneSet <- suppressMessages(
  readr::read_tsv(file = file_genes, col_names = c("chr", "start", "end", "geneId", "score", "strand"))
) %>% 
  dplyr::select(geneId)

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
  dplyr::filter(comparison %in% polIIComparison)

polIIDiffPairs <- purrr::transpose(polIIRatioConf)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))

polII1 <- polIIDiffPairs[[polIIComparison]]$group1
polII2 <- polIIDiffPairs[[polIIComparison]]$group2

##################################################################################

## read the experiment sample details and select only those which are to be plotted

tfData <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = unique(map(tfDiffPairs, "samples") %>% unlist() %>% unname()),
  dataPath = TF_dataPath, profileMatrixSuffix = matrixType
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

expressionData %>% 
  dplyr::select(geneId, starts_with("hasPeak")) %>% 
  dplyr::group_by_at(.vars = vars(starts_with("hasPeak"))) %>% 
  dplyr::summarise(n = n())

##################################################################################
## peak data
peakExpDf <- dplyr::filter_at(
  .tbl = expressionData, .vars = unname(tfCols$hasPeak),
  .vars_predicate = any_vars(. == "TRUE")
) %>% 
  dplyr::filter_at(
    .vars = unname(polIICols$is_expressed[c(polII1, polII2)]),
    .vars_predicate = any_vars(. == TRUE)
  ) %>% 
  dplyr::select(
    geneId, starts_with("hasPeak"), starts_with("peakCoverage"),
    !!! map(polIIDiffPairs, function(x) as.name(x[["comparison"]])) %>% unname()
  )


## group genes into polII fold change bins
peakExpDf <- dplyr::mutate(
  peakExpDf,
  group = case_when(
    !! as.name(polIIDiffPairs[[polIIComparison]]$comparison) >= 2 ~ "G1: LFC >= 2",
    !! as.name(polIIDiffPairs[[polIIComparison]]$comparison) >= 1 ~ "G2: 1 <= LFC < 2",
    !! as.name(polIIDiffPairs[[polIIComparison]]$comparison) >= 0.5 ~ "G3: 0.5 <= LFC < 1",
    !! as.name(polIIDiffPairs[[polIIComparison]]$comparison) >= 0 ~ "G4: 0 <= LFC < 0.5",
    !! as.name(polIIDiffPairs[[polIIComparison]]$comparison) <= -2 ~ "G8: -2 >= LFC",
    !! as.name(polIIDiffPairs[[polIIComparison]]$comparison) <= -1 ~ "G7: -2 < LFC <= -1",
    !! as.name(polIIDiffPairs[[polIIComparison]]$comparison) <= -0.5 ~ "G6: -1 < LFC <= -0.5",
    !! as.name(polIIDiffPairs[[polIIComparison]]$comparison) < 0 ~ "G5: -0.5 < LFC < 0",
    TRUE ~ "0"
  )
) %>% 
  dplyr::group_by(group) %>% 
  dplyr::mutate(
    group = paste(group, "\n(n=", n(), ")", sep = "")
  ) %>% 
  dplyr::ungroup()


## box plots
longDf <- tidyr::pivot_longer(
  data = peakExpDf,
  cols = c(starts_with("hasPeak."), starts_with("peakCoverage.")),
  names_to = c(".value", "sampleId"),
  names_sep = "\\."
)

longDf$sampleId <- forcats::fct_relevel(.f = longDf$sampleId, names(tfCols$hasPeak))

## add TF, time info columns
longDf <- dplyr::left_join(
  x = longDf, y = dplyr::select(tfData, sampleId, TF, timepoint), by = "sampleId"
) %>% 
  dplyr::arrange(geneId, TF)

## make timepoint wise pair columns for hasPeak and peakCoverage at 20h and 48h
pairDf <- tidyr::pivot_wider(
  data = longDf,
  id_cols = c(geneId, TF, polII_untagged.48h_vs_20h, group),
  names_from = c(timepoint),
  values_from = c(hasPeak, peakCoverage),
  names_sep = "."
) %>% 
  dplyr::filter_at(.vars = vars(starts_with("hasPeak.")), .vars_predicate = all_vars(. == TRUE))

pltDf <- tidyr::pivot_longer(
  data = pairDf,
  cols = c(starts_with("hasPeak."), starts_with("peakCoverage.")),
  names_to = c(".value", "timepoint"),
  names_sep = "\\."
) %>% 
  dplyr::mutate(
    timepoint = forcats::fct_relevel(.f = timepoint, c("20h", "48h")),
    TF = forcats::fct_relevel(.f = TF, c("kdmB", "sntB", "ecoA", "rpdA")),
    group = forcats::fct_relevel(group, sort(unique(pairDf$group)))
  )


pt <- ggplot(data = pltDf, mapping = aes(x = group, y = peakCoverage)) +
  geom_boxplot(mapping = aes(fill = timepoint), alpha = 1, size = 0.5) +
  geom_point(mapping = aes(color = timepoint), shape = 16, size = 1, position=position_jitterdodge(0.2)) +
  scale_color_manual(values = c("20h" = "#b35806", "48h" = "#542788")) +
  scale_fill_manual(values = c("20h" = "#fee0b6", "48h" = "#d8daeb")) +
  # scale_y_continuous(trans = "log2") +
  labs(title = "KERS peak coverage comparison across the bins of 48h/20h ratios for untagged polII") +
  facet_wrap(. ~ TF, nrow = 1, ncol = 4, scales = "free_y") + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_blank(),
    strip.text.x = element_text(size = 14),
    legend.key.size = unit(1, "cm"),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold")
  )


pdf(file = paste(outPrefix, ".box.pdf", sep = ""), width = 20, height = 10)
pt
dev.off()



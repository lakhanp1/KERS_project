library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(scales)
library(ggplot2)


## This script compares the TF signal across polII fold change groups for TFs of interest


rm(list = ls())


##################################################################################
## main configuration
comparisonName <- "tf_signal_polII_correlation"
outPrefix <- here::here("kdmB_analysis", "combined_analysis", comparisonName)

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


## polII signal fold change pairs
polIIDiffPairs <- list(
  p1 = list(
    name = "untagged_48h_vs_untagged_20h_polII",
    title = "polII log2(untagged_48h \n vs untagged_20h)",
    samples = c("An_untagged_20h_polII_1", "An_untagged_48h_polII_1")
  ),
  p2 = list(
    name = "kdmB_del_48h_vs_kdmB_del_20h_polII",
    title = "polII log2(kdmB_del_48h \n vs kdmB_del_20h)",
    samples =  c("An_kdmB_del_20h_polII_1", "An_kdmB_del_48h_polII_1")
  ),
  p3 = list(
    name = "kdmB_del_20h_vs_untagged_20h_polII",
    title = "polII log2(kdmB_del_20h \n vs untagged_20h_polII)",
    samples = c("An_untagged_20h_polII_1", "An_kdmB_del_20h_polII_1")
  ),
  p4 = list(
    name = "kdmB_del_48h_vs_untagged_48h_polII",
    title = "polII log2(kdmB_del_48h \n vs untagged_48h_polII)",
    samples = c("An_untagged_48h_polII_1", "An_kdmB_del_48h_polII_1")
  )
)

mainPolIIPair <- "p1"

polII1 <- polIIDiffPairs[[mainPolIIPair]]$samples[1]
polII2 <- polIIDiffPairs[[mainPolIIPair]]$samples[2]

# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
matrixType <- "normalizedmatrix_5kb"
matrixDim = c(500, 200, 100, 10)

## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")


orgDb <- org.Anidulans.FGSCA4.eg.db

showExpressionHeatmap <- FALSE



##################################################################################

## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- suppressMessages(
  AnnotationDbi::select(x = orgDb, keys = geneSet$gene, columns = "DESCRIPTION", keytype = "GID")
)

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))

## gene information annotations: cluster and TF and polII expression values
geneInfo <- add_gene_info(file = file_geneInfo, clusterDf = geneSet)

head(geneInfo)

##################################################################################

## read the experiment sample details and select only those which are to be plotted

tfData <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = unique(map(tfDiffPairs, "samples") %>% unlist() %>% unname()),
                                 dataPath = TF_dataPath,
                                 matrixSource = matrixType)

polIIData <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = unique(map(polIIDiffPairs, "samples") %>% unlist() %>% unname()),
                                    dataPath = polII_dataPath,
                                    matrixSource = matrixType)

exptData <- dplyr::bind_rows(tfData, polIIData)

polII_ids <- exptData$sampleId[which(exptData$IP_tag == "polII")]
tfIds <- exptData$sampleId[which(exptData$IP_tag %in% c("HA", "MYC", "TAP") & exptData$TF != "untagged")]


polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)


tfCols <- sapply(c("hasPeak", "pval", "peakType", "tesPeakType", "peakCoverage", "peakDist", "summitDist", "upstreamExpr", "peakExpr", "relativeDist"),
                 FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
                 simplify = F, USE.NAMES = T)


expressionData <- get_polII_expressions(exptInfo = exptData,
                                        genesDf = geneInfo)

## add fold change columns
for (i in names(polIIDiffPairs)) {
  expressionData <- get_fold_change(df = expressionData,
                                    nmt = polIIDiffPairs[[i]]$samples[2],
                                    dmt = polIIDiffPairs[[i]]$samples[1],
                                    newCol = polIIDiffPairs[[i]]$name,
                                    isExpressedCols = polIICols$is_expressed)
}


expressionData <- get_TF_binding_data(exptInfo = tfData,
                                      genesDf = expressionData)


expressionData %>% 
  dplyr::select(gene, starts_with("hasPeak")) %>% 
  dplyr::group_by_at(.vars = vars(starts_with("hasPeak"))) %>% 
  dplyr::summarise(n = n())

##################################################################################
## peak data
peakExpDf <- dplyr::filter_at(.tbl = expressionData,
                              .vars = unname(tfCols$hasPeak),
                              .vars_predicate = any_vars(. == "TRUE")) %>% 
  dplyr::filter_at(.vars = unname(polIICols$is_expressed[c(polII1, polII2)]),
                   .vars_predicate = any_vars(. == TRUE)) %>% 
  dplyr::select(gene, starts_with("hasPeak"), starts_with("peakCoverage"),
                !!! map(polIIDiffPairs, function(x) as.name(x[["name"]])) %>% unname())


## group genes into polII fold change bins
peakExpDf <- dplyr::mutate(peakExpDf,
                           group = case_when(
                             !! as.name(polIIDiffPairs$p1$name) >= 2 ~ "G1: LFC >= 2",
                             !! as.name(polIIDiffPairs$p1$name) >= 1 ~ "G2: 1 <= LFC < 2",
                             !! as.name(polIIDiffPairs$p1$name) >= 0.5 ~ "G3: 0.5 <= LFC < 1",
                             !! as.name(polIIDiffPairs$p1$name) >= 0 ~ "G4: 0 <= LFC < 0.5",
                             !! as.name(polIIDiffPairs$p1$name) <= -2 ~ "G8: -2 >= LFC",
                             !! as.name(polIIDiffPairs$p1$name) <= -1 ~ "G7: -2 < LFC <= -1",
                             !! as.name(polIIDiffPairs$p1$name) <= -0.5 ~ "G6: -1 < LFC <= -0.5",
                             !! as.name(polIIDiffPairs$p1$name) < 0 ~ "G5: -0.5 < LFC < 0",
                             TRUE ~ "0"
                           ))

## convert dataframe to long format
longDf <- data.table::melt(
  data = as.data.table(peakExpDf),
  measure.vars = list(tfCols$hasPeak, tfCols$peakCoverage),
  variable.name = "sample",
  value.name = c("hasPeak", "peakCoverage")
)

longDf$sample <- forcats::fct_recode(
  .f = longDf$sample,
  structure(.Data = levels(longDf$sample), names = names(tfCols$hasPeak))
)

longDf$sample <- as.character(longDf$sample)

## add TF, time info columns
longDf <- dplyr::left_join(longDf, dplyr::select(tfData, sampleId, TF, timepoint), by = c("sample" = "sampleId")) %>% 
  as.data.table()

## make timepoint wise pair columns for hasPeak and peakCoverage
pairDf <- dcast.data.table(
  data = as.data.table(longDf),
  formula = gene + TF + untagged_48h_vs_untagged_20h_polII + group ~ timepoint,
  value.var = c("hasPeak", "peakCoverage")) %>% 
  as.tibble() %>% 
  dplyr::filter_at(.vars = vars(starts_with("hasPeak_")), .vars_predicate = all_vars(. == TRUE))


pltDf <- melt.data.table(
  data = as.data.table(pairDf),
  measure.vars = list(c("hasPeak_20h", "hasPeak_48h"), c("peakCoverage_20h", "peakCoverage_48h")),
  variable.name = "timepoint",
  value.name = c("hasPeak", "peakCoverage"))

pltDf$timepoint <- forcats::fct_recode(
  .f = pltDf$timepoint,
  structure(levels(pltDf$timepoint), names = c("20h", "48h"))
)


pltDf$TF <- factor(pltDf$TF, levels = c("kdmB", "sntB", "ecoA", "rpdA"))
pltDf$group <- factor(pltDf$group, levels = sort(unique(pltDf$group)))
# pltDf$timepoint <- as.character(pltDf$timepoint)


pt <- ggplot(data = pltDf, mapping = aes(x = group, y = peakCoverage)) +
  geom_boxplot(mapping = aes(fill = timepoint), alpha = 1, size = 0.5) +
  geom_point(mapping = aes(color = timepoint), shape = 16, size = 1, position=position_jitterdodge(0.2)) +
  scale_color_manual(values = c("20h" = "#b35806", "48h" = "#542788")) +
  scale_fill_manual(values = c("20h" = "#fee0b6", "48h" = "#d8daeb")) +
  # scale_y_continuous(trans = "log2") +
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


pdf(file = paste(outPrefix, "_box.pdf", sep = ""), width = 20, height = 10)
pt
dev.off()



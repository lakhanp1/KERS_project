suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(gghalves))
suppressPackageStartupMessages(library(ggpubr))

## KERS binding groups vs polII DEG ratio

rm(list = ls())

##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
timePoint <- "20h"
analysisName <- sprintf(fmt = "KERS_groups_vs_polII_DEG_ratios_%s", timePoint)

outDir <- here::here("analysis", "14_KERS_groups_vs_data", analysisName)
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_plotDegs <- paste(outDir, "/", "polII_ratio_ids.txt", sep = "")
file_kersBinding <- sprintf(
  fmt = here::here("analysis", "03_KERS_complex", "KERS_complex_%s",
                   "KERS_complex_%s.peaks_data.tab"),
  timePoint, timePoint
)


## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")
file_polIIRatioConf <- here::here("data", "reference_data", "polII_ratio.config.tab")

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db


##################################################################################
## polII ratio configuration
polIIDiffPairs <- suppressMessages(readr::read_tsv(file = file_plotDegs))
polIIRatioConf <- dplyr::left_join(
  x = polIIDiffPairs,
  y = suppressMessages(readr::read_tsv(file = file_polIIRatioConf)),
  by = "comparison"
)


polIIDiffPairs <- purrr::transpose(polIIRatioConf)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))


polII_ids <- unique(purrr::map(polIIDiffPairs, `[`, c("group2", "group1")) %>% unlist() %>% unname())

polIIData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = polII_ids,
  dataPath = polII_dataPath
) %>% 
  dplyr::mutate(
    TF = forcats::fct_relevel(TF, "untagged", "kdmB_del"),
    timepoint = forcats::fct_relevel(timepoint, "20h", "48h"),
  ) %>% 
  dplyr::arrange(timepoint, TF)

tfData <- NULL
inputData <- NULL
histData <- NULL
exptData <- dplyr::bind_rows(tfData, inputData, histData, polIIData)

exptDataList <- purrr::transpose(exptData)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)



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

geneInfo <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID"))


kersGenes <- suppressMessages(readr::read_tsv(file = file_kersBinding)) %>% 
  dplyr::select(geneId, bindingGroup)

geneInfo <- dplyr::left_join(x = geneInfo, y = kersGenes, by = "geneId") %>% 
  tidyr::replace_na(replace = list(bindingGroup = "0000"))

expressionData <- get_polII_expressions(
  exptInfo = polIIData, genesDf = geneInfo
)

glimpse(expressionData)


##################################################################################

rowId <- 1
## add fold change columns
for (rowId in names(polIIDiffPairs)) {
  expressionData <- get_fold_change(
    df = expressionData,
    nmt = polIIDiffPairs[[rowId]]$group1,
    dmt = polIIDiffPairs[[rowId]]$group2,
    newCol = polIIDiffPairs[[rowId]]$comparison,
    lfcLimit = 0,
    isExpressedCols = polIICols$is_expressed
  )
}

glimpse(expressionData)



##################################################################################
plotData <- dplyr::filter(
  expressionData, bindingGroup %in% c("KERS", "00RS")
) %>% 
  dplyr::mutate(
    bindingGroup = forcats::fct_relevel(.f = bindingGroup, "KERS")
  )


ggplotDf <- dplyr::select(
  plotData, geneId, bindingGroup, polIIRatioConf$comparison
) %>% 
  tidyr::pivot_longer(
    cols = !c(geneId, bindingGroup),
    names_to = "comparison",
    values_to = "lfc"
  ) %>% 
  dplyr::mutate(
    comparison = forcats::fct_relevel(.f = comparison, polIIRatioConf$comparison),
    bindingGroup = forcats::fct_relevel(.f = bindingGroup, "KERS")
  )

ggpubr::compare_means(
  formula = lfc ~ bindingGroup, data = ggplotDf, group.by = "comparison",
  method = "wilcox.test", paired = FALSE
)



pt_violin <- ggplot2::ggplot(
  data = ggplotDf,
  mapping = aes(x = comparison, y = lfc,  group = bindingGroup)
) +
  geom_quasirandom(mapping = aes(color = bindingGroup), dodge.width = 1) +
  geom_boxplot(
    position = position_dodge(1), width=0.3, 
    color = "black", fill = NA, outlier.shape = NA
  ) +
  stat_compare_means(paired = FALSE, label = "p.format", size = 6) +
  facet_grid(cols = vars(comparison), scales = "free") +
  guides(color = guide_legend(override.aes = list(size = 5) ) ) +
  labs(
    title = "KERS binding groups vs polII log2(mutant/WT)",
    y = "log2(polII-ChIPseq-FPKM)"
  ) +
  scale_color_manual(
    name = "binding",
    values = c("KERS" = "#E69F00", "00RS" = "#56B4E9")
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
  plot = pt_violin, width = 12, height = 8
)



##################################################################################




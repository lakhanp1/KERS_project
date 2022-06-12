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

source("D:/work_lakhan/github/omics_utils/02_RNAseq_scripts/s02_DESeq2_functions.R")

##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
timePoint <- "20h"
analysisName <- sprintf(fmt = "KERS_groups_vs_polII_DEseq2_%s", timePoint)

outDir <- here::here("analysis", "14_KERS_groups_vs_data", analysisName)
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_degIds <- paste(outDir, "/", "DESeq2_ids.txt", sep = "")
file_kersBinding <- sprintf(
  fmt = here::here("analysis", "03_KERS_complex", "KERS_complex_%s",
                   "KERS_complex_%s.peaks_data.tab"),
  timePoint, timePoint
)

diffDataPath <- here::here("analysis", "02_polII_DEGs")
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_conf.txt")

## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")
file_polIIRatioConf <- here::here("data", "reference_data", "polII_ratio.config.tab")

orgDb <- org.Anidulans.FGSCA4.eg.db


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


degIds <- suppressMessages(readr::read_tsv(file = file_degIds))

rnaseqInfo <- get_diff_info(degInfoFile = file_RNAseq_info, dataPath = diffDataPath) %>% 
  dplyr::filter(comparison %in% degIds$comparison)

## get DEG data for each SM_TF cluster from its own polII_DEG set
rowId <- 1
degData <- NULL


for (rowId in 1:nrow(rnaseqInfo)) {
  
  tmpDf <- suppressMessages(readr::read_tsv(file = rnaseqInfo$deg[rowId])) %>% 
    dplyr::select(geneId, log2FoldChange, shrinkLog2FC, pvalue, padj) %>% 
    dplyr::mutate(
      comparison = rnaseqInfo$comparison[rowId]
    )
  
  subData <- dplyr::left_join(
    x = geneInfo, y = tmpDf, by = "geneId"
  )
  
  degData <- dplyr::bind_rows(degData, subData)
  
}


##################################################################################
plotData <- dplyr::filter(
  degData, bindingGroup %in% c("KERS", "00RS")
) %>% 
  dplyr::mutate(
    bindingGroup = forcats::fct_relevel(.f = bindingGroup, "KERS")
  )


ggplotDf <- dplyr::mutate(
  plotData,
  comparison = forcats::fct_relevel(.f = comparison, degIds$comparison),
  bindingGroup = forcats::fct_relevel(.f = bindingGroup, "KERS")
)

ggpubr::compare_means(
  formula = log2FoldChange ~ bindingGroup, data = ggplotDf, group.by = "comparison",
  method = "wilcox.test", paired = FALSE
)



pt_violin <- ggplot2::ggplot(
  data = ggplotDf,
  mapping = aes(x = comparison, y = log2FoldChange,  group = bindingGroup)
) +
  geom_quasirandom(mapping = aes(color = bindingGroup), dodge.width = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_boxplot(
    position = position_dodge(1), width=0.3, 
    color = "black", fill = NA, outlier.shape = NA
  ) +
  stat_compare_means(paired = FALSE, label = "p.format", size = 6, label.y = 4.5) +
  coord_cartesian(
    ylim = c(-5, 5)
  ) +
  facet_grid(cols = vars(comparison), scales = "free") +
  guides(color = guide_legend(override.aes = list(size = 5) ) ) +
  labs(
    title = "KERS binding groups vs polII DESeq2 log2(mutant/WT)",
    y = "log2(mutant/WT) from DESeq2"
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




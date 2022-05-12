suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggpubr))

## log2(fold-change) heatmaps for mutant/WT

rm(list = ls())

##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
analysisName <- "geneset1"
workDir <- here::here("analysis", "13_geneset_analysis", analysisName)
outPrefix <- paste(workDir, "/", analysisName, sep = "")

file_plotDegs <- paste(workDir, "/", "polII_ratio_ids.txt", sep = "")
file_geneSubset <- paste(workDir, "/", "geneList.txt", sep = "")


## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")
file_polIIRatioConf <- here::here("data", "reference_data", "polII_ratio.config.tab")

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF


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
geneSet <- data.table::fread(
  file = file_genes, header = F,
  col.names = c("chr", "start", "end", "geneId", "score", "strand"))

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, keytype = "GID",
                                  columns = c("GENE_NAME", "DESCRIPTION")) %>% 
  dplyr::rename(geneId = GID)

geneInfo <- geneDesc

expressionData <- get_polII_expressions(
  exptInfo = polIIData, genesDf = geneInfo
)

glimpse(expressionData)

geneSubset <- suppressMessages(
  readr::read_tsv(file = file_geneSubset, comment = "#")
)



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


geneSubset <- dplyr::left_join(
  x = geneSubset, y = expressionData, by = "geneId"
)

##################################################################################

## polII signal matrix
polIIMat <- as.matrix(log2(expressionData[, polII_ids] + 1))
rownames(polIIMat) <- expressionData$geneId

quantile(polIIMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.997, 0.999, 0.9999, 1), na.rm = T)


polII_color <- colorRamp2(
  breaks = c(0, quantile(polIIMat, c(0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99, 0.995, 0.999))),
  colors = c("white", RColorBrewer::brewer.pal(n = 9, name = "RdPu"))
)


lfc_color <- colorRamp2(
  breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2),
  colors = rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))
)

##################################################################################
genesetPolIIMat <- dplyr::select(.data = geneSubset, geneName, !!!polII_ids) %>% 
  tibble::column_to_rownames(var = "geneName") %>% 
  as.matrix() %>% 
  log2()


ht_signal <- Heatmap(
  matrix = genesetPolIIMat,
  col = polII_color,
  column_split = stringr::str_replace(
    string = colnames(genesetPolIIMat), pattern = "An.*_(\\d+h)_polII_\\d", replacement = "\\1"
  ),
  column_gap = unit(c(4), "mm"),
  # column_title = "polII ChIPseq signal",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  # cell_fun = function(j, i, x, y, width, height, fill) {
  #   grid.text(sprintf("%.1f", rep1Mat[i, j]), x, y, gp = gpar(fontsize = 10))
  # },
  column_labels = stringr::str_replace(
    string = colnames(genesetPolIIMat), pattern = "An_(.*)_(\\d+h)_polII_\\d", replacement = "\\1"
  ),
  heatmap_legend_param = list(
    title = "log2(polII signal)",
    # at = colAt,
    legend_height = unit(4, "cm"),
    labels_gp = gpar(fontsize = 12),
    title_gp = gpar(fontsize = 12, fontface = "bold")
  ),
  # heatmap_width = unit(5, "cm"),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  row_names_gp = gpar(fontsize = 16, fontface = "bold"),
  row_names_max_width = unit(20, "cm")
)



genesetLfcMat <- dplyr::select(.data = geneSubset, geneName, starts_with("cmp.")) %>% 
  tibble::column_to_rownames(var = "geneName") %>% 
  as.matrix()

ht_lfc <- Heatmap(
  matrix = genesetLfcMat,
  col = lfc_color,
  column_split = stringr::str_replace(
    string = colnames(genesetLfcMat), pattern = "cmp.*\\.(\\d+h)", replacement = "\\1"
  ),
  column_gap = unit(c(4), "mm"),
  # column_title = "mutant / WT",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  # cell_fun = function(j, i, x, y, width, height, fill) {
  #   grid.text(sprintf("%.1f", lfcMat[i, j]), x, y, gp = gpar(fontsize = 10))
  # },
  column_labels = stringr::str_replace(
    string = colnames(genesetLfcMat), pattern = "cmp.(.*)\\.\\d+h", replacement = "\\1"
  ),
  heatmap_legend_param = list(
    title = "log2(fold change)",
    legend_height = unit(4, "cm"),
    labels_gp = gpar(fontsize = 12),
    title_gp = gpar(fontsize = 12, fontface = "bold")
  ),
  # heatmap_width = unit(5, "cm"),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  row_names_gp = gpar(fontsize = 16, fontface = "bold"),
  row_names_max_width = unit(10, "cm")
)


htList <- ht_lfc + ht_signal

pdf(file = paste(outPrefix, ".lfc_heatmap.pdf", sep = ""), width = 6, height = 8)
draw(
  object = ht_lfc,
  column_title = analysisName,
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
)
dev.off()


##################################################################################




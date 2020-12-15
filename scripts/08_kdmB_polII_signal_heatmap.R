suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))

## this script plots the polII signal heatmap and polII fold change heatmap
## for the kdmB complex members

rm(list = ls())

##################################################################################
analysisName <- "kdmB_del_sntB_del_polII"
outDir <- here::here("analysis", "02_polII_signal")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_polIIsamples <-paste(outDir, "/kdmB_del_sntB_del_samples.txt", sep = "")
file_factors <- here::here("data", "reference_data", "factor_names.list")

file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db

file_factorSignal <- paste(polII_dataPath, "/factors_polII_signal.tab", sep = "")
file_factorLfc <- paste(polII_dataPath, "/factors_polII_signal_LFC.tab", sep = "")

##################################################################################
exptInfo <- suppressMessages(readr::read_tsv(file = file_exptInfo))

polIISamples <- suppressMessages(readr::read_tsv(file = file_polIIsamples)) %>% 
  dplyr::left_join(y = exptInfo, by = "sampleId")


factorSignal <- suppressMessages(readr::read_tsv(file = file_factorSignal)) %>% 
  dplyr::select(sample, kdmB, rpdA, sntB, ecoA)

dt <- dplyr::left_join(x = polIISamples, y = factorSignal, by = c("sampleId" = "sample")) %>% 
  dplyr::mutate(sampleId = gsub(pattern = "_[12]$", replacement = "", x = sampleId, perl = T))


##################################################################################
## replicate 1 polII signal heatmap

rep1Df <- dplyr::filter(dt, rep == 1) %>% 
  dplyr::select(sampleId, colnames(factorSignal)[-1]) %>% 
  dplyr::mutate_at(.vars = vars(!!!colnames(factorSignal)[-1]),
                   .funs = funs(if_else(condition = . < 1, true = 1, false = .)))


rep1Mat <- log2(as.matrix(tibble::column_to_rownames(rep1Df, "sampleId")))


quantile(rep1Mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)

rep1Col <- colorRamp2(breaks = seq(0, quantile(rep1Mat, 0.999), length.out = 5),
                      colors = RColorBrewer::brewer.pal(n = 5, name = "RdPu"))


rep1Ht <- Heatmap(
  matrix = rep1Mat,
  col = rep1Col,
  column_title = "polII signal for deletion mutants of kdmB complex",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  # cell_fun = function(j, i, x, y, width, height, fill) {
  #   grid.text(sprintf("%.1f", rep1Mat[i, j]), x, y, gp = gpar(fontsize = 10))
  # },
  heatmap_legend_param = list(
    title = "log2(polII signal)",
    # at = colAt,
    legend_height = unit(4, "cm"),
    labels_gp = gpar(fontsize = 12),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    title_position = "topcenter"
  ),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  row_names_gp = gpar(fontsize = 16),
  row_names_max_width = unit(20, "cm")
)


pdf(file = paste(outPrefix, "_kdmB_rep1_FPKM_heatmap.pdf", sep = ""), width = 8, height = 6)
draw(rep1Ht)
dev.off()

##################################################################################
## polII signal fold change heatmap

lfcDataAll <- suppressMessages(readr::read_tsv(file = file_factorLfc)) %>% 
  dplyr::select(lfcCol, sampleId1, sampleId2, kdmB, rpdA, sntB, ecoA)

lfcDf <- dplyr::left_join(x = polIISamples, y = lfcDataAll, by = c("sampleId" = "sampleId1")) %>% 
  dplyr::filter(! is.na(lfcCol)) %>% 
  dplyr::select(!! colnames(lfcDataAll)[-2:-3])


lfcMat <- tibble::column_to_rownames(lfcDf, var = "lfcCol") %>% 
  as.matrix()

rownames(lfcMat) <- gsub(pattern = "(lfc.|An_|_polII)", replacement = "", x = rownames(lfcMat)) %>% 
  gsub(pattern = "_vs_", replacement = " / ", x = .)


lfc_color <- colorRamp2(
  breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2),
  # colors = RColorBrewer::brewer.pal(n = 7, name = "PuOr"),
  colors = c("#b35806", "#f1a340", "#f7f7f7", "#f7f7f7", "#f7f7f7", "#998ec3", "#542788")
)


htLfc <- Heatmap(
  matrix = lfcMat,
  col = lfc_color,
  column_title = "polII signal fold change for deletion mutants of kdmB complex",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  # cell_fun = function(j, i, x, y, width, height, fill) {
  #   grid.text(sprintf("%.1f", lfcMat[i, j]), x, y, gp = gpar(fontsize = 10))
  # },
  heatmap_legend_param = list(
    title = "log2(fold change)",
    legend_height = unit(4, "cm"),
    labels_gp = gpar(fontsize = 12),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    title_position = "topcenter"
  ),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_side = "left",
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  row_names_gp = gpar(fontsize = 12),
  row_names_max_width = unit(10, "cm")
)


pdf(file = paste(outPrefix, "_lfc_heatmap.pdf", sep = ""), width = 8, height = 6)
draw(htLfc)
dev.off()






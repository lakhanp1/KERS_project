suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(GGally))

## this script performs following tasks
## 1) extract the polII signal for all samples
## 2) prepare polII signal matrix for gene of interest (kdmB complex members) and plot heatmap
## 3) calculate the log2 fold change for each polII sample using appropriate control data


rm(list = ls())

##################################################################################
analysisName <- "polII_signal"
outPrefix <- here::here("..", "data", "A_nidulans", "polII_data", analysisName)

file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db

file_factors <- here::here("data", "reference_data", "factor_names.list")

file_polIIsamples <- paste(polII_dataPath, "/sample_polII.list", sep = "")
file_polIICtrlPairs <- paste(polII_dataPath, "/polII_sample_control_pairs.txt", sep = "")
file_polIITimeDiffPairs <- paste(polII_dataPath, "/polII_48h_vs_20h_diff_pairs.txt", sep = "")
file_polIIMat <- paste(polII_dataPath, "/polII_signal_matrix.tab", sep = "")

##################################################################################
## extract polII signal matrix for all the genes
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "geneId", "score", "strand")) %>% 
  dplyr::select(-score)

geneDesc <- select(x = orgDb, keys = geneSet$geneId, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID"))

polIISamples <- suppressMessages(readr::read_tsv(file = file_polIIsamples, col_names = "id"))

polII_info <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = polIISamples$id,
  dataPath = polII_dataPath
)

polIICols <- list(
  exp = structure(polIISamples$id, names = polIISamples$id),
  is_expressed = structure(paste("is_expressed", ".", polIISamples$id, sep = ""), names = polIISamples$id)
)


polIIMat <- suppressMessages(readr::read_tsv(file = file_polIIMat))


##################################################################################
## polII signal matrix for genes of interest

goi <- suppressMessages(readr::read_tsv(file = file_factors))

goiPolII <- dplyr::left_join(goi, polIIMat, by = "geneId")


exprDf <- dplyr::select(goiPolII, -geneId, -chr, -start, -end, -strand, -length, -DESCRIPTION,
                        -starts_with("is_expressed")) %>% 
  tidyr::gather(key = sample, value = expression, -name, factor_key = TRUE) %>% 
  tidyr::spread(key = name, value = expression) %>% 
  dplyr::mutate_if(.predicate = is.factor, .funs = as.character) %>% 
  dplyr::select(sample, goi$name) %>% 
  as.data.frame()

readr::write_tsv(x = exprDf, file = paste(polII_dataPath, "/factors_polII_signal.tab", sep = ""))

## for log2 calculations, set the values to 1 if they are < 1
exprDf <- dplyr::mutate_at(
  .tbl = exprDf,
  .vars = vars(!!!goi$name),
  .funs = list(~ if_else(condition = . < 1, true = 1, false = .))
)


exprMat <- log2(as.matrix(exprDf[, goi$name]))
rownames(exprMat) <- exprDf$sample

quantile(exprMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)

htCol <- colorRamp2(breaks = quantile(exprMat, c(0, 0.1, 0.9, 0.99, 0.999)),
                    colors = c("white", "#f7f4f9", "#ce1256", "#980043", "#67001f"))

# htCol <- colorRamp2(breaks = seq(0, quantile(exprMat, 0.999), length.out = 8),
#            colors = RColorBrewer::brewer.pal(n = 8, name = "RdPu"))

colAt <- as.numeric(sprintf("%.0f", c(seq(quantile(exprMat, 0.1), quantile(exprMat, 0.99), length = 4),
                                      quantile(exprMat, 0.995))))

ht <- Heatmap(
  matrix = exprMat,
  col = htCol,
  column_title = "polII signal for various factors in mutants",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", exprMat[i, j]), x, y, gp = gpar(fontsize = 8))
  },
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
  row_names_gp = gpar(fontsize = 12),
  row_names_max_width = unit(10, "cm")
)


pdf(file = paste(polII_dataPath, "/factors_polII_signal.pdf", sep = ""), width = 8, height = 10)
draw(ht)
dev.off()


##################################################################################
## polII data fold change results: mutant-wt pairs
polIICtrlPairs <- suppressMessages(readr::read_tsv(file = file_polIICtrlPairs)) %>% 
  dplyr::mutate(lfcCol = paste("lfc.", sampleId1, "_vs_", sampleId2, sep = ""))

lfcMat <- polIIMat

# i <- 1

for (i in 1:nrow(polIICtrlPairs)) {
  
  lfcMat <- get_fold_change(
    df = lfcMat,
    nmt = polIICtrlPairs$sampleId1[i],
    dmt = polIICtrlPairs$sampleId2[i],
    newCol = polIICtrlPairs$lfcCol[i],
    isExpressedCols = polIICols$is_expressed
  )
  
}


lfcMat <- dplyr::select(lfcMat, geneId, chr, start, end, strand, length, DESCRIPTION, starts_with("lfc."))

readr::write_tsv(x = lfcMat, file = paste(polII_dataPath, "/polII_LFC_mut_vs_wt.tab", sep = ""))

pt <- ggpairs(data = dplyr::select(lfcMat, starts_with("lfc.")),
              upper = list(continuous = wrap("points", size = 0.1)),
              lower = list(continuous = wrap("cor", size = 4)),
              diag = list(continuous = "densityDiag")) +
  theme_bw() +
  theme(
    strip.text.y = element_text(size = 8, angle = 0, hjust = 0),
    strip.text.x = element_text(size = 8, angle = 90, hjust = 0)
  )

png(filename = paste(polII_dataPath, "/polII_LFC_mut_vs_wt.png", sep = ""),
    width = 6000, height = 6000, res = 250)

pt
dev.off()


##################################################################################
## polII data fold change results: 48h-20h pairs
polIITimeDiffPairs <- suppressMessages(readr::read_tsv(file = file_polIITimeDiffPairs)) %>% 
  dplyr::mutate(lfcCol = paste("lfc.", sampleId1, "_vs_", sampleId2, sep = ""))

lfcMat2 <- polIIMat

# i <- 1

for (i in 1:nrow(polIITimeDiffPairs)) {
  
  lfcMat2 <- get_fold_change(
    df = lfcMat2,
    nmt = polIITimeDiffPairs$sampleId1[i],
    dmt = polIITimeDiffPairs$sampleId2[i],
    newCol = polIITimeDiffPairs$lfcCol[i],
    isExpressedCols = polIICols$is_expressed
  )
  
}


lfcMat2 <- dplyr::select(lfcMat2, geneId, chr, start, end, strand, length, DESCRIPTION, starts_with("lfc."))

readr::write_tsv(x = lfcMat2, file = paste(polII_dataPath, "/polII_LFC_40h_vs_20h.tab", sep = ""))

pt2 <- ggpairs(data = dplyr::select(lfcMat2, starts_with("lfc.")),
               upper = list(continuous = wrap("points", size = 0.1)),
               lower = list(continuous = wrap("cor", size = 4)),
               diag = list(continuous = "densityDiag")) +
  theme_bw() +
  theme(
    strip.text.y = element_text(size = 8, angle = 0, hjust = 0),
    strip.text.x = element_text(size = 8, angle = 90, hjust = 0)
  )

png(filename = paste(polII_dataPath, "/polII_LFC_40h_vs_20h.png", sep = ""),
    width = 6000, height = 6000, res = 250)

pt2
dev.off()

##################################################################################
## fold change heatmap for genes of interest
goi <- read_tsv(file = file_factors)

goiLfc <- dplyr::left_join(goi, lfcMat, by = "geneId")


lfcDf <- dplyr::select(goiLfc, name, starts_with("lfc.")) %>% 
  tidyr::gather(key = lfcCol, value = expression, -name, factor_key = TRUE) %>% 
  tidyr::spread(key = name, value = expression) %>% 
  dplyr::mutate_if(.predicate = is.factor, .funs = as.character) %>% 
  dplyr::left_join(y = polIICtrlPairs, by = c("lfcCol" = "lfcCol")) %>% 
  dplyr::select(lfcCol, sampleId1, sampleId2, goi$name) %>% 
  as.data.frame()


fwrite(x = lfcDf, file = paste(polII_dataPath, "/factors_polII_signal_LFC.tab", sep = ""),
       sep = "\t", col.names = T, quote = F, row.names = F)


lfcMat <- dplyr::select(lfcDf, -sampleId1, -sampleId2) %>% 
  tibble::column_to_rownames(var = "lfcCol") %>% 
  as.matrix()

rownames(lfcMat) <- gsub(pattern = "(lfc.|An_|_polII_\\d)", replacement = "", x = rownames(lfcMat)) %>% 
  gsub(pattern = "_vs_", replacement = " / ", x = .)

lfc_color <- colorRamp2(breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2),
                        # colors = RColorBrewer::brewer.pal(n = 7, name = "PuOr"),
                        colors = c("#b35806", "#f1a340", "#f7f7f7", "#f7f7f7", "#f7f7f7", "#998ec3", "#542788"))


htLfc <- Heatmap(
  matrix = lfcMat,
  col = lfc_color,
  column_title = "polII signal fold change for various factors in mutants",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
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


pdf(file = paste(polII_dataPath, "/factors_polII_signal_LFC.pdf", sep = ""), width = 10, height = 12)
draw(htLfc)
dev.off()





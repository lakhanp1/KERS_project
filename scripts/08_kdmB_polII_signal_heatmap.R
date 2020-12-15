library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)

## this script plots the polII signal heatmap and polII fold change heatmap
## for the kdmB complex members

rm(list = ls())

##################################################################################
analysisName <- "kdmB_del_sntB_del_polII"
outPrefix <- here::here("kdmB_analysis", "polII_signal", analysisName)

file_polIIsamples <- here::here("kdmB_analysis", "polII_signal", "kdmB_del_sntB_del_samples.txt")

## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")


orgDb <- org.Anidulans.FGSCA4.eg.db

file_factors <- here::here("data", "referenceData/factor_names.list")

orgDb <- org.Anidulans.FGSCA4.eg.db

file_factorSignal <- paste(polII_dataPath, "/factors_polII_signal.tab", sep = "")
file_factorLfc <- paste(polII_dataPath, "/factors_polII_signal_LFC.tab", sep = "")

geneCdsFile <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_version_s10-m04-r03_CDS_Unique.bed"


##################################################################################
exptInfo <- readr::read_tsv(file = file_exptInfo, col_names = T, na = character())

polIISamples <- fread(file = file_polIIsamples, sep = "\t", header = F,
                      stringsAsFactors = F, col.names = c("id"), data.table = F) %>% 
  dplyr::left_join(y = exptInfo, by = c("id" = "sampleId"))


factorSignal <- readr::read_tsv(file = file_factorSignal, col_names = T) %>% 
  dplyr::select(sample, kdmB, rpdA, sntB, ecoA)

dt <- dplyr::left_join(x = polIISamples, y = factorSignal, by = c("id" = "sample")) %>% 
  dplyr::mutate(sample = gsub(pattern = "_[12]$", replacement = "", x = id, perl = T))


##################################################################################
## replicate 1 polII signal heatmap

rep1Df <- dplyr::filter(dt, rep == 1) %>% 
  dplyr::select(id, colnames(factorSignal)[-1]) %>% 
  dplyr::mutate_at(.vars = vars(!!!colnames(factorSignal)[-1]),
                   .funs = funs(if_else(condition = . < 1, true = 1, false = .)))


rep1Mat <- log2(as.matrix(tibble::column_to_rownames(rep1Df, "id")))


quantile(rep1Mat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)

rep1Col <- colorRamp2(breaks = seq(0, quantile(rep1Mat, 0.999), length.out = 5),
                    colors = RColorBrewer::brewer.pal(n = 5, name = "RdPu"))


rep1Ht <- Heatmap(matrix = rep1Mat,
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
              row_names_gp = gpar(fontsize = 12),
              row_names_max_width = unit(20, "cm")
)


pdf(file = paste(outPrefix, "_kdmB_rep1_FPKM_heatmap.pdf", sep = ""), width = 8, height = 6)
draw(rep1Ht)
dev.off()

##################################################################################
## polII signal fold change heatmap

lfcDataAll <- readr::read_tsv(file = file_factorLfc, col_names = T) %>% 
  dplyr::select(lfcCol, sampleId, control, kdmB, rpdA, sntB, ecoA)


lfcDf <- dplyr::left_join(x = polIISamples, y = lfcDataAll, by = c("id" = "sampleId")) %>% 
  dplyr::filter(! is.na(lfcCol)) %>% 
  dplyr::select(!! colnames(lfcDataAll)[-2:-3])


lfcMat <- tibble::column_to_rownames(lfcDf, var = "lfcCol") %>% 
  as.matrix()

rownames(lfcMat) <- gsub(pattern = "(lfc.|An_|_polII)", replacement = "", x = rownames(lfcMat)) %>% 
  gsub(pattern = "_vs_", replacement = " / ", x = .)


lfc_color <- colorRamp2(breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2),
                        # colors = RColorBrewer::brewer.pal(n = 7, name = "PuOr"),
                        colors = c("#b35806", "#f1a340", "#f7f7f7", "#f7f7f7", "#f7f7f7", "#998ec3", "#542788"))


htLfc <- Heatmap(matrix = lfcMat,
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






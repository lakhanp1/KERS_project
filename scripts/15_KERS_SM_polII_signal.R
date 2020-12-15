library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(esquisse)
library(summarytools)
library(here)
library(ggpubr)



rm(list = ls())

##################################################################################
analysisName <- "KERS_20h_SM"
outDir <- here::here("kdmB_analysis", "SM_analysis", analysisName)

outPrefix <- paste(outDir, "/", analysisName, sep = "")

tfIds <- c("An_kdmB_20h_HA_1", "An_ecoA_20h_HA_1", "An_rpdA_20h_HA_1", "An_sntB_20h_HA_1")
polII_ids <- c("An_untagged_20h_polII_1")

file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")

orgDb <- org.Anidulans.FGSCA4.eg.db

if(!dir.exists(outDir)){
  dir.create(outDir)
}

##################################################################################
## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, columns = c("SM_CLUSTER"), keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))

##################################################################################

tfInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = tfIds,
                                 dataPath = TF_dataPath,
                                 matrixSource = "normalizedmatrix")

polII_info <- get_sample_information(exptInfoFile = file_exptInfo,
                                     samples = polII_ids,
                                     dataPath = polII_dataPath,
                                     matrixSource = "normalizedmatrix")

exptData <- dplyr::bind_rows(tfInfo, polII_info)
exptDataList <- purrr::transpose(exptData) %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

polIICols <- list(
  exp = structure(polII_info$sampleId, names = polII_info$sampleId),
  is_expressed = structure(paste("is_expressed", ".", polII_info$sampleId, sep = ""), names = polII_info$sampleId)
)

tfCols <- sapply(
  X = c("peakDist", "featureCovFrac", "hasPeak", "peakCoverage", "peakPosition", "peakId", "peakType",
        "peakPval", "peakEnrichment", "preference", "peakCategory"),
  FUN = function(x){ structure(paste(x, ".", tfInfo$sampleId, sep = ""), names = tfInfo$sampleId) },
  simplify = F, USE.NAMES = T)

##################################################################################
## prepare data for plotting
tfBindingMat <- peak_target_matrix(sampleInfo = tfInfo, position = "best")
polIIData <- get_polII_expressions(genesDf = geneSet, exptInfo = polII_info) %>% 
  dplyr::mutate_at(.vars = polIICols$exp, .funs = list(~log2(pmax(., 1))))


quantile(unlist(polIIData[polIICols$exp]),
         c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.997, 0.999, 0.9999, 1), na.rm = T)


polIIColor <- colorRamp2(
  breaks = c(0, quantile(unlist(polIIData[polIICols$exp]),
                         c(0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99, 0.995, 0.999)) ),
  colors = c("white", RColorBrewer::brewer.pal(n = 9, name = "RdPu")))


mergedData <- dplyr::left_join(x = polIIData, y = tfBindingMat, by = "gene") %>% 
  dplyr::mutate_at(.vars = vars(starts_with("hasPeak.")), .funs = list(~if_else(is.na(.), FALSE, .))) %>% 
  dplyr::filter(!is.na(SM_CLUSTER)) %>% 
  dplyr::mutate(
    SM_CLUSTER = gsub(pattern = "SM_cluster_", replacement = "", x = SM_CLUSTER, fixed = TRUE)
  ) %>% 
  dplyr::group_by(SM_CLUSTER) %>% 
  dplyr::mutate(index = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(gene, SM_CLUSTER, index, starts_with("hasPeak"), starts_with("peakPval"),
                unname(polIICols$exp), unname(polIICols$is_expressed))



##################################################################################
ptTheme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        strip.text.y = element_text(hjust = 0.5, size = 14, face = "bold", angle = 180),
        strip.background = element_rect(fill="white", size = 0.2),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = unit(rep(0.5, 4), "cm"))


##################################################################################
## heatbox with KERS targets

# mergedData <- dplyr::filter(mergedData, SM_CLUSTER %in% c("01", "02", "03", "04", "05"))
plotTitle <- paste("SM genes | KERS targets |", polII_ids)
## prepare peak data
hasPeak <- as.data.table(mergedData) %>% 
  data.table::melt(id.vars = c("gene", "SM_CLUSTER", "index"),
                   measure.vars = list(hasPeak = tfCols$hasPeak, peakPval = tfCols$peakPval),
                   variable.name = "sampleId") %>% 
  as_tibble() %>% 
  dplyr::filter(hasPeak == TRUE)

hasPeak$sampleId <- forcats::fct_recode(
  .f = hasPeak$sampleId,
  structure(levels(hasPeak$sampleId), names = tfInfo$sampleId)
)

tfPointPosition <- structure(seq(0, 1, length.out = length(tfInfo$sampleId)+2)[2:(1+length(tfInfo$sampleId))],
                             names = tfInfo$sampleId)



hasPeak$tfPosition <- tfPointPosition[hasPeak$sampleId]

## prepare polII signal data
polIISignal <- as.data.table(mergedData) %>% 
  data.table::melt(id.vars = c("gene", "SM_CLUSTER", "index"),
                   measure.vars = list(polII = polIICols$exp, is_expressed = polIICols$is_expressed),
                   variable.name = "sampleId") %>% 
  as_tibble() %>% 
  dplyr::mutate(polII = if_else(is_expressed == FALSE, 0, polII))

polIISignal$sampleId <- forcats::fct_recode(
  .f = polIISignal$sampleId,
  structure(levels(polIISignal$sampleId), names = polII_info$sampleId)
)

tfColor <- purrr::map_chr(.x = exptDataList[tfInfo$sampleId],
               .f = function(x) unlist(strsplit(x = x$color, split = ","))[2])

pt <- ggplot() +
  geom_tile(data = polIISignal, mapping = aes(x = 0.5, y = sampleId, fill = polII)) +
  geom_point(data = hasPeak,
             mapping = aes(x = tfPosition, y = "tf", color = sampleId),
             size = 2) +
  scale_fill_gradientn(
    name = "log2(polII-signal)",
    colours = c("white", RColorBrewer::brewer.pal(n = 9, name = "RdPu")),
    values = rescale(
      c(0, quantile(unlist(polIIData[polIICols$exp]),
                    c(0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99, 0.995, 0.999), names =F)
      )
    )
  ) +
  scale_color_manual(name = "KERS peak", values = tfColor) +
  scale_x_continuous(expand = expand_scale(add = c(0.0, 0.0))) +
  ggtitle(plotTitle) +
  facet_wrap(facets =  vars(gene), ncol = 20, strip.position = "left", dir = "v") +
  ptTheme + theme(strip.text.y = element_blank())


pdf(file = paste(outPrefix, ".tile_plot.pdf", sep = ""), width = 10, height = 10)
print(pt)
dev.off()










suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggpubr))

## Tile plot showing KERS binding and polII FPKM signal for each of the SM clusters gene

rm(list = ls())

##################################################################################
analysisName <- "polII_FPKM_KERS_binding"
outDir <- here::here("analysis", "07_SM_analysis", "KERS_20h_SM")

outPrefix <- paste(outDir, "/", analysisName, sep = "")

tfIds <- c("An_kdmB_20h_HA_1", "An_ecoA_20h_HA_1", "An_rpdA_20h_HA_1", "An_sntB_20h_HA_1")
polII_ids <- c("An_untagged_20h_polII_1")

## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/annotation_resources/geneid2go.ANidulans.topGO.map"

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db

# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
matrixType <- "normalizedmatrix_5kb"
matrixDim = c(500, 200, 100, 10)

##################################################################################
if(!dir.exists(outDir)) dir.create(path = outDir, recursive = TRUE)

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

##################################################################################

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfIds,
  dataPath = TF_dataPath
)

polII_info <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = polII_ids,
  dataPath = polII_dataPath
)

exptInfo <- dplyr::bind_rows(tfInfo, polII_info)
exptInfoList <- purrr::transpose(exptInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))


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

##################################################################################
## prepare data for plotting
tfBindingMat <- peak_target_matrix(sampleInfo = tfInfo, position = "best")
polIIData <- get_polII_expressions(genesDf = geneInfo, exptInfo = polII_info)

quantile(unlist(polIIData[polIICols$exp]),
         c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.997, 0.999, 0.9999, 1), na.rm = T)


polIIColor <- colorRamp2(
  breaks = c(0, quantile(unlist(polIIData[polIICols$exp]),
                         c(0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99, 0.995, 0.999)) ),
  colors = c("white", RColorBrewer::brewer.pal(n = 9, name = "RdPu")))


mergedData <- dplyr::left_join(x = polIIData, y = tfBindingMat, by = "geneId") %>% 
  dplyr::mutate_at(.vars = vars(starts_with("hasPeak.")), .funs = list(~if_else(is.na(.), FALSE, .))) %>% 
  dplyr::filter(!is.na(SM_ID)) %>% 
  dplyr::group_by(SM_CLUSTER) %>% 
  dplyr::mutate(index = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(geneId, SM_CLUSTER, index, starts_with("hasPeak"), starts_with("peakPval"),
                unname(polIICols$exp), unname(polIICols$is_expressed))

##################################################################################
## heatbox with KERS targets

# mergedData <- dplyr::filter(mergedData, SM_CLUSTER %in% c("01", "02", "03", "04", "05"))
plotTitle <- paste("SM genes | KERS targets |", polII_ids)

## prepare peak data
hasPeak <- pivot_longer(
  data = mergedData,
  cols = c(starts_with("hasPeak"), starts_with("peakPval")),
  names_to = c(".value", "sampleId"),
  names_sep = "\\.",
  values_drop_na = TRUE
)  %>% 
  dplyr::filter(hasPeak == TRUE) %>% 
  dplyr::select(geneId, SM_CLUSTER, index, sampleId, hasPeak, peakPval) %>% 
  dplyr::mutate(
    sampleId = forcats::fct_relevel(.f = sampleId, tfIds)
  )


tfPointPosition <- structure(
  .Data = seq(0, 1, length.out = length(tfInfo$sampleId)+2)[2:(1+length(tfInfo$sampleId))],
  names = tfInfo$sampleId
)

hasPeak$tfPosition <- tfPointPosition[hasPeak$sampleId]

tfColor <- purrr::map_chr(
  .x = exptInfoList[tfInfo$sampleId],
  .f = function(x) unlist(strsplit(x = x$color, split = ","))[2]
)

## prepare polII signal data
polIISignal <- as.data.table(mergedData) %>% 
  data.table::melt(
    id.vars = c("geneId", "SM_CLUSTER", "index"),
    measure.vars = list(polII = polIICols$exp, is_expressed = polIICols$is_expressed),
    variable.name = "sampleId"
  ) %>% 
  as_tibble() %>% 
  dplyr::mutate(polII = if_else(is_expressed == FALSE, 0, polII))

## prepare polII FC data
pivotSpec <- tibble::tibble(
  .name = c(polIICols$exp, polIICols$is_expressed),
  .value = c(rep("polII", length(polIICols$exp)), rep("is_expressed", length(polIICols$is_expressed))),
  sampleId = names(c(polIICols$exp, polIICols$is_expressed))
)

polIILfc <- dplyr::select(mergedData, !starts_with(c("hasPeak.", "peakPval"))) %>% 
  pivot_longer_spec(
    spec = pivotSpec
  )  %>% 
  dplyr::mutate(
    sampleId = forcats::fct_relevel(.f = sampleId, polII_ids)
  )

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
        legend.direction = "horizontal",
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = unit(rep(0.5, 4), "cm"))

pt <- ggplot() +
  geom_tile(data = polIISignal, mapping = aes(x = 0.5, y = sampleId, fill = polII)) +
  geom_point(data = hasPeak,
             mapping = aes(x = tfPosition, y = "tf", color = sampleId),
             size = 2) +
  scale_fill_gradientn(
    name = "log2(polII-signal)",
    colours = c("white", RColorBrewer::brewer.pal(n = 9, name = "RdPu")),
    values = scales::rescale(
      c(0, quantile(unlist(polIIData[polIICols$exp]),
                    c(0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99, 0.995, 0.999), names =F)
      )
    )
  ) +
  scale_color_manual(name = "KERS peak", values = tfColor) +
  scale_x_continuous(expand = expansion(add = c(0.0, 0.0))) +
  ggtitle(plotTitle) +
  facet_wrap(facets =  vars(geneId), ncol = 30, strip.position = "left", dir = "v") +
  ptTheme + theme(strip.text.y = element_blank())


pdf(file = paste(outPrefix, ".tile_plot.pdf", sep = ""), width = 15, height = 10)
print(pt)
dev.off()










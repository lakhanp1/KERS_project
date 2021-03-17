suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggpubr))

## Tile plot showing KERS binding and log2(mutant/WT) for each of the SM clusters gene

rm(list = ls())

##################################################################################
analysisName <- "binding_vs_del_polII_LFC_20h"
outDir <- here::here("analysis", "07_SM_analysis", "KERS_20h_SM")

outPrefix <- paste(outDir, "/", analysisName, sep = "")

tfIds <- c("An_kdmB_20h_HA_1", "An_ecoA_20h_HA_1", "An_rpdA_20h_HA_1", "An_sntB_20h_HA_1")
polIIComparison <- "polII_20h.kdmB_del_vs_untagged"


file_polIIRatioConf <- here::here("data", "reference_data", "polII_ratio.config.tab")

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

## polII ratio configuration
polIIRatioConf <- suppressMessages(readr::read_tsv(file = file_polIIRatioConf)) %>% 
  dplyr::filter(comparison %in% polIIComparison)

polIIDiffPairs <- purrr::transpose(polIIRatioConf)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))


##################################################################################

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfIds,
  dataPath = TF_dataPath
)

polII_info <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = unique(purrr::map(polIIDiffPairs, `[`, c("group2", "group1")) %>% unlist() %>% unname()),
  dataPath = polII_dataPath
)

exptInfo <- dplyr::bind_rows(tfInfo, polII_info)
exptInfoList <- purrr::transpose(exptInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

polII_ids <- exptInfo$sampleId[which(exptInfo$IP_tag == "polII")]
tfIds <- exptInfo$sampleId[which(exptInfo$IP_tag %in% c("HA", "MYC", "TAP") & exptInfo$TF != "untagged")]


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

i <- 1
## add fold change columns
for (i in names(polIIDiffPairs)) {
  polIIData <- get_fold_change(
    df = polIIData,
    nmt = polIIDiffPairs[[i]]$group1,
    dmt = polIIDiffPairs[[i]]$group2,
    newCol = polIIDiffPairs[[i]]$comparison,
    lfcLimit = 0.585,
    isExpressedCols = polIICols$is_expressed
  )
}


mergedData <- dplyr::left_join(x = polIIData, y = tfBindingMat, by = "geneId") %>% 
  dplyr::mutate_at(.vars = vars(starts_with("hasPeak.")), .funs = list(~if_else(is.na(.), FALSE, .))) %>% 
  dplyr::filter(!is.na(SM_ID)) %>% 
  dplyr::mutate(
    SM_ID = gsub(pattern = "cluster_", replacement = "", x = SM_ID, fixed = TRUE)
  ) %>% 
  dplyr::group_by(SM_CLUSTER) %>% 
  dplyr::mutate(index = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(geneId, SM_CLUSTER, index, starts_with("hasPeak"), starts_with("peakPval"),
                unname(polIICols$exp), unname(polIICols$is_expressed),
                unname(purrr::map_chr(polIIDiffPairs, "comparison")))


##################################################################################
## polII FC box with KERS targets

# mergedData <- dplyr::filter(mergedData, SM_CLUSTER %in% c("01", "02", "03", "04", "05"))

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


## prepare polII FC data
polIILfc <- pivot_longer(
  data = mergedData,
  cols = unname(purrr::map_chr(polIIDiffPairs, "comparison")),
  names_to = "pair",
  values_to = "lfc",
  values_drop_na = TRUE
)  %>% 
  dplyr::select(geneId, SM_CLUSTER, index, pair, lfc) %>% 
  dplyr::mutate(
    pair = forcats::fct_relevel(.f = pair, unname(purrr::map_chr(polIIDiffPairs, "comparison")))
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

plotTitle <- paste("SM genes | KERS targets |", polIIComparison)

pt <- ggplot() +
  geom_tile(data = polIILfc, mapping = aes(x = 0.5, y = pair, fill = lfc)) +
  geom_point(data = hasPeak,
             mapping = aes(x = tfPosition, y = "tf", color = sampleId),
             size = 2) +
  scale_fill_gradient2(
    name = "log2(polII-fold-change)",
    low = "#B35806", mid = "#F7F7F7", high = "#542788", midpoint = 0
  ) +
  scale_color_manual(name = "KERS peak", values = tfColor) +
  scale_x_continuous(expand = expansion(add = c(0.0, 0.0))) +
  scale_y_discrete(limits = c(levels(polIILfc$pair), "tf")) +
  ggtitle(plotTitle) +
  facet_wrap(facets =  vars(geneId), ncol = 35, strip.position = "left", dir = "v") +
  ptTheme + theme(strip.text.y = element_blank())


pdf(file = paste(outPrefix, ".tile_plot.pdf", sep = ""), width = 15, height = 10)
print(pt)
dev.off()










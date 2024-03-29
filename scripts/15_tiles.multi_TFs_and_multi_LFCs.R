suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggpubr))

## Tile plot with facets of gene groups showing
## 1) binding of multiple TFs using geom_point
## 2) LFC of different comparisons using geom_tiles

rm(list = ls())

##################################################################################
analysisName <- "kdmB"
outDir <- here::here("analysis", "07_SM_analysis")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

tfIds <- c("An_kdmB_20h_HA_1", "An_kdmB_48h_HA_1")

## polII signal fold change pairs
polIIDiffPairs <- c("polII_untagged.48h_vs_20h", "polII_kdmB_del.48h_vs_20h",
                    "polII_20h.kdmB_del_vs_untagged", "polII_48h.kdmB_del_vs_untagged")

## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")
file_polIIRatioConf <- here::here("data", "reference_data", "polII_ratio.config.tab")

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")

orgDb <- org.Anidulans.FGSCA4.eg.db


##################################################################################
if(!dir.exists(outDir)) dir.create(path = outDir, recursive = TRUE)

## genes to read
geneSet <- suppressMessages(
  readr::read_tsv(file = file_genes, col_names = c("chr", "start", "end", "geneId", "score", "strand"))
) %>% 
  dplyr::select(geneId, chr, start, end, strand)

geneDesc <- suppressMessages(
  AnnotationDbi::select(x = orgDb, keys = geneSet$geneId,
                        columns = c("DESCRIPTION"), keytype = "GID")
)

smInfo <- suppressMessages(
  AnnotationDbi::select(x = orgDb, keys = keys(orgDb, keytype = "SM_CLUSTER"),
                        columns = c("GID", "SM_ID"), keytype = "SM_CLUSTER")) 

geneInfo <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID")) %>% 
  dplyr::left_join(y = smInfo, by = c("geneId" = "GID"))

glimpse(geneInfo)

## polII ratio configuration
polIIRatioConf <- suppressMessages(readr::read_tsv(file = file_polIIRatioConf)) %>% 
  dplyr::filter(comparison %in% polIIDiffPairs)

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
  dplyr::group_by(SM_CLUSTER) %>% 
  dplyr::arrange(chr, start, .by_group = TRUE) %>% 
  dplyr::mutate(index = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(
    geneId, SM_ID, SM_CLUSTER, index, starts_with("hasPeak"), starts_with("peakPval"),
    unname(polIICols$exp), unname(polIICols$is_expressed),
    unname(purrr::map_chr(polIIDiffPairs, "comparison"))
  )


glimpse(mergedData)


##################################################################################

## SM cluster polII LFC signal plot

pltTitle <- paste(c("SM cluster polII fold change:", analysisName), collapse = " ")

## prepare peak data
hasPeak <- pivot_longer(
  data = mergedData,
  cols = c(starts_with("hasPeak"), starts_with("peakPval")),
  names_to = c(".value", "sampleId"),
  names_sep = "\\.",
  values_drop_na = TRUE
)  %>% 
  dplyr::filter(hasPeak == TRUE) %>%
  dplyr::select(geneId, SM_ID, SM_CLUSTER, index, sampleId, hasPeak, peakPval) %>% 
  dplyr::mutate(
    sampleId = forcats::fct_relevel(.f = sampleId, tfIds)
  )

tfColor <- purrr::map_chr(
  .x = exptInfoList[tfIds],
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

yOrder <- unname(c(purrr::map_chr(polIIDiffPairs, "comparison"), tfIds))

##################################################################################

ptTheme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        strip.background = element_rect(fill="white", size = 0.1),
        strip.text.x = element_text(size = 12, hjust = 0, margin = margin(1,1,1,1)),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = unit(rep(0.2, 4), "cm"))


pt2 <- ggplot() +
  geom_point(
    data = hasPeak,
    mapping = aes(x = index , y = sampleId, color = sampleId),
    size = 2, shape = 16
  ) +
  geom_tile(
    data = polIILfc,
    mapping = aes(x = index, y = pair, fill = lfc),
    color = "black", size = 0.2, height = 1) +
  scale_fill_gradient2(
    name = "log2(polII-fold-change)",
    low = "#B35806", mid = "#F7F7F7", high = "#542788", midpoint = 0
  ) +
  scale_colour_discrete(
    # type = tfColor,
    breaks = tfIds, name = ""
  ) +
  scale_x_continuous(expand = expansion(add = c(0.0, 0.0))) +
  scale_y_discrete(
    limits = yOrder,
    expand = expansion(add = c(0.0, 0.5))
  ) +
  facet_wrap(
    facets = . ~ SM_CLUSTER, scales = "free_y",
    ncol = 6, dir = "v"
  ) +
  ggtitle(
    label = pltTitle,
    subtitle = paste(
      "Each tile is a gene representing", length(yOrder), "values (bottom to top):",
      paste(yOrder, collapse = ", ")
    )
  ) +
  ptTheme +
  theme(
    legend.position = "bottom",
    # legend.justification = c("right", "bottom"),
    legend.key.width = unit(2, "cm")
  )


pdf(file = paste(outPrefix, ".binding_and_polII_LFC.tile_plot.pdf", sep = ""), width = 18, height = 10)
pt2
dev.off()


##################################################################################






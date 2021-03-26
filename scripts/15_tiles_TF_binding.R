suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggpubr))

## 1) plots SM cluster wise genes binding signal plot as geom_tiles

rm(list = ls())

##################################################################################
analysisName <- "kdmB"
outDir <- here::here("analysis", "07_SM_analysis", "kdmB")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

tfIds <- c("An_kdmB_20h_HA_1", "An_kdmB_48h_HA_1")


## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

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

##################################################################################

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfIds,
  dataPath = TF_dataPath
)


exptInfo <- dplyr::bind_rows(tfInfo)
exptInfoList <- purrr::transpose(exptInfo)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

tfIds <- exptInfo$sampleId[which(exptInfo$IP_tag %in% c("HA", "MYC", "TAP") & exptInfo$TF != "untagged")]

tfCols <- sapply(
  c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
    "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
    "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
  simplify = F, USE.NAMES = T)

##################################################################################
## prepare data for plotting
tfBindingMat <- peak_target_matrix(sampleInfo = tfInfo, position = "best")

mergedData <- dplyr::left_join(x = geneInfo, y = tfBindingMat, by = "geneId") %>% 
  dplyr::mutate_at(.vars = vars(starts_with("hasPeak.")), .funs = list(~if_else(is.na(.), FALSE, .))) %>% 
  dplyr::filter(!is.na(SM_ID)) %>% 
  dplyr::group_by(SM_CLUSTER) %>% 
  dplyr::arrange(chr, start, .by_group = TRUE) %>% 
  dplyr::mutate(index = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(
    geneId, SM_ID, SM_CLUSTER, index, starts_with("hasPeak"), starts_with("peakPval")
  )


glimpse(mergedData)

##################################################################################

ptTheme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        strip.background = element_rect(fill="white", size = 0.1),
        strip.text.x = element_text(size = 10, hjust = 0, margin = margin(1,1,1,1)),
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = unit(rep(0.5, 4), "cm"))

## SM cluster binding plot

## prepare peak data
hasPeak <- pivot_longer(
  data = mergedData,
  cols = c(starts_with("hasPeak"), starts_with("peakPval")),
  names_to = c(".value", "sampleId"),
  names_sep = "\\.",
  values_drop_na = TRUE
)  %>% 
  # dplyr::filter(hasPeak == TRUE) %>% 
  dplyr::select(geneId, SM_ID, SM_CLUSTER, index, sampleId, hasPeak, peakPval) %>% 
  dplyr::mutate(
    sampleId = forcats::fct_relevel(.f = sampleId, tfIds)
  )


pltTitle <- paste("SM cluster binding comparison:", analysisName)

pt1 <- ggplot(data = hasPeak) +
  geom_tile(
    mapping = aes(x = index, y = sampleId, fill = hasPeak, color = sampleId),
    size = 0.75, height = 0.9
  ) +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "white"),
    guide = FALSE
  ) +
  scale_color_discrete(
    breaks = tfIds,
    # values = structure(c("#B35806", "#542788"), names = unname(tfCols$hasPeak))
    limits = tfIds,
    name = ""
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(
    facets = . ~ SM_CLUSTER, scales = "free_y",
    ncol = 5, dir = "v"
  ) +
  ggtitle(label = pltTitle) +
  ptTheme


pdf(file = paste(outPrefix, ".SM_clusters_binding.tile_plot.pdf", sep = ""), width = 15, height = 8)
pt1
dev.off()



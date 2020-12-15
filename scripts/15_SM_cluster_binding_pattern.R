library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(esquisse)
library(summarytools)
library(here)
library(ggpubr)

## 1) plots SM cluster wise genes binding signal plot as geom_tiles
## 2) plots cluster wise binding signal and polII signal fold change plot using geom_tile
## 3) plots the binding signal and polII fold change values in two groups: bound and unbound

rm(list = ls())

##################################################################################
analysisName <- "kdmB_del"
outPrefix <- here::here("kdmB_analysis", "SM_analysis", analysisName)

tfIds <- c("An_kdmB_20h_HA_1", "An_kdmB_48h_HA_1")
polII1 <- "An_untagged_20h_polII_1"
polII2 <- "An_untagged_48h_polII_1"

otherPolIIs <- c("An_kdmB_del_20h_polII_1", "An_kdmB_del_48h_polII_1")

## polII signal fold change pairs
polIIDiffPairs <- list(
  p1 = list(
    name = "48h_vs_20h_untagged_polII",
    title = "polII log2(untagged_48h \n vs untagged_20h)",
    samples = c(polII1, polII2)
  ),
  p2 = list(
    name = "48h_vs_20h_kdmB_del_polII",
    title = "polII log2(kdmB_del_48h \n vs kdmB_del_20h)",
    samples = c(otherPolIIs[1], otherPolIIs[2])
  ),
  p3 = list(
    name = "kdmB_del_vs_untagged_20h_polII",
    title = "polII log2(kdmB_del \n vs kdmB_untagged 20h)",
    samples = c("An_untagged_20h_polII_1", "An_kdmB_del_20h_polII_1")
  ),
  p4 = list(
    name = "kdmB_del_vs_untagged_48h_polII",
    title = "polII log2(kdmB_del \n vs kdmB_untagged 48h)",
    samples = c("An_untagged_48h_polII_1", "An_kdmB_del_48h_polII_1")
  )
)

polIIPairId <- "p1"

## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")


orgDb <- org.Anidulans.FGSCA4.eg.db

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
                                     samples = c(polII1, polII2, otherPolIIs),
                                     dataPath = polII_dataPath,
                                     matrixSource = "normalizedmatrix")

exptData <- dplyr::bind_rows(tfInfo, polII_info)

polII_ids <- exptData$sampleId[which(exptData$IP_tag == "polII")]
tfIds <- exptData$sampleId[which(exptData$IP_tag %in% c("HA", "MYC", "TAP") & exptData$TF != "untagged")]


polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)

tfCols <- sapply(
  c("peakDist", "featureCovFrac", "hasPeak", "peakCoverage", "peakPosition", "peakId", "peakType",
    "peakPval", "peakEnrichment", "preference", "peakCategory"),
  FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
  simplify = F, USE.NAMES = T)


##################################################################################

chipData <- get_polII_expressions(genesDf = geneSet, exptInfo = polII_info)

## add fold change columns
for (i in names(polIIDiffPairs)) {
  chipData <- get_fold_change(df = chipData,
                              nmt = polIIDiffPairs[[i]]$samples[2],
                              dmt = polIIDiffPairs[[i]]$samples[1],
                              newCol = polIIDiffPairs[[i]]$name,
                              lfcLimit = 0.58,
                              isExpressedCols = polIICols$is_expressed)
}

## add TF binding information
chipData <- get_TF_binding_data(genesDf = chipData, exptInfo = tfInfo, allColumns = FALSE) %>% 
  dplyr::filter(! is.na(SM_CLUSTER))

# dplyr::filter(chipData, ! is.na(SM_CLUSTER)) %>%
#   dplyr::group_by_at(.vars = vars(unname(polIICols$is_expressed[polIIDiffPairs$p1$samples]))) %>%
#   dplyr::summarise(n = n())
# 
# dplyr::group_by_at(chipData, .vars = vars(starts_with("hasPeak."))) %>%
#   dplyr::summarise(n = n())

chipData <- dplyr::mutate(
  chipData,
  SM_CLUSTER = gsub(pattern = "SM_cluster_", replacement = "", x = SM_CLUSTER, fixed = TRUE)) %>% 
  dplyr::group_by(SM_CLUSTER) %>% 
  dplyr::arrange(start, .by_group = TRUE) %>% 
  dplyr::mutate(index = 1:n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(gene, SM_CLUSTER, starts_with("hasPeak."), index, !! unname(purrr::map_chr(polIIDiffPairs, "name"))) %>% 
  as.data.frame()

view(dfSummary(chipData))



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
        # legend.direction = "vertical",
        legend.title = element_text(size = 13, face = "bold"),
        plot.margin = unit(rep(0.5, 4), "cm"))

##################################################################################


## SM cluster binding plot
pltDf1 <- tidyr::gather(chipData, key = "sample", value = "hasPeak",
                        -gene, -SM_CLUSTER, -index, - ends_with("_polII"))

pltDf1$sample <- factor(pltDf1$sample, levels = c(unname(tfCols$hasPeak)))

pltTitle <- paste("SM cluster binding comparison:", analysisName)

pt1 <- ggplot(data = pltDf1) +
  geom_tile(mapping = aes(x = index, y = sample, fill = hasPeak, color = sample), size = 0.75, height = 0.9) +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "white"),
    guide = FALSE
  ) +
  scale_color_discrete(
    breaks = unname(tfCols$hasPeak),
    # values = structure(c("#B35806", "#542788"), names = unname(tfCols$hasPeak))
    name = ""
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(facets = SM_CLUSTER ~ ., scales = "free_y", ncol = 2, strip.position = "left", dir = "v") +
  ggtitle(pltTitle) +
  ptTheme


pdf(file = paste(outPrefix, "_cluster_binding_cmp.pdf", sep = ""), width = 10, height = 10)
pt1
dev.off()

##################################################################################
## SM cluster polII signal plot

pltTitle <- paste(c("SM cluster polII fold change:", analysisName), collapse = " ")


pt2 <- ggplot() +
  geom_point(
    data = tidyr::gather(chipData, key = "sample", value = "hasPeak",
                         -gene, -SM_CLUSTER, -index, - ends_with("_polII")) %>%
      dplyr::filter(hasPeak != 0) %>%
      dplyr::mutate(hasPeak = as.character(hasPeak)),
    mapping = aes(x = index , y = sample, color = sample),
    size = 2, shape = 16
  ) +
  geom_tile(
    data = tidyr::gather(chipData, key = "sample", value = "polII",
                         -gene, -SM_CLUSTER, -index, - starts_with("hasPeak.")),
    mapping = aes(x = index, y = sample, fill = polII),
    color = "black", size = 0.2, height = 1) +
  scale_fill_gradient2(
    name = paste("log2(", "polII fold change", ")", sep = ""),
    low = "#B35806", mid = "#F7F7F7", high = "#542788", midpoint = 0
  ) +
  scale_colour_manual(
    name = "",
    values = structure(c("red", "blue"), names = tfCols$hasPeak),
    labels = names(tfCols$hasPeak),
    breaks = unname(tfCols$hasPeak)) +
  scale_x_continuous(expand = expand_scale(add = c(0.0, 0.0))) +
  scale_y_discrete(
    limits = unname(c(rev(purrr::map_chr(polIIDiffPairs, "name")), tfCols$hasPeak[2], tfCols$hasPeak[1])),
    expand = expand_scale(add = c(0.0, 0.5))
  ) +
  facet_wrap(facets = SM_CLUSTER ~ ., scales = "free_y", ncol = 4, strip.position = "left", dir = "v",
             labeller = labeller(vs = label_both, am = label_value)) +
  ggtitle(pltTitle) +
  ptTheme



# png(filename = paste(outPrefix, "_SM_cluster_polII_diff.png"), width = 4000, height = 6000, res = 350)
pdf(file = paste(outPrefix, "_SM_cluster_polII_diff.pdf", sep = ""), width = 18, height = 10)
pt2
dev.off()


##################################################################################
## plot tiles in two groups, showing kdmB peak and not showing peak
pt3Df <- dplyr::mutate(.data = chipData,
                       hasPeak = !! as.name(tfCols$hasPeak[1]) | !! as.name(tfCols$hasPeak[2]) ) %>% 
  dplyr::arrange_at(.vars = unname(tfCols$hasPeak), .funs = funs(desc(.)))

## convert gene column to factor 
pt3Df$gene <- factor(pt3Df$gene, levels = pt3Df$gene)

# df = pt3Df

## function to generate tile plot for a df
SM_tile_plot <- function(df, nrow){
  pt = ggplot() +
    # geom_point(
    #   data = tidyr::gather(df, key = "tf", value = "tfPeak",
    #                        -gene, -SM_CLUSTER, -index, -hasPeak, -ends_with("_polII")) %>% 
    #     dplyr::mutate(tfPeak = as.logical(tfPeak)) %>% 
    #     dplyr::mutate(peak = if_else(tfPeak == TRUE, tf, "FALSE")),
    #   mapping = aes(x = 1, y = tf, color = peak),
    #   size = 2, shape = 16
    # ) +
    geom_tile(
      data = tidyr::gather(df, key = "polII", value = "lfc",
                           -gene, -SM_CLUSTER, -index, -hasPeak, -starts_with("hasPeak.")),
      mapping = aes(x = 1, y = polII, fill = lfc),
      color = "black", size = 0.2, height = 1
    ) +
    scale_fill_gradient2(
      name = paste("log2(", "polII fold change", ")", sep = ""),
      low = "#B35806", mid = "#F7F7F7", high = "#542788", midpoint = 0
    ) +
    # scale_colour_manual(
    #   name = "",
    #   breaks = unname(tfCols$hasPeak),
    #   values = structure(c("red", "blue"), names = tfCols$hasPeak),
    #   labels = names(tfCols$hasPeak)) +
    scale_y_discrete(
      limits = unname(c(rev(purrr::map_chr(polIIDiffPairs, "name")))),
      # limits = unname(c(rev(purrr::map_chr(polIIDiffPairs, "name")), tfCols$hasPeak[2], tfCols$hasPeak[1])),
      expand = expand_scale(add = c(0.0, 0.5))
    ) +
    scale_x_continuous(expand = expand_scale(add = c(0.0, 0.0))) +
    facet_wrap(facets = gene ~ ., scales = "free", nrow = nrow, strip.position = "left", dir = "v") +
    theme_bw() +
    ptTheme + theme(strip.text.y = element_blank())
  
  return(pt)
}


statTab <- dplyr::group_by_at(pt3Df, .vars = vars(unname(tfCols$hasPeak))) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange_at(.vars = vars(unname(tfCols$hasPeak)), desc) %>% 
  dplyr::mutate_if(.predicate = is.logical, .funs = as.character)

pt3Table <- ggtexttable(statTab, rows = NULL, theme = ttheme("lCyanWhite"))

pt3.1 <- SM_tile_plot(
  df = dplyr::filter(.data = pt3Df, !! as.name(tfCols$hasPeak[1]) == TRUE, !! as.name(tfCols$hasPeak[2]) == TRUE),
  nrow = 16
)



pt3.2 <- SM_tile_plot(
  df = dplyr::filter(.data = pt3Df, !! as.name(tfCols$hasPeak[1]) == TRUE, !! as.name(tfCols$hasPeak[2]) == FALSE),
  nrow = 16
)

pt3.3 <- SM_tile_plot(
  df = dplyr::filter(.data = pt3Df, !! as.name(tfCols$hasPeak[1]) == FALSE, !! as.name(tfCols$hasPeak[2]) == TRUE),
  nrow = 16
)

pt3.4 <- SM_tile_plot(
  df = dplyr::filter(.data = pt3Df, !! as.name(tfCols$hasPeak[1]) == FALSE, !! as.name(tfCols$hasPeak[2]) == FALSE),
  nrow = 16
)

# c(72, 32, 60, 411)
# c(5, 3, 4, 23)
pt3 <- ggpubr::ggarrange(
  ggarrange(
    pt3.1, pt3.2, pt3.3, pt3.4, nrow = 1, ncol = 4, common.legend = TRUE, legend = "bottom", widths = c(5, 3, 4, 23)
  ),
  pt3Table,
  ncol = 1, nrow = 2, heights = c(8, 2)
)


pdf(file = paste(outPrefix, "_SM_cluster_hasPeak_pairs.pdf", sep = ""), width = 12, height = 12)
pt3
dev.off()







suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(summarytools))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))

## generate profile plots for KERS complex
## group the genes into different categories based on kdmB complex member binding status
## plot the profile plots for groups
## plot average signal for the groups

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")
##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
analysisName <- "KERS_complex_20h"
outDir <- here::here("analysis", "03_KERS_complex", analysisName)
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_plotSamples <- paste(outDir, "/", "samples.txt", sep = "")

matrixType <- "normalizedMatrix_5kb"
up <- 5000
body <- 2000
down <- 1000
binSize <- 10
matrixDim = c(c(up, body, down)/binSize, binSize)

showExpressionHeatmap = TRUE

## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/annotation_resources/geneid2go.ANidulans.topGO.map"

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db

anLables <- list()

##################################################################################

sampleList <- suppressMessages(readr::read_tsv(file = file_plotSamples, comment = "#"))

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

## read the experiment sample details and select only those which are to be plotted
tempSInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = sampleList$sampleId,
  dataPath = TF_dataPath,
  profileMatrixSuffix = matrixType
)


polII_ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "polII")]
tfIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag %in% c("HA", "MYC", "TAP") & tempSInfo$TF != "untagged")]
inputIds <- tempSInfo$sampleId[which(! tempSInfo$IP_tag %in% c("polII", "HIST") & tempSInfo$TF == "untagged")]
histIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST")]


## read the experiment sample details and select only those which are to be plotted
tfData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = tfIds,
  dataPath = TF_dataPath, profileMatrixSuffix = matrixType
)


inputData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = inputIds,
  dataPath = TF_dataPath, profileMatrixSuffix = matrixType
)

polIIData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = polII_ids,
  dataPath = polII_dataPath, profileMatrixSuffix = matrixType
)

histData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = histIds,
  dataPath = hist_dataPath, profileMatrixSuffix = matrixType
)

exptData <- dplyr::bind_rows(tfData, inputData, histData, polIIData)

exptDataList <- purrr::transpose(exptData) %>% 
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

peakTargetMat <- peak_target_matrix(sampleInfo = tfData, position = "best")

peakTargetMat <- dplyr::mutate_at(
  .tbl = peakTargetMat,
  .vars = vars(starts_with("hasPeak.")),
  .funs = ~ factor(., levels = c(TRUE, FALSE))
)

expressionData <- get_polII_expressions(exptInfo = exptData,
                                        genesDf = geneInfo)

expressionData <- dplyr::left_join(x = expressionData, y = peakTargetMat, by = "geneId")

hasPeakDf <- expressionData %>% 
  dplyr::filter_at(.vars = vars(starts_with("hasPeak")), .vars_predicate = any_vars(. == TRUE)) %>% 
  dplyr::mutate(group = group_indices(., !!! lapply(unname(tfCols$hasPeak), as.name))) %>% 
  dplyr::mutate(group = sprintf(fmt = "%02d", group))

dplyr::group_by_at(hasPeakDf, .vars = vars(starts_with("hasPeak."))) %>% 
  dplyr::summarise(n = n_distinct(geneId))


readr::write_tsv(x = hasPeakDf, file = paste(outPrefix, ".peaks_data.tab", sep = ""))

dplyr::filter_at(
  hasPeakDf, .vars = vars(starts_with("hasPeak.")), .vars_predicate = all_vars(. == TRUE)
) %>% 
  dplyr::select(geneId) %>% 
  readr::write_tsv(
    file = paste(outPrefix, ".common_genes.tab")
  )

## group label data mean profile facet plot
groupLabelDf <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(starts_with("hasPeak."), group)) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(groupLabels = paste(group, ": ", n, " genes", sep = ""))

groupLabels <- structure(groupLabelDf$groupLabels, names = groupLabelDf$group)

##################################################################################
## binding stats
bindingMat = combinatorial_binding_matrix(sampleInfo = tfData, summitRegion = 100)

readr::write_tsv(x = bindingMat, file = paste(outPrefix, ".binding_matrix.tab", sep = ""))

##################################################################################
## topGO enrichment
goEnrich <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(starts_with("hasPeak."), group)) %>%
  do(
    topGO_enrichment(
      genes = .$geneId, orgdb = orgDb, inKeytype = "GID",
      type = "BP", goNodeSize = 5, genenameKeytype = "GENE_NAME"
    )
  )


readr::write_tsv(x = goEnrich, file = paste(outPrefix, ".peakGroups.topGO.tab", sep = ""))


## pathway enrichment
keggEnr <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(starts_with("hasPeak."), group)) %>%
  do(
    tibble::as_tibble(
      fgsea_kegg_overrepresentation(
        genes = .$geneId, orgdb = orgDb, keggOrg = "ani", inKeytype = "GID",
        keggKeytype = "KEGG_ID", genenameKeytype = "GENE_NAME", pvalueCutoff = 0.05
      )
    )
  )


readr::write_tsv(x = keggEnr, file = paste(outPrefix, ".peakGroups.KEGG_enrichment.tab", sep = ""))


##################################################################################

## profile matrix of the genes which show binding
matList <- import_profiles(
  exptInfo = exptData, geneList = unique(geneInfo$geneId), source = matrixType,
  up = matrixDim[1], target = matrixDim[2], down = matrixDim[3]
)

## get average signal over all factors to select color
tfMeanProfile <- getSignalsFromList(lt = matList[tfIds])
quantile(tfMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)

# tfMeanColor <- colorRamp2(quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T), c("white", "red"))
tfColorList <- sapply(
  X = c(tfIds, inputIds),
  FUN = function(x){
    return(
      colorRamp2(
        breaks = quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T),
        colors = unlist(strsplit(x = exptDataList[[x]]$color, split = ","))
      )
    )
  }
)

colorList <- tfColorList

##################################################################################
## polII signal matrix
polIIMat <- data.matrix(log2(expressionData[, polII_ids] + 1))
rownames(polIIMat) <- expressionData$geneId

quantile(polIIMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.997, 0.999, 0.9999, 1), na.rm = T)


polII_color <- colorRamp2(breaks = c(0, quantile(polIIMat, c(0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99, 0.995, 0.999))),
                          colors = c("white", RColorBrewer::brewer.pal(n = 9, name = "RdPu")))


##################################################################################
ylimList <- sapply(exptData$sampleId, function(x){return(c(0,25))}, simplify = FALSE)

profiles_peaks <- multi_profile_plots(
  exptInfo = exptData[exptData$sampleId %in% c(tfIds, inputIds), ],
  genesToPlot = hasPeakDf$geneId,
  clusters = NULL,
  profileColors = colorList,
  matSource = matrixType,
  matBins = matrixDim,
  ylimFraction = ylimList,
  column_title_gp = gpar(fontsize = 14)
)

anGl_peaks <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = hasPeakDf$geneId,
  axis_param = list(at = c(0, 2000, 4000), labels = c("0kb", "2kb", ">4kb"))
)

peaks_htlist <- profiles_peaks$heatmapList + anGl_peaks$an

pdfWd <- 3 + (length(profiles_peaks$heatmapList@ht_list) * 2) +
  (length(polII_ids) * 0.25 * showExpressionHeatmap) + 3

title_peak= "Transcription factor binding profile: kdmB complex TF bound genes"

# draw Heatmap and add the annotation name decoration
pdf(file = paste(outPrefix, ".profilePlot_unclustered.pdf", sep = ""), width = pdfWd, height = 12)

draw(peaks_htlist,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_peak,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "right",
     gap = unit(7, "mm"),
     padding = unit(rep(0.3, times = 4), "cm")
)

dev.off()



##################################################################################

## profile matrix with selected genes only
newClusters <- dplyr::select(hasPeakDf, geneId, group) %>% 
  dplyr::rename(cluster = group)

newClusters$cluster <- factor(x = newClusters$cluster, levels = sort(unique(newClusters$cluster)))

profiles_peakGroups <- multi_profile_plots(
  exptInfo = exptData[exptData$sampleId %in% c(tfIds, inputIds), ],
  genesToPlot = hasPeakDf$geneId,
  clusters = newClusters,
  profileColors = colorList,
  matSource = matrixType,
  matBins = matrixDim,
  column_title_gp = gpar(fontsize = 12)
)

## polII signal heatmap
polIIMat_peakGroups <- polIIMat[hasPeakDf$geneId, , drop=FALSE]

polIIht_peakExp <- signal_heatmap(
  log2_matrix = polIIMat_peakGroups,
  htName = "polII_signal",
  col_title = polII_ids,
  legend_title = "log2(polII_singal)",
  color = polII_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE)


## heatmap of binary assignment of samples to different group
htMat <- dplyr::select(hasPeakDf, geneId, starts_with("hasPeak")) %>% 
  dplyr::mutate_if(.predicate = is.logical, .funs = as.character) %>% 
  tibble::column_to_rownames("geneId")

## column name as annotation for Heatmap
colNameAnn <- HeatmapAnnotation(
  colName = anno_text(
    x = unname(tfCols$hasPeak),
    rot = 90, just = "left",
    location = unit(1, "mm"),
    gp = gpar(fontsize = 10)),
  annotation_height = unit.c(max_text_width(unname(tfCols$hasPeak)))
)

grp_ht <- Heatmap(
  as.matrix(htMat),
  col = c("TRUE" = "black", "FALSE" = "white"),
  heatmap_legend_param = list(title = "Peak detected"),
  # column_names_side = "top",
  show_column_names = FALSE,
  top_annotation = colNameAnn,
  cluster_columns = FALSE, cluster_rows = FALSE,
  width = unit(3, "cm"),
  show_row_names = FALSE
)

## gene length annotation
anGl_peakGroups <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = hasPeakDf$geneId,
  axis_param = list(at = c(0, 2000, 4000), labels = c("0kb", "2kb", ">4kb"))
)


peakGroups_htlist <- profiles_peakGroups$heatmapList + polIIht_peakExp + grp_ht + anGl_peakGroups$an

pdfWd <- 3 + 
  (length(profiles_peakGroups$heatmapList@ht_list) * 2) +
  (length(polII_ids) * 0.5 * showExpressionHeatmap) + 3

title_peak= "Transcription factor binding profile: kdmB complex TF bound genes"

# draw Heatmap and add the annotation name decoration
pdf(file = paste(outPrefix, ".profilePlot_groups.pdf", sep = ""), width = pdfWd, height = 12)

draw(peakGroups_htlist,
     main_heatmap = exptData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_peak,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     heatmap_legend_side = "right",
     gap = unit(5, "mm"),
     padding = unit(rep(0.5, times = 4), "cm")
)


dev.off()


##################################################################################
## average line plots for peak groups
testDf <- dplyr::filter(hasPeakDf, group == "07")

lineColors = purrr::map_chr(.x = exptDataList[c(tfIds, inputIds)],
                            .f = function(x) unlist(strsplit(x = x$color, split = ","))[2])

lineShape = structure(.Data = c(1, 1, 1, 1, 6),
                      names = exptData$sampleId[exptData$sampleId %in% c(tfIds, inputIds)])


draw_avg_profile_plot(
  exptInfo = exptData[exptData$sampleId %in% c(tfIds, inputIds), ],
  profileMats = matList,
  genes = testDf$geneId,
  lineColors = lineColors,
  lineShape = lineShape
)


ap <- geneset_average_profile(
  exptInfo = exptData[exptData$sampleId %in% c(tfIds, inputIds), ],
  profileMats = matList,
  genes = testDf$geneId,
  cluster = "group_1"
)


## average profile data for each group
groupMeanProfiles <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(starts_with("hasPeak."), group)) %>%
  dplyr::do(
    geneset_average_profile(exptInfo = exptData[exptData$sampleId %in% c(tfIds, inputIds), ],
                            profileMats = matList,
                            genes = .$geneId,
                            cluster = unique(.$group))
  ) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()

## decide the axis labels
axisBrk <- NULL
axisLab <- NULL
targetEnd <- NULL
if(attributes(matList[[1]])$target_is_single_point){
  axisBrk <- c(
    attributes(matList[[1]])$upstream_index[1],
    attributes(matList[[1]])$target_index[1],
    tail(attributes(matList[[1]])$downstream_index, 1)
  )
  
  axisLab <- c(
    -attributes(matList[[1]])$extend[1],
    attributes(matList[[1]])$target_name,
    attributes(matList[[1]])$extend[2]
  )
  
  targetEnd <- axisBrk[2] + 1
  
} else if(! attributes(matList[[1]])$target_is_single_point){
  axisBrk <- c(
    attributes(matList[[1]])$upstream_index[1],
    attributes(matList[[1]])$target_index[1],
    tail(attributes(matList[[1]])$target_index, 1),
    tail(attributes(matList[[1]])$downstream_index, 1)
  )
  
  axisLab <- c(
    -attributes(matList[[1]])$extend[1],
    "START", "END",
    attributes(matList[[1]])$extend[2]
  )
  
  targetEnd <- axisBrk[3]
}

## plot the groups using facets
p = ggplot(data = groupMeanProfiles) +
  geom_line(mapping = aes(x = bin, y = mean, group = sample, color = sample, linetype = sample),
            size = 0.8, alpha = 0.8) +
  geom_hline(yintercept = 0, size = 2, color = "grey70") +
  geom_segment(mapping = aes(x = axisBrk[2], y = 0, xend = targetEnd, yend = 0),
               size = 10, lineend = "butt", color = "grey50") +
  scale_color_manual(values = lineColors) +
  scale_linetype_manual(values = lineShape) +
  scale_x_continuous(breaks = axisBrk, labels = axisLab) +
  ggtitle(paste("Average profile for different binding patterns of kdmB complex members: ", analysisName)) +
  ylab("Read coverage") +
  facet_wrap(group ~ ., ncol = 5, scales = "free_y", labeller = labeller(group = groupLabels)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 15),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, 'cm'),
        # legend.justification = c(1.1, 1.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        plot.margin = unit(rep(0.5, 4), "cm")
  ) +
  guides(
    color = guide_legend(
      label.position = "bottom")
  )


pdf(file = paste(outPrefix, ".mean_profiles.pdf", sep = ""), width = 18, height = 10)
p
dev.off()


##################################################################################
## upset plot
cm = make_comb_mat(
  dplyr::select(hasPeakDf, geneId, tfCols$hasPeak) %>% 
    dplyr::filter_at(.vars = vars(!geneId), .vars_predicate = any_vars(. == TRUE))
)

ht_upset <- UpSet(
  cm,
  pt_size = unit(4, "mm"), lwd = 2,
  set_order = tfData$sampleId,
  row_labels = tfData$TF,
  row_names_gp = gpar(fontsize = 14, fontface = "italic"),
  column_title = "KERS bound genes statistics",
  column_title_gp = gpar(fontsize = 20, fontface = "bold"),
  # bg_col = stringr::str_replace(string = tfData$color, pattern = "white,", replacement = ""),
  top_annotation = HeatmapAnnotation(
    # "combSize" = anno_text(
    #   x = paste("(", comb_size(cm), ")", sep = ""),
    #   location = unit(0, "npc"), just = "left", rot = 90,
    #   gp = gpar(fontface = "bold", color = "black")
    # ),
    "Intersection\nsize" = anno_barplot(
      x = comb_size(cm), 
      border = FALSE, 
      gp = gpar(fill = "black"),
      height = unit(4, "cm"),
      axis = TRUE,
      axis_param = list(at = c(0, 500, 1000, 1500, 2000)),
      add_numbers = TRUE, numbers_rot = 90,
      numbers_gp = gpar(fontface = "bold", color = "black")
    ),
    annotation_name_side = "left", 
    annotation_name_rot = 0
  ),
  right_annotation = upset_right_annotation(
    m = cm, 
    gp = gpar(
      fill = stringr::str_replace(string = tfData$color, pattern = "white,", replacement = "")
    ),
    annotation_name_side = "top",
    axis_param = list(side = "top"),
    axis = FALSE,
    add_numbers = TRUE, numbers_rot = 0,
    numbers_gp = gpar(fontface = "bold", color = "black")
  ),
  width = unit(12, "cm"), height = unit(3, "cm")
)

pdf(file = paste(outPrefix, ".peak_upset.pdf", sep = ""), width = 12, height = 8)
ht_upset
dev.off()

##################################################################################
statsDf <- dplyr::group_by_at(hasPeakDf, .vars = vars(starts_with("hasPeak."), "group")) %>% 
  dplyr::summarise(count = n_distinct(geneId)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(position = cumsum(count)) %>% 
  tidyr::pivot_longer(cols = starts_with("hasPeak."),
                      names_to = "sampleId", values_to = "peak") %>% 
  dplyr::mutate(
    sampleId = stringr::str_replace(
      string = sampleId, pattern = "hasPeak.", replacement = "")
  )

statsDf$sampleId <- factor(statsDf$sampleId, levels = unique(statsDf$sampleId))
statsDf$id = as.numeric(statsDf$sampleId)



pt_stats <- ggplot(data = statsDf) +
  # geom_rect(
  #   mapping = aes(ymin = id - 0.5, ymax = id + 0.5,
  #                 xmin = position - count, xmax = position,
  #                 fill = peak)) +
  # scale_y_discrete(
  #   limits = as.character(1:4),
  #   labels = structure(levels(statsDf$sampleId), names = as.character(1:4)),
  #   expand = expand_scale(add = c(0.05, 0))) +
  # scale_fill_manual(values = c("TRUE" = "#3c4ca0ff", "FALSE" = "white")) +
  geom_segment(
    mapping = aes(x = position - count, xend = position,
                  y = sampleId, yend = sampleId, color = peak),
    size = 49
  ) +
  scale_color_manual(values = c("TRUE" = "#3c4ca0ff", "FALSE" = "white")) +
  scale_y_discrete(
    limits = rev(levels(statsDf$sampleId)),
    expand = expand_scale(add = c(0.55, 0.5))
  ) +
  scale_x_continuous(
    labels = seq(from = 0, to = 5000, by = 1000),
    expand = expand_scale(add = c(0,100))) +
  labs(title = "common and unique targets") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_blank(),
    legend.position = "none",
    axis.line.x = element_line(size = 1.5),
    axis.ticks = element_line(size = 1.5),
    axis.ticks.length = unit(2, "mm"),
    plot.title = element_text(size = 24)
  )


# pdf(file = paste(outPrefix, ".peak_stats.pdf", sep = ""), width = 18, height = 6)
png(file = paste(outPrefix, ".peak_stats.png", sep = ""), width = 5000, height = 2000, res = 300)
pt_stats
dev.off()

##################################################################################










suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(summarytools))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))

## KERS complex at gene level

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")
##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
analysisName <- "KERS_complex_48h"
outDir <- here::here("analysis", "03_KERS_complex", analysisName)
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_plotSamples <- paste(outDir, "/", "samples.txt", sep = "")

regionMatType <- "deeptools"
regionMatSuffix <- "normalized_profile"
regionMatDim = c(200, 200, 100, 10)

tssMatType <- "normalizedmatrix"
tssMatSuffix <- "normalizedmatrix_3kbATG3kb"
tssMatDim <- c(300, 1, 300, 10)

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
  dataPath = TF_dataPath
)


polII_ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "polII")]
tfIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag %in% c("HA", "MYC", "TAP") & tempSInfo$TF != "untagged")]
inputIds <- tempSInfo$sampleId[which(! tempSInfo$IP_tag %in% c("polII", "HIST") & tempSInfo$TF == "untagged")]
histIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST")]
histH3Ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST" & tempSInfo$TF == "H3")]


## read the experiment sample details and select only those which are to be plotted
tfData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = tfIds,
  dataPath = TF_dataPath, profileMatrixSuffix = tssMatSuffix
)

inputData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = inputIds,
  dataPath = TF_dataPath, profileMatrixSuffix = tssMatSuffix
)

polIIData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = polII_ids,
  dataPath = polII_dataPath, profileMatrixSuffix = regionMatSuffix
)

histData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = histIds,
  dataPath = hist_dataPath, profileMatrixSuffix = tssMatSuffix
)

h3Data <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = histH3Ids,
  dataPath = hist_dataPath, profileMatrixSuffix = tssMatSuffix
)


regionProfileData <- dplyr::bind_rows(polIIData)
tssProfileData <- dplyr::bind_rows(tfData, inputData, histData, h3Data)


exptData <- dplyr::bind_rows(tfData, inputData, histData, h3Data, polIIData)

exptDataList <- purrr::transpose(exptData)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)


tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
  simplify = F, USE.NAMES = T)

##################################################################################
peakTargetMat <- peak_target_matrix(sampleInfo = tfData, position = "best")

peakTargetMat <- dplyr::mutate_at(
  .tbl = peakTargetMat,
  .vars = vars(starts_with("hasPeak.")),
  .funs = ~ factor(., levels = c(TRUE, FALSE))
)

expressionData <- get_polII_expressions(exptInfo = exptData,
                                        genesDf = geneInfo)

expressionData <- dplyr::left_join(x = expressionData, y = peakTargetMat, by = "geneId")


## genes with binding by at least one of KERS
hasPeakDf <- expressionData %>% 
  dplyr::filter_at(.vars = vars(starts_with("hasPeak")), .vars_predicate = any_vars(. == TRUE)) %>% 
  dplyr::mutate(group = group_indices(., !!! lapply(unname(tfCols$hasPeak), as.name))) %>% 
  dplyr::mutate(group = sprintf(fmt = "%02d", group))

dplyr::group_by_at(hasPeakDf, .vars = vars(starts_with("hasPeak."))) %>% 
  dplyr::summarise(n = n_distinct(geneId))

readr::write_tsv(x = hasPeakDf, file = paste(outPrefix, ".peaks_data.tab", sep = ""))

## genes bound by all the KERS members
commonPeaks <- dplyr::filter_at(
  .tbl = hasPeakDf, .vars = vars(starts_with("hasPeak.")),
  .vars_predicate = all_vars(. == TRUE)
)

dplyr::filter_at(
  hasPeakDf, .vars = vars(starts_with("hasPeak.")), .vars_predicate = all_vars(. == TRUE)
) %>% 
  dplyr::select(geneId) %>% 
  readr::write_tsv(
    file = paste(outPrefix, ".common_genes.tab", sep = "")
  )

## group label data mean profile facet plot
groupLabelDf <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(starts_with("hasPeak."), group)) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(groupLabels = paste(group, ": ", n, " genes", sep = ""))

groupLabels <- structure(groupLabelDf$groupLabels, names = groupLabelDf$group)


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

## color list
tssMatList <- NULL
regionMatList <- NULL

regionMatList <- import_profiles(
  exptInfo = regionProfileData,
  geneList = expressionData$geneId,
  source = regionMatType,
  up = regionMatDim[1], target = regionMatDim[2], down = regionMatDim[3],
  targetType = "region", targetName = "gene"
)

tssMatList <- import_profiles(
  exptInfo = tssProfileData,
  geneList = expressionData$geneId,
  up = tssMatDim[1], target = tssMatDim[2], down = tssMatDim[3],
  targetType = "point", targetName = "ATG"
)

matList <- c(regionMatList, tssMatList)

##################################################################################

## tf colors
tfMeanProfile <- NULL
tfColorList <- NULL
if(length(c(tfIds)) == 1){
  tfMeanProfile <- matList[[tfIds]]
} else{
  tfMeanProfile <- getSignalsFromList(lt = matList[tfIds])
}

quantile(tfMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# tfMeanColor <- colorRamp2(quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T), c("white", "red"))
tfColorList <- sapply(
  X = c(tfIds, inputIds),
  FUN = function(x){
    return(
      colorRamp2(breaks = quantile(tfMeanProfile, c(0.50, 0.99), na.rm = T),
                 colors = unlist(strsplit(x = exptDataList[[x]]$color, split = ","))))
  }
)

## polII colors
polIIMeanProfile <- NULL
polIIColorList <- NULL
if(nrow(polIIData) == 1){
  polIIMeanProfile <- matList[[polIIData$sampleId]]
} else{
  polIIMeanProfile <- getSignalsFromList(lt = matList[polIIData$sampleId])
}
quantile(polIIMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
polIIMeanColor <- colorRamp2(
  breaks = quantile(polIIMeanProfile, c(0.01, 0.5, 0.99), na.rm = T),
  colors = c("blue", "white", "red")
)
polIIColorList <- sapply(X = polIIData$sampleId, FUN = function(x){return(polIIMeanColor)})

colorList <- unlist(list(tfColorList, polIIColorList))


##################################################################################
## polII signal matrix
polIIMat <- data.matrix(log2(expressionData[, polII_ids] + 1))
rownames(polIIMat) <- expressionData$geneId

quantile(polIIMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.997, 0.999, 0.9999, 1), na.rm = T)


polIIHt_color <- colorRamp2(
  breaks = c(0, quantile(polIIMat, c(0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99, 0.995, 0.999))),
  colors = c("white", RColorBrewer::brewer.pal(n = 9, name = "RdPu"))
)

##################################################################################
ylimList <- list()
ylimList <- sapply(tssProfileData$sampleId, function(x){return(c(0,30))}, simplify = FALSE)


regionProfiles_peak <- multi_profile_plots(
  exptInfo = regionProfileData,
  genesToPlot = commonPeaks$geneId,
  matSource = regionMatType,
  matBins = regionMatDim,
  drawClusterAn = FALSE,
  profileColors = colorList,
  expressionColor = NULL,
  # plotExpression = showExpressionHeatmap,
  column_title_gp = gpar(fontsize = 12),
  ylimFraction = ylimList,
  showAnnotation = FALSE
)

tssProfiles_peak <- multi_profile_plots(
  exptInfo = tssProfileData,
  genesToPlot = commonPeaks$geneId,
  targetType = "point",
  targetName = "ATG",
  matSource = tssMatType,
  matBins = tssMatDim,
  clusters = NULL,
  profileColors = colorList,
  column_title_gp = gpar(fontsize = 12),
  ylimFraction = ylimList,
  showAnnotation = FALSE
)



## gene length annotation
anGl_peaks <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = commonPeaks$geneId,
  # pointSize = unit(4, "mm"),
  axis_param = list(at = c(2000, 4000), labels = c("2kb", "> 4kb"))
)

peaks_htlist <- NULL
peaks_htlist <- peaks_htlist + anGl_peaks$an
peaks_htlist <- peaks_htlist + tssProfiles_peak$heatmapList
peaks_htlist <- peaks_htlist + regionProfiles_peak$heatmapList


## make sure that the order of genes in the heatmap list and in the dataframe is same
if(all(rownames(peaks_htlist@ht_list[[ tfData$profileName[1] ]]@matrix) == commonPeaks$geneId)){
  
  rowOrd_peaks <- order(commonPeaks[[ tfCols$peakDist[[1]] ]], decreasing = TRUE)
  
}

pdfWd <- 2 + 
  (length(peaks_htlist) * 2) +
  (length(polII_ids) * 0.25 * showExpressionHeatmap) +
  2

# wd <- 500 + (nrow(exptData) * 700) + (length(polII_ids) * 500 * showExpressionHeatmap)
title_peak= paste(analysisName, ": genes bound by each KERS member", collapse = "")


pdf(file = paste(outPrefix, ".common_targets.profiles.pdf", sep = ""), width = pdfWd, height = 13)

draw(
  peaks_htlist,
  main_heatmap = exptData$profileName[1],
  # annotation_legend_list = list(profile1$legend),
  column_title = title_peak,
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_sub_title_side = "left",
  heatmap_legend_side = "right",
  gap = unit(7, "mm"),
  # row_order = rowOrd_peaks,
  # split = rep(1, nrow(commonPeaks)),
  padding = unit(rep(0.5, times = 4), "cm")
)


dev.off()



##################################################################################

## profile matrix with selected genes only
newClusters <- dplyr::select(hasPeakDf, geneId, group) %>% 
  dplyr::rename(cluster = group)

newClusters$cluster <- factor(x = newClusters$cluster, levels = sort(unique(newClusters$cluster)))

# profiles_peakGroups <- multi_profile_plots(
#   exptInfo = exptData[exptData$sampleId %in% c(tfIds, inputIds), ],
#   genesToPlot = hasPeakDf$geneId,
#   clusters = newClusters,
#   profileColors = colorList,
#   matSource = matrixType,
#   matBins = matrixDim,
#   column_title_gp = gpar(fontsize = 12)
# )


## KERS binding profile plot surrounding ATG
tssProfiles_peakGroups <- multi_profile_plots(
  exptInfo = tssProfileData,
  genesToPlot = hasPeakDf$geneId,
  targetType = "point",
  targetName = "ATG",
  matSource = tssMatType,
  matBins = tssMatDim,
  clusters = newClusters,
  profileColors = colorList,
  column_title_gp = gpar(fontsize = 12),
  # ylimFraction = ylimList
)

## polII profile over gene body
regionProfiles_peakGroups <- multi_profile_plots(
  exptInfo = regionProfileData,
  genesToPlot = hasPeakDf$geneId,
  matSource = regionMatType,
  matBins = regionMatDim,
  drawClusterAn = FALSE,
  profileColors = colorList,
  expressionColor = NULL,
  # plotExpression = showExpressionHeatmap,
  column_title_gp = gpar(fontsize = 12)
)

## polII signal heatmap
polIIMat_peakGroups <- polIIMat[hasPeakDf$geneId, , drop=FALSE]

polIIht_peakGroups <- signal_heatmap(
  log2_matrix = polIIMat_peakGroups,
  htName = "polII_signal",
  col_title = polII_ids,
  legend_title = "log2(polII_singal)",
  color = polIIHt_color,
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

peakGroups_htlist <- NULL
peakGroups_htlist <- peakGroups_htlist + tssProfiles_peakGroups$heatmapList
peakGroups_htlist <- peakGroups_htlist + polIIht_peakGroups
peakGroups_htlist <- peakGroups_htlist + grp_ht
peakGroups_htlist <- peakGroups_htlist + anGl_peakGroups$an


pdfWd <- 3 + 
  (length(tssProfiles_peakGroups$heatmapList@ht_list) * 2) +
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
testDf <- dplyr::filter(hasPeakDf, group == "01")

lineColors = purrr::map_chr(.x = exptDataList[c(tfIds, inputIds)],
                            .f = function(x) unlist(strsplit(x = x$color, split = ","))[2])

lineShape = structure(.Data = c(1, 1, 1, 1, 6),
                      names = exptData$sampleId[exptData$sampleId %in% c(tfIds, inputIds)])


draw_avg_profile_plot(
  exptInfo = exptData[exptData$sampleId %in% c(tfIds, inputIds), ],
  profileMats = matList[c(tfIds, inputIds)],
  genes = testDf$geneId,
  lineColors = lineColors,
  lineShape = lineShape
)



## average profile data for each group
groupMeanProfiles <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(starts_with("hasPeak."), group)) %>%
  dplyr::do(
    geneset_average_profile(
      exptInfo = exptData[exptData$sampleId %in% c(tfIds, inputIds), ],
      profileMats = matList[c(tfIds, inputIds)],
      genes = .$geneId,
      cluster = unique(.$group)
    )
  ) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()

## decide the axis labels
axisBrk <- NULL
axisLab <- NULL
targetEnd <- NULL
profileAttributes <- attributes(matList[[exptData$sampleId[1]]])

if(profileAttributes$target_is_single_point){
  axisBrk <- c(
    profileAttributes$upstream_index[1],
    profileAttributes$target_index[1],
    tail(profileAttributes$downstream_index, 1)
  )
  
  axisLab <- c(
    -profileAttributes$extend[1],
    profileAttributes$target_name,
    profileAttributes$extend[2]
  )
  
  targetEnd <- axisBrk[2] + 1
  
} else if(! profileAttributes$target_is_single_point){
  axisBrk <- c(
    profileAttributes$upstream_index[1],
    profileAttributes$target_index[1],
    tail(profileAttributes$target_index, 1),
    tail(profileAttributes$downstream_index, 1)
  )
  
  axisLab <- c(
    -profileAttributes$extend[1],
    "START", "END",
    profileAttributes$extend[2]
  )
  
  targetEnd <- axisBrk[3]
}

## plot the groups using facets
pt_grpLines = ggplot(data = groupMeanProfiles) +
  geom_line(
    mapping = aes(x = bin, y = mean, group = sampleId, color = sampleId, linetype = sampleId),
    size = 0.8, alpha = 0.8
  ) +
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
pt_grpLines
dev.off()


##################################################################################
## upset plot
cm = make_comb_mat(
  dplyr::select(hasPeakDf, tfCols$hasPeak) %>% 
    dplyr::mutate_all(.funs = ~ as.logical(.))
)

ht_upset <- UpSet(
  cm,
  pt_size = unit(4, "mm"), lwd = 2,
  set_order = tfData$sampleId,
  row_labels = tfData$TF,
  row_names_gp = gpar(fontsize = 14, fontface = "italic"),
  column_title = paste("KERS bound genes statistics", analysisName),
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

pdf(file = paste(outPrefix, ".peak_stats_upset.pdf", sep = ""), width = 12, height = 8)
ht_upset
dev.off()


##################################################################################
## peak annotation stats

peakAnnotations <- dplyr::select(commonPeaks, geneId, contains(!!tfIds[1])) %>% 
  dplyr::left_join(
    y = suppressMessages(readr::read_tsv(file = exptDataList[[tfIds[1]]]$peakAnno)) %>% 
      dplyr::select(peakId, summitDist),
    by = setNames(c("peakId"), tfCols$peakId[tfIds[1]])
  ) %>% 
  dplyr::mutate(
    newCategory = dplyr::case_when(
      !!sym(tfCols$peakType[tfIds[1]]) == "5UTR" ~ "promoter",
      !!sym(tfCols$peakType[tfIds[1]]) == "tx_start" ~ "promoter",
      !!sym(tfCols$peakType[tfIds[1]]) == "3UTR" ~ "nearEnd",
      !!sym(tfCols$peakType[tfIds[1]]) == "tx_end" ~ "nearEnd",
      !!sym(tfCols$peakType[tfIds[1]]) == "EXON" ~ "Exon/Intron",
      !!sym(tfCols$peakType[tfIds[1]]) == "INTRON" ~ "Exon/Intron",
      !!sym(tfCols$peakType[tfIds[1]]) == "include_tx" ~ "featureInPeak",
      !!sym(tfCols$peakType[tfIds[1]]) == "upstream" ~ "upstream",
      TRUE ~ !!sym(tfCols$peakType[tfIds[1]])
    ),
    newCategory = dplyr::case_when(
      summitDist > -1000 & summitDist < 0 ~ "promoter",
      TRUE ~ newCategory
    )
  ) %>% 
  dplyr::mutate(
    newCategory = forcats::fct_relevel(
      .f = newCategory, "upstream", "promoter", "Exon/Intron", "nearEnd"
    )
  )


# peakAnnotations <- dplyr::left_join(
#   x = peakAnnotations,
#   y = suppressMessages(readr::read_tsv(file = exptDataList[[tfIds[1]]]$peakAnno)) %>% 
#     dplyr::select(peakId, summitDist),
#   by = setNames(c("peakId"), tfCols$peakId[tfIds[1]])
# )


pt_distDensity <- dplyr::filter(peakAnnotations, newCategory %in% c("promoter", "upstream")) %>% 
  dplyr::mutate(
    summitDist = pmax(summitDist, -2000)
  ) %>% 
  ggplot2::ggplot(
    mapping = aes(x = summitDist)
  ) +
  geom_density(size = 2) +
  labs(
    title = "Summit distance density distribution for promoter and upstream peaks",
    x = "Distance to TSS",
    y = "Density"
  ) +
  coord_cartesian(xlim = c(-2000, 100)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    panel.grid = element_blank()
  )

pdf(file = paste(outPrefix, ".summit_dist_density.pdf", sep = ""), width = 7, height = 6)
pt_distDensity
dev.off()


pieDf <- dplyr::group_by(peakAnnotations, newCategory) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::arrange(desc(newCategory)) %>% 
  dplyr::mutate(
    cumsum = rev(cumsum(rev(count))),
    pos = count/2 + lead(cumsum, 1),
    pos = if_else(is.na(pos), count/2, pos)
  ) %>% 
  dplyr::mutate(
    label = round(count / sum(count), digits = 4),
    ypos = cumsum(count)- 0.5*count
  )

pt_pie <- ggplot2::ggplot(
  data = pieDf,
  mapping = aes(x = 1, y = count, fill = newCategory)
) +
  geom_bar(stat = "identity", color = "black", size = 1) +
  geom_text_repel(
    mapping = aes(x = 1.4, y = ypos, label = scales::percent(label)),
    hjust = 0.5, direction = "x", xlim = 1,
    segment.size = unit(1, "cm"),
    nudge_x = 0.3, size = 8, show.legend = FALSE
  ) +
  scale_fill_discrete(
    breaks = c("upstream", "promoter", "Exon/Intron", "nearEnd"),
    labels = c("upstream (< -1kb)", "promoter (-1kb to START)", "Exon/Intron", "after STOP")
  ) +
  scale_x_continuous(limits = c(0, 1.5), expand = expansion(add = c(0, 0))) +
  coord_polar(theta = "y", start = 0, direction = -1) +
  labs(title = "Peak annotation distribution") +
  theme_void() +
  theme(
    plot.title = element_text(size=30, face="bold", hjust = 0.5),
    legend.key.size = unit(1, "cm"),
    # legend.position = "bottom",
    # legend.direction = "vertical",
    legend.text = element_text(size=20, face="bold"),
    legend.title = element_blank()
  )

pdf(file = paste(outPrefix, ".peak_annotation_pie.pdf", sep = ""), width = 10, height = 7)
pt_pie
dev.off()


##################################################################################

suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(summarytools))


## 1) compares the binding of two TF pairs
## 2) perform GO enrichment, GO grouping, KEGG enrichment for different groups
## 3) generate TF binding change plot for TF pairs w.r.t. polII fold change in 48h vs 20h data


rm(list = ls())

source(file = "https://raw.githubusercontent.com/lakhanp1/omics_utils/master/04_GO_enrichment/s01_enrichment_functions.R")

##################################################################################
## main configuration
comparisonName <- "kdmB_48h_vs_20h.2"
outDir <- here::here("analysis", "05_KERS_48h_vs_20h", comparisonName)
outPrefix <- paste(outDir, "/", comparisonName, sep = "")

file_polIIRatioConf <- here::here("data", "reference_data", "polII_ratio.config.tab")

mainPolIIPair <- "polII_untagged.48h_vs_20h"
polIIPairs <- c(
  mainPolIIPair, "polII_kdmB_del.48h_vs_20h", "polII_20h.kdmB_del_vs_untagged",
  "polII_48h.kdmB_del_vs_untagged"
)

tfDiffPairs <- list(
  p1 = list(
    name = "kdmB_48h_vs_20h",
    samples = c("An_kdmB_20h_HA_1", "An_kdmB_48h_HA_1")
    # ),
    # p2 = list(
    #   name = "rpdA_48h_vs_20h",
    #   samples = c("An_rpdA_20h_HA_1", "An_rpdA_48h_HA_1")
    # ),
    # p3 = list(
    #   name = "sntB_48h_vs_20h",
    #   samples = c("An_sntB_20h_HA_1", "An_sntB_48h_HA_1")
    # ),
    # p4 = list(
    #   name = "ecoA_48h_vs_20h",
    #   samples = c("An_ecoA_20h_HA_1", "An_ecoA_48h_HA_1")
  )
)

mainTfPair <- "p1"

otherHist <- c("An_H3_20h_HIST_1", "An_H3_48h_HIST_1")

# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
matrixType <- "normalizedmatrix_5kb"
matrixDim = c(500, 200, 100, 10)

clusterStorePath <- paste(outPrefix, ".profile_kmeans_clusters.txt", sep = "")

## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/annotation_resources/geneid2go.ANidulans.topGO.map"

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db

showExpressionHeatmap <- FALSE

outPrefix_tfPolII <- paste(outPrefix, ".tfPolII", sep = "")
outPrefix_tfSpecific <- paste(outPrefix, ".specific_binding", sep = "")
outPrefix_peakExp <- paste(outPrefix, ".pkExpGenes", sep = "")

anLables <- list()
anLables[["is_SM_gene"]] <- "SM gene"
anLables[["is_TF"]] <- "Transcription Factor"
anLables[["gene_length"]] <- "Gene Length"


##################################################################################

## genes to read and annotations
geneSet <- suppressMessages(
  readr::read_tsv(file = file_genes, col_names = c("chr", "start", "end", "geneId", "score", "strand"))
)

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
  dplyr::filter(comparison %in% polIIPairs)

polIIDiffPairs <- purrr::transpose(polIIRatioConf)  %>% 
  purrr::set_names(nm = purrr::map(., "comparison"))

polII1 <- polIIDiffPairs[[mainPolIIPair]]$group1
polII2 <- polIIDiffPairs[[mainPolIIPair]]$group2


tf1 <- tfDiffPairs[[mainTfPair]]$samples[1]
tf2 <- tfDiffPairs[[mainTfPair]]$samples[2]


##################################################################################

## read the experiment sample details and select only those which are to be plotted

tfData <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = unique(purrr::map(tfDiffPairs, "samples") %>% unlist() %>% unname()),
  dataPath = TF_dataPath,
  profileMatrixSuffix = matrixType
)

polIIData <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = unique(purrr::map(polIIDiffPairs, `[`, c("group2", "group1")) %>% unlist() %>% unname()),
  dataPath = polII_dataPath,
  profileMatrixSuffix = "normalized_profile"
)

histData <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = otherHist,
  dataPath = hist_dataPath,
  profileMatrixSuffix = matrixType
)

exptData <- dplyr::bind_rows(tfData, polIIData, histData)
exptDataList <- purrr::transpose(exptData)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

polII_ids <- exptData$sampleId[which(exptData$IP_tag == "polII")]
tfIds <- exptData$sampleId[which(exptData$IP_tag %in% c("HA", "MYC", "TAP") & exptData$TF != "untagged")]


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


expressionData <- get_polII_expressions(exptInfo = exptData,
                                        genesDf = geneInfo)

## add fold change columns
for (i in names(polIIDiffPairs)) {
  expressionData <- get_fold_change(
    df = expressionData,
    nmt = polIIDiffPairs[[i]]$group1,
    dmt = polIIDiffPairs[[i]]$group2,
    newCol = polIIDiffPairs[[i]]$comparison,
    isExpressedCols = polIICols$is_expressed)
}


expressionData <- get_TF_binding_data(exptInfo = tfData,
                                      genesDf = expressionData)


expressionData %>% 
  dplyr::select(geneId, starts_with("hasPeak")) %>% 
  dplyr::group_by_at(.vars = vars(starts_with("hasPeak"))) %>% 
  dplyr::summarise(n = n())


expressionData <- dplyr::group_by_at(expressionData, unname(tfCols$hasPeak[c(tf1, tf2)])) %>%
  dplyr::mutate(
    group = paste(sprintf(fmt = "%02d", cur_group_id()), " (n=", n(), ")", sep = "")
  ) %>% 
  dplyr::ungroup()


readr::write_tsv(x = expressionData, file = paste(outPrefix, ".data.tab", sep = ""))

glimpse(expressionData)

## peak data
hasPeakDf <- dplyr::filter_at(
  .tbl = expressionData,
  .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)])),
  .vars_predicate = any_vars(. == "TRUE")
)

dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)]), group)) %>%
  dplyr::summarise(n = n())

## check the MA plot
plot_MA_gg(df = expressionData, s1 = polII2, s2 = polII1, title = "MA plot", colorCol = "group")
# plot_MA_gg(df = expressionData, s1 = polII1, s2 = polII2, title = "MA plot", colorCol = "group")
plot_MA_gg(df = expressionData, s1 = polIIDiffPairs[[2]]$group2,
           s2 = polIIDiffPairs[[2]]$group1, title = "MA plot", colorCol = "group")


plot_scatter(df = expressionData, s1 = polII1, s2 = polII2, title = "Scatter plot", colorCol = "group")
plot_scatter(df = expressionData, s1 = polII2, s2 = polII1, title = "Scatter plot", colorCol = "group")

peakCovDf <- tidyr::gather(data = hasPeakDf, key = "sampleId", value = "coverage", starts_with("peakCoverage.")) %>% 
  dplyr::select(geneId, sampleId, coverage) %>% 
  dplyr::mutate(
    sampleId = gsub(pattern = "peakCoverage.", replacement = "", x = sampleId, fixed = T),
    tf = gsub(pattern = "An_(\\w+)_(20h|48h)_HA_1", replacement = "\\1", x = sampleId, perl = T),
    time = gsub(pattern = "An_(\\w+)_(20h|48h)_HA_1", replacement = "\\2", x = sampleId, perl = T)) %>% 
  as_tibble()

ggplot(data = peakCovDf,
       mapping = aes(x = coverage, group = sampleId, color = tf, linetype = time)) +
  # geom_density(size = 1) +
  geom_line(stat = "density", size = 1) +
  coord_cartesian(xlim = c(0, 100)) +
  facet_wrap(tf ~ ., nrow = 2, ncol = 2)


##################################################################################

## join the profile matrices and do clustering
# mergedKm = merged_profile_matrix_cluster(name = comparisonName,
#                                          exptInfo = exptData,
#                                          genes = geneSet$geneId,
#                                          clusterStorePath = clusterStorePath,
#                                          k = 12)

# clusterInfo <- fread(file = clusterStorePath, sep = "\t", header = T, stringsAsFactors = F)
# 
# expressionData <- dplyr::left_join(x = expressionData, y = clusterInfo, by = c("geneId" = "geneId"))


##################################################################################
## topGO enrichment
goEnrich <- dplyr::group_by_at(
  .tbl = hasPeakDf, .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)]), group)
) %>%
  do(
    topGO_enrichment(
      goMapFile = file_topGoMap, genes = .$geneId, goNodeSize = 5,
      orgdb = orgDb, geneNameColumn = "GENE_NAME"
    )
  )

readr::write_tsv(x = goEnrich, file = paste(outPrefix, ".topGO.tab", sep = ""))


## pathway enrichment
keggEnr <- dplyr::group_by_at(
  .tbl = hasPeakDf, .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)]), group)
) %>% 
  do(
    keggprofile_enrichment(
      genes = .$geneId, orgdb = orgDb, keytype = "GID", keggIdColumn = "KEGG_ID",
      keggOrg = "ani", pvalCut = 0.05, geneNameColumn = "GENE_NAME"
    )
  )

readr::write_tsv(x = keggEnr, file = paste(outPrefix, ".KEGG.tab", sep = ""))


##################################################################################

## polII signal matrix
polIIMat <- as.matrix(log2(expressionData[, polII_ids] + 1))
rownames(polIIMat) <- expressionData$geneId

quantile(polIIMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.997, 0.999, 0.9999, 1), na.rm = T)


polII_color <- colorRamp2(
  breaks = c(0, quantile(polIIMat, c(0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99, 0.995, 0.999))),
  colors = c("white", RColorBrewer::brewer.pal(n = 9, name = "RdPu"))
)


## polII signal fold change matrix
lfcMat <- as.matrix(expressionData[, purrr::map_chr(polIIDiffPairs, "comparison"), drop = FALSE])
rownames(lfcMat) <- expressionData$geneId

lfc_color <- colorRamp2(
  breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2),
  colors = RColorBrewer::brewer.pal(n = 7, name = "PuOr")
)

##################################################################################
## colors for profile matrix
matList <- import_profiles(
  exptInfo = exptData,
  geneList = geneInfo$geneId,
  source = matrixType,
  up = matrixDim[1], target = matrixDim[2], down = matrixDim[3]
)


## tf colors
tfMeanProfile <- NULL
if(length(c(tfIds)) == 1){
  tfMeanProfile <- matList[[tfIds]]
} else{
  tfMeanProfile <- getSignalsFromList(lt = matList[tfData$sampleId])
}

quantile(tfMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
tfMeanColorList <- sapply(
  X = tfData$sampleId,
  FUN = function(x){
    return(
      colorRamp2(
        breaks = quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T),
        colors = unlist(strsplit(x = exptDataList[[x]]$color, split = ","))
      )
    )
  }
)

tfWiseColors <- sapply(
  X = tfData$sampleId,
  FUN = function(x){
    cat(x, "\n")
    print(quantile(
      matList[[x]],
      c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T))
    
    return(
      colorRamp2(
        breaks = quantile(matList[[x]], c(0.50, 0.99), na.rm = T),
        colors = unlist(strsplit(x = exptDataList[[x]]$color, split = ","))
      )
    )
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
polIIMeanColor <- colorRamp2(quantile(polIIMeanProfile, c(0.01, 0.5, 0.995), na.rm = T), c("blue", "white", "red"))
polIIColorList <- sapply(X = polIIData$sampleId, FUN = function(x){return(polIIMeanColor)})

## histone colors
histMeanProfile <- NULL
histColorList <- NULL
if(nrow(histData) == 1){
  histMeanProfile <- matList[[histData$sampleId]]
} else{
  histMeanProfile <- getSignalsFromList(lt = matList[histData$sampleId])
}

histColorList <- sapply(
  X = histData$sampleId,
  FUN = function(x){
    return(
      colorRamp2(
        breaks = quantile(histMeanProfile, c(0.20, 0.995), na.rm = T),
        colors = unlist(strsplit(x = exptDataList[[x]]$color, split = ","))
      )
    )
  }
)

## LFC(TF2/TF1) matrix
## TF1 scalled matrix
tf1ScalledMat <- scale(x = matList[[tf1]], center = TRUE, scale = TRUE)
quantile(tf1ScalledMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
tf1ScalledColor <- colorRamp2(
  breaks = quantile(tf1ScalledMat, c(0.50, 0.99), na.rm = T),
  colors = unlist(strsplit(x = exptDataList[[tf1]]$color, split = ","))
)

## TF2 scalled matrix
tf2ScalledMat <- scale(x = matList[[tf2]], center = TRUE, scale = TRUE)
quantile(tf2ScalledMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
tf2ScalledColor <- colorRamp2(
  breaks = quantile(tf2ScalledMat, c(0.50, 0.99), na.rm = T),
  colors = unlist(strsplit(x = exptDataList[[tf2]]$color, split = ","))
)

## Difference between TF2 and TF1 scalled matrix
scalledTfDiffMat <- tf2ScalledMat - tf1ScalledMat
plot(density(scalledTfDiffMat))
quantile(
  x = scalledTfDiffMat,
  c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.05, seq(0.1, 0.9, by = 0.1),
    0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T
)


scalledTfDiffColor <- colorRamp2(
  breaks = c(-3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3),
  colors = RColorBrewer::brewer.pal(n = 11, name = "PuOr")
)

##################################################################################
## genes which have peak in both TFs and polII signal in either one or both TFs

peakExpDf <- dplyr::filter_at(
  .tbl = expressionData, .vars = unname(tfCols$hasPeak[c(tf1, tf2)]),
  .vars_predicate = all_vars(. == TRUE)
) %>% 
  dplyr::filter_at(
    .vars = unname(polIICols$is_expressed[c(polII1, polII2)]),
    .vars_predicate = any_vars(. == TRUE)
  ) 

## group genes into polII fold change bins
peakExpDf <- dplyr::mutate(
  peakExpDf,
  group = case_when(
    !! as.name(polIIDiffPairs[[mainPolIIPair]]$comparison) >= 2 ~ "G1: LFC >= 2",
    !! as.name(polIIDiffPairs[[mainPolIIPair]]$comparison) >= 1 ~ "G2: 1 <= LFC < 2",
    !! as.name(polIIDiffPairs[[mainPolIIPair]]$comparison) >= 0.5 ~ "G3: 0.5 <= LFC < 1",
    !! as.name(polIIDiffPairs[[mainPolIIPair]]$comparison) >= 0 ~ "G4: 0 <= LFC < 0.5",
    !! as.name(polIIDiffPairs[[mainPolIIPair]]$comparison) <= -2 ~ "G8: -2 >= LFC",
    !! as.name(polIIDiffPairs[[mainPolIIPair]]$comparison) <= -1 ~ "G7: -2 < LFC <= -1",
    !! as.name(polIIDiffPairs[[mainPolIIPair]]$comparison) <= -0.5 ~ "G6: -1 < LFC <= -0.5",
    !! as.name(polIIDiffPairs[[mainPolIIPair]]$comparison) < 0 ~ "G5: -0.5 < LFC < 0",
    TRUE ~ "0"
  )
) %>% 
  dplyr::group_by(group) %>% 
  dplyr::mutate(
    group = paste(group, "\n(n=", n(), ")", sep = "")
  ) %>% 
  dplyr::ungroup()
  

plot_scatter(df = peakExpDf, s1 = polII1, s2 = tfCols$peakCoverage[tf2], title = "Scatter plot", colorCol = "group")
plot_MA_gg(df = peakExpDf, s1 = tfCols$peakCoverage[tf2], s2 = tfCols$peakCoverage[tf1], colorCol = "group")
plot_scatter(df = peakExpDf, s1 = tfCols$peakCoverage[tf1], s2 = tfCols$peakCoverage[tf2], title = "Scatter plot", colorCol = "group")

dplyr::group_by(peakExpDf, group) %>% 
  dplyr::summarise(n = n()) %>% 
  readr::write_tsv(file = paste(outPrefix_peakExp, ".stats.tab", sep = ""))


peakExp_clusters <- dplyr::select(peakExpDf, geneId, cluster = group)


## profile heatmap: no need to include input
multiProf_peakExp <- multi_profile_plots(
  exptInfo = tfData[tfData$sampleId %in% tfIds, ],
  genesToPlot = peakExpDf$geneId,
  # matSource = matrixType,
  matBins = matrixDim,
  clusters = peakExp_clusters,
  profileColors = tfWiseColors,
  # ylimFraction = ylimList,
  column_title_gp = gpar(fontsize = 12)
)

## polII profile plots
polIIProf_peakExp <- multi_profile_plots(
  exptInfo = polIIData,
  genesToPlot = peakExpDf$geneId,
  clusters = peakExp_clusters,
  drawClusterAn = FALSE,
  matSource = "deeptools",
  matBins = c(200, 200, 100, 10),
  profileColors = polIIColorList,
  column_title_gp = gpar(fontsize = 12)
)

## Scalled TF diff profile
scalledTfDiffProf_peakExp <- profile_heatmap(
  profileMat = scalledTfDiffMat[peakExpDf$geneId, ],
  signalName = comparisonName,
  profileColor = scalledTfDiffColor,
  column_title_gp = gpar(fontsize = 12),
  geneGroups = peakExp_clusters,
  ylimFraction = c(-2.3, 1.7)
)


## polII signal heatmap
polIIMat_peakExp <- polIIMat[peakExpDf$geneId, polII_ids]

polIIht_peakExp <- signal_heatmap(
  log2_matrix = polIIMat_peakExp,
  htName = "polII_exp",
  col_title = polII_ids,
  legend_title = "log2(polII_singal)",
  color = polII_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE
)

## polII fold change heatmap
lfc_peakExp <- lfcMat[peakExpDf$geneId, ]

lfcHt_peakExp <- signal_heatmap(
  log2_matrix = lfc_peakExp,
  htName = "polII_lfc",
  legend_title = "log2(fold change)",
  color = lfc_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE
)

## gene length annotation
gl_peakExp <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = peakExpDf$geneId,
  axis_param = list(at = c(0, 2000, 4000), labels = c("0kb", "2kb", ">4kb"))
)

htlist_peakExp <- gl_peakExp$an +
  polIIht_peakExp +
  lfcHt_peakExp + 
  multiProf_peakExp$heatmapList +
  scalledTfDiffProf_peakExp$heatmap +
  polIIProf_peakExp$heatmapList



if( all(rownames(multiProf_peakExp$profileHeatmaps[[1]]$heatmap@matrix) == peakExpDf$geneId) ){
  ## row order by polII LFC
  rowOrd_peakExp <- order(peakExpDf[[ polIIDiffPairs[[mainPolIIPair]]$comparison ]], decreasing = TRUE)
}


title_peakExp <- paste("TF binding vs polII transcription: genes bound both at 20h and 48h and atleast one timepoint showing polII signal (n=", nrow(peakExpDf), ")", sep = "")

pdfWd <- 4 + (length(multiProf_peakExp$profileHeatmaps) * 3) + 3 +
  (length(polIIProf_peakExp$profileHeatmaps) * 3) +
  (length(polII_ids) * 0.5 * showExpressionHeatmap) + 
  (length(polIIDiffPairs) * 0.5)

# draw Heatmap and add the annotation name decoration
pdf(file = paste(outPrefix_peakExp, ".profile.pdf", sep = ""), width = pdfWd, height = 12)

draw(htlist_peakExp,
     main_heatmap = tfData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_peakExp,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     gap = unit(6, "mm"),
     # heatmap_legend_side = "bottom",
     # annotation_legend_side = "bottom",
     row_order = rowOrd_peakExp,
     padding = unit(rep(0.5, times = 4), "cm")
)


dev.off()


## box plots
boxDt <- tidyr::pivot_longer(
  data = peakExpDf,
  cols = c(starts_with("hasPeak."), starts_with("peakCoverage.")),
  names_to = c(".value", "sampleId"),
  names_sep = "\\."
)

boxDt$sampleId <- forcats::fct_relevel(.f = boxDt$sampleId, names(tfCols$hasPeak))

boxDt <- dplyr::left_join(
  x = boxDt, y = dplyr::select(tfData, sampleId, TF, timepoint), by = "sampleId"
) %>% 
  dplyr::group_by(sampleId) %>% 
  dplyr::mutate(coverageRank = rank(peakCoverage)) %>% 
  dplyr::ungroup() %>% 
  as.data.frame()


pt <- ggplot(
  data = boxDt, mapping = aes(x = group, y = peakCoverage)) +
  geom_boxplot(mapping = aes(fill = timepoint), alpha = 1, size = 0.5) +
  geom_point(mapping = aes(color = timepoint), shape = 16, size = 1, position=position_jitterdodge(0.2)) +
  scale_color_manual(values = c("20h" = "#b35806", "48h" = "#542788")) +
  scale_fill_manual(values = c("20h" = "#fee0b6", "48h" = "#d8daeb")) +
  # scale_y_continuous(trans = "log2") +
  facet_wrap(. ~ TF, nrow = 1, ncol = 4) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
    axis.text.y = element_text(size = 14)
  )


pdf(file = paste(outPrefix_peakExp, ".box.pdf", sep = ""), width = 8, height = 12)
pt
dev.off()


##################################################################################
## profile plot with genes showing peak specific and common between the conditions
tfSpecificDf <- hasPeakDf

clusters_tfSpecific <- dplyr::select(tfSpecificDf, geneId, cluster = group)

ylimList <- sapply(X = tfData$sampleId, FUN = function(x){c(0, 30)}, USE.NAMES = T, simplify = F)

## TF profile plot
multiProf_tfSpecific <- multi_profile_plots(
  exptInfo = tfData[tfData$sampleId %in% tfIds, ],
  genesToPlot = tfSpecificDf$geneId,
  matSource = matrixType,
  matBins = matrixDim,
  profileColors = tfMeanColorList,
  clusters = clusters_tfSpecific,
  column_title_gp = gpar(fontsize = 12),
  ylimFraction = ylimList
)

## polII profile plots
polIIProf_tfSpecific <- multi_profile_plots(
  exptInfo = polIIData,
  genesToPlot = tfSpecificDf$geneId,
  clusters = clusters_tfSpecific,
  drawClusterAn = FALSE,
  matSource = "deeptools",
  matBins = c(200, 200, 100, 10),
  profileColors = polIIColorList,
  column_title_gp = gpar(fontsize = 12)
)

## Scalled TF diff profile
scalledTfDiffProf_tfSpecific <- profile_heatmap(
  profileMat = scalledTfDiffMat[tfSpecificDf$geneId, ],
  signalName = comparisonName,
  profileColor = scalledTfDiffColor,
  column_title_gp = gpar(fontsize = 12),
  geneGroups = clusters_tfSpecific
)

## histone profile plots
histProfiles <- multi_profile_plots(
  exptInfo = histData,
  genesToPlot = tfSpecificDf$geneId,
  matSource = matrixType,
  matBins = matrixDim,
  profileColors = histColorList,
  clusters = clusters_tfSpecific,
  column_title_gp = gpar(fontsize = 12),
  drawClusterAn = FALSE
)

## polII signal heatmap
polIIMat_tfSpecific <- polIIMat[tfSpecificDf$geneId, polII_ids]

polIIHt_tfSpecific <- signal_heatmap(
  log2_matrix = polIIMat_tfSpecific,
  htName = "polII_exp",
  col_title = polII_ids,
  legend_title = "log2(polII_singal)",
  color = polII_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE
)

## polII signal fold change heatmap
lfc_tfSpecific <- lfcMat[tfSpecificDf$geneId, ]

lfc_heatmap <- signal_heatmap(
  log2_matrix = lfc_tfSpecific,
  htName = "polII_lfc",
  legend_title = "log2(fold change)",
  color = lfc_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE
)


gl_tfSpecific <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = tfSpecificDf$geneId,
  axis_param = list(at = c(0, 2000, 4000), labels = c("0kb", "2kb", ">4kb"))
)

htlist_tfSpecific <- gl_tfSpecific$an + 
  multiProf_tfSpecific$heatmapList + 
  scalledTfDiffProf_tfSpecific$heatmap +
  polIIHt_tfSpecific + 
  lfc_heatmap +
  polIIProf_tfSpecific$heatmapList +
  histProfiles$heatmapList



if( all(rownames(multiProf_tfSpecific$profileHeatmaps[[1]]$heatmap@matrix) == tfSpecificDf$geneId) ){
  rowOrd_tfSpecific <- order(tfSpecificDf[[ tfCols$peakDist[tf1] ]], tfSpecificDf[[ tfCols$peakDist[tf2] ]],
                             decreasing = TRUE)
}


title_tfSpecific <- "Differential binding of kdmB at 20h and 48h"

pdfWd <- 2 + (length(multiProf_tfSpecific$profileHeatmaps) * 2) + 3 +
  (length(polII_ids) * 0.5 * showExpressionHeatmap) + 
  (length(polIIDiffPairs) * 0.5) + 
  (length(polIIProf_tfSpecific$profileHeatmaps) * 2) +
  (length(histProfiles$profileHeatmaps) * 2) + 3

# draw Heatmap and add the annotation name decoration
pdf(file = paste0(outPrefix_tfSpecific, ".profile.pdf", collapse = ""), width = pdfWd, height = 12)
draw(htlist_tfSpecific,
     main_heatmap = tfData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_tfSpecific,
     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     row_sub_title_side = "left",
     gap = unit(6, "mm"),
     row_order = rowOrd_tfSpecific,
     padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()


##################################################################################
## genes with either peak detected or with top 10% polII signal 

tfPolIIDf <- dplyr::select(expressionData, -group) %>% 
  dplyr::filter_at(
    .vars = unname(c(tfCols$hasPeak[c(tf1, tf2)], polIICols$is_expressed[c(polII1, polII2)])),
    .vars_predicate = any_vars(. == TRUE)
  ) %>% 
  dplyr::group_by_at(
    .vars = vars(unname(c(tfCols$hasPeak[c(tf1, tf2)], polIICols$is_expressed[c(polII1, polII2)])))
  ) %>% 
  dplyr::mutate(
    group = paste(sprintf(fmt = "%02d", cur_group_id()), " (n=", n(), ")", sep = "")
  ) %>% 
  dplyr::ungroup()


tfpolII_clusters <- dplyr::select(tfPolIIDf, geneId, cluster = group)

tfPolIIDf %>% 
  dplyr::group_by_at(
    .vars = vars(unname(c(tfCols$hasPeak[c(tf1, tf2)], polIICols$is_expressed[c(polII1, polII2)])), "group")
  ) %>% 
  dplyr::summarise(n = n())

## profile heatmap
multiProf_tfPolII <- multi_profile_plots(
  exptInfo = tfData[tfData$sampleId %in% tfIds, ],
  genesToPlot = tfPolIIDf$geneId,
  matSource = matrixType,
  matBins = matrixDim,
  profileColors = tfMeanColorList,
  column_title_gp = gpar(fontsize = 12),
  clusters = tfpolII_clusters
)

## polII profile plots
polIIProf_tfPolII <- multi_profile_plots(
  exptInfo = polIIData,
  genesToPlot = tfPolIIDf$geneId,
  clusters = tfpolII_clusters,
  drawClusterAn = FALSE,
  matSource = "deeptools",
  matBins = c(200, 200, 100, 10),
  profileColors = polIIColorList,
  column_title_gp = gpar(fontsize = 12)
)

## Scalled TF diff profile
scalledTfDiffProf_tfPolII <- profile_heatmap(
  profileMat = tfLfcMat[tfPolIIDf$geneId, ],
  signalName = comparisonName,
  profileColor = scalledTfDiffColor,
  geneGroups = tfpolII_clusters
)

## polII signal heatmap
polIIMat_tfPolII <- polIIMat[tfPolIIDf$geneId, polII_ids]

polIIht_tfPolII <- signal_heatmap(
  log2_matrix = polIIMat_tfPolII,
  htName = "polII_exp",
  col_title = polII_ids,
  legend_title = "log2(polII_singal)",
  color = polII_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE
)

## polII fold change heatmap
lfc_tfPolII <- lfcMat[tfPolIIDf$geneId, ]

lfcHt_tfPolII <- signal_heatmap(
  log2_matrix = lfc_tfPolII,
  htName = "polII_lfc",
  legend_title = "log2(fold change)",
  color = lfc_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE
)

## gene length annotation
gl_tfPolII <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = tfPolIIDf$geneId,
  axis_param = list(at = c(0, 2000, 4000), labels = c("0kb", "2kb", ">4kb"))
)

htlist_tfPolII <- gl_tfPolII$an +
  multiProf_tfPolII$heatmapList +
  scalledTfDiffProf_tfPolII$heatmap +
  polIIht_tfPolII + lfcHt_tfPolII +
  polIIProf_tfPolII$heatmapList



if( all(rownames(multiProf_tfPolII$profileHeatmaps[[1]]$heatmap@matrix) == tfPolIIDf$geneId) ){
  ## order by peak distance
  rowOrd_tfPolII <- order(
    tfPolIIDf[[ tfCols$peakDist[tf1] ]],
    tfPolIIDf[[ polIIDiffPairs[[mainPolIIPair]]$comparison ]],
    decreasing = TRUE
  )
}


title_tfPolII <- "Differential binding of TF at 20h and 48h: genes with macs2 peaks or top 10% polII signal"


pdfWd <- 2 + (length(multiProf_tfPolII$profileHeatmaps) * 3.5) + 3 +
  (length(polII_ids) * 0.5 * showExpressionHeatmap) + 
  (length(polIIProf_tfPolII$profileHeatmaps) * 3.5) +
  (length(polIIDiffPairs) * 0.5)

# draw Heatmap and add the annotation name decoration
pdf(file = paste0(outPrefix_tfPolII, ".profile.pdf", collapse = ""), width = pdfWd, height = 12)

draw(htlist_tfPolII,
     main_heatmap = tfData$profileName[1],
     column_title = title_tfPolII,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     gap = unit(6, "mm"),
     row_order = rowOrd_tfPolII,
     padding = unit(rep(0.5, times = 4), "cm")
)

dev.off()

##################################################################################



suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(BSgenome.Anidulans.FGSCA4.AspGD))
suppressPackageStartupMessages(library(viridis))


## This script plots the profile heatmaps for multiple samples. If the expression values are present, it
## also plots a simple heatmap using these expression values

rm(list = ls())

source("D:/work_lakhan/github/omics_utils/04_GO_enrichment/s01_enrichment_functions.R")
##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
comparisonName <- "KERS_20h_PIC"
outDir <- here::here("analysis", "12_KERS_PIC", comparisonName)

file_plotSamples <- file.path(outDir, "samples.txt")

outPrefix <- file.path(outDir, comparisonName)

# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
regionMatType <- "deeptools"
regionMatSuffix <- "normalized_profile"
regionMatDim = c(200, 200, 100, 10)

tssMatType <- "normalizedmatrix"
tssMatSuffix <- "normalizedmatrix_3kbATG3kb"
tssMatDim <- c(300, 1, 300, 10)

showExpressionHeatmap <- FALSE

## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF
bsGnom <- BSgenome.Anidulans.FGSCA4.AspGD


## colors
colList <- list()


##################################################################################
sampleList <- suppressMessages(readr::read_tsv(file = file_plotSamples, comment = "#"))

tempSInfo <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = sampleList$sampleId,
  dataPath = TF_dataPath, profileMatrixSuffix = regionMatSuffix
)

tfIds <- tempSInfo$sampleId[which(
  tempSInfo$IP_tag %in% c("HA", "MYC", "TAP") & tempSInfo$TF != "untagged"
)]
kersIds <- tempSInfo$sampleId[which(
  tempSInfo$IP_tag %in% c("HA", "MYC", "TAP") &
    tempSInfo$TF != "untagged" & tempSInfo$complex == "kdmB"
)]
picIds <- tempSInfo$sampleId[which(
  tempSInfo$IP_tag %in% c("HA", "MYC", "TAP") & tempSInfo$complex == "PIC"
)]
inputIds <- tempSInfo$sampleId[which(! tempSInfo$IP_tag %in% c("polII", "HIST") & tempSInfo$TF == "untagged")]
polII_ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "polII")]


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


regionProfileData <- dplyr::bind_rows(polIIData)
tssProfileData <- dplyr::bind_rows(tfData, inputData)

exptData <- dplyr::bind_rows(tfData, inputData)

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

## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "geneId", "score", "strand"))

kmClust <- dplyr::left_join(
  x = suppressMessages(readr::read_tsv(file = tfData$clusterFile[1])),
  y = geneSet, by = c("geneId" = "geneId")
)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, keytype = "GID",
                                  columns = c("GENE_NAME", "DESCRIPTION"))

geneInfo <- dplyr::left_join(x = kmClust, y = geneDesc, by = c("geneId" = "GID"))

glimpse(geneInfo)

expressionData <- get_TF_binding_data(
  exptInfo = dplyr::filter(.data = tfData),
  genesDf = geneInfo
)

expressionData <- get_polII_expressions(
  exptInfo = polIIData, genesDf = expressionData
) %>% 
  dplyr::mutate_at(
    .vars = vars(unname(tfCols$hasPeak)),
    .funs = ~replace_na(data = ., replace = FALSE)
  )

# glimpse(expressionData)
readr::write_tsv(x = expressionData, file = paste(outPrefix, ".data.tab", sep = ""))
##################################################################################

expressionData <- dplyr::mutate(
  .data = expressionData,
  kersPeak = dplyr::if_else(
    condition = dplyr::if_all(
      .cols = all_of(unname(tfCols$hasPeak[kersIds])),
      .fns = ~ . == TRUE
    ),
    true = TRUE, false = FALSE
  )
) %>% 
  dplyr::mutate_at(
    .vars = vars(c(kersPeak, starts_with("hasPeak."))),
    .funs = ~ factor(., levels = c(TRUE, FALSE))
  )

m = make_comb_mat(
  dplyr::select(expressionData, geneId, kersPeak, tfCols$hasPeak[picIds]) %>% 
    dplyr::filter_at(.vars = vars(!geneId), .vars_predicate = any_vars(. == TRUE))
)
UpSet(m)

motifSeqs <- NULL

for (rowId in 1:nrow(tfData)) {
  
  summitSeq500 <- get_peak_summit_seq(
    file = tfData$peakFile[rowId], sampleId = tfData$sampleId[rowId],
    genome = bsGnom, length = 500, column_name_prefix = FALSE) %>% 
    dplyr::rename(
      summitSeqRegion_500 = summitSeqRegion,
      summitSeq_500 = summitSeq
    )
  
  summitSeq200 <- get_peak_summit_seq(
    file = tfData$peakFile[rowId], sampleId = tfData$sampleId[rowId],
    genome = bsGnom, length = 200, column_name_prefix = FALSE) %>% 
    dplyr::rename(
      summitSeqRegion_200 = summitSeqRegion,
      summitSeq_200 = summitSeq
    )
  
  summitSeq <- dplyr::left_join(
    x = summitSeq500, y = summitSeq200, by = "name"
  ) %>% 
    dplyr::rename(peakId = name) %>% 
    dplyr::mutate(sampleId = tfData$sampleId[rowId])
  
  motifSeqs <- dplyr::bind_rows(motifSeqs, summitSeq)
}

## extract summit sequences 
commonSummitSeqs <- dplyr::filter_at(
  .tbl = expressionData, .vars = vars(tfCols$hasPeak), .vars_predicate = all_vars(. == TRUE)
) %>% 
  dplyr::select(peakId = peakId.An_kdmB_20h_HA_1) %>% 
  dplyr::distinct() %>% 
  dplyr::left_join(y = motifSeqs, by = "peakId")

seqSet500 <- NULL
seqSet200 <- NULL
seqSet500 <- DNAStringSet(x = commonSummitSeqs$summitSeq_500)
names(seqSet500) <- commonSummitSeqs$peakId
writeXStringSet(
  x = seqSet500, filepath = paste(outDir, "/KERS_PIC.summitSeqs.500bp.fasta", sep = "")
)

seqSet200 <- DNAStringSet(x = commonSummitSeqs$summitSeq_200)
names(seqSet200) <- commonSummitSeqs$peakId
writeXStringSet(
  x = seqSet200, filepath = paste(outDir, "/KERS_PIC.summitSeqs.200bp.fasta", sep = "")
)


## PIC specific peak summit seq
picSpecificSummitSeqs <- dplyr::filter(.data = expressionData, !kersPeak) %>% 
  dplyr::filter_at(
    .vars = vars(tfCols$hasPeak[picIds]),
    .vars_predicate = all_vars(. == TRUE)
  ) %>% 
  dplyr::select(peakId = peakId.veA_wt_TBP_15h_mix22_2) %>% 
  dplyr::distinct() %>% 
  dplyr::left_join(y = motifSeqs, by = "peakId")

seqSet500 <- NULL
seqSet200 <- NULL
seqSet500 <- DNAStringSet(x = picSpecificSummitSeqs$summitSeq_500)
names(seqSet500) <- picSpecificSummitSeqs$peakId
writeXStringSet(
  x = seqSet500, filepath = paste(outDir, "/PIC_specific.summitSeqs.500bp.fasta", sep = "")
)

seqSet200 <- DNAStringSet(x = picSpecificSummitSeqs$summitSeq_200)
names(seqSet200) <- picSpecificSummitSeqs$peakId
writeXStringSet(
  x = seqSet200, filepath = paste(outDir, "/PIC_specific.summitSeqs.200bp.fasta", sep = "")
)

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
      colorRamp2(breaks = quantile(tfMeanProfile, c(0.40, 0.99), na.rm = T),
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
polIIMeanColor <- colorRamp2(quantile(polIIMeanProfile, c(0.01, 0.5, 0.99), na.rm = T), c("blue", "white", "red"))
polIIColorList <- sapply(X = polIIData$sampleId, FUN = function(x){return(polIIMeanColor)})

## histone colors
histMeanProfile <- NULL
histColorList <- NULL
# if(nrow(histData) == 1){
#   histMeanProfile <- matList[[histIds]]
# } else{
#   histMeanProfile <- getSignalsFromList(lt = matList[histIds])
# }
# quantile(histMeanProfile, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# # histMeanColor <- colorRamp2(quantile(histMeanProfile, c(0.30, 0.995), na.rm = T), c("black", "yellow"))
# histColorList <- sapply(
#   X = c(histIds, histH3Ids),
#   FUN = function(x){
#     return(colorRamp2(breaks = quantile(histMeanProfile, c(0.60, 0.90), na.rm = T),
#                       colors =  unlist(strsplit(x = exptDataList[[x]]$color, split = ","))))
#   }
# )

colorList <- unlist(list(tfColorList, polIIColorList, histColorList))
ylimList <- list()
# ylimList <- sapply(c(histIds, histH3Ids), function(x){return(c(5, 13))}, simplify = FALSE)
# ylimList <- append(x = ylimList,
#                    values = sapply(c(tfIds, inputIds), function(x){return(c(0, 25))}, simplify = FALSE))


##################################################################################
# plot genes which has TF peak in any of the TF samples
kersPeaksDf <- filter_at(
  .tbl = expressionData,
  .vars = unname(tfCols$hasPeak[kersIds]),
  .vars_predicate = all_vars(. == "TRUE")
) %>% 
  dplyr::mutate(group = "peaks")

nGenes <- nrow(kersPeaksDf)

## genes which do not have peak and no polII signal
noPeakPolIIDf <-dplyr::filter_at(
  .tbl = expressionData,
  .vars = unname(c(tfCols$hasPeak, polIICols$is_expressed)),
  .vars_predicate = all_vars(. == FALSE)
) %>% 
  dplyr::slice_min(
    order_by = !!sym(unname(polIICols$exp)), n = nGenes, with_ties = FALSE
  ) %>%
  # dplyr::slice_sample(n = nGenes) %>% 
  dplyr::mutate(group = "control")

plotData <- dplyr::bind_rows(kersPeaksDf, noPeakPolIIDf) %>% 
  dplyr::mutate(
    group = forcats::fct_relevel(.f = group, "peaks", "control")
  )

regionProfiles_peak <- multi_profile_plots(
  exptInfo = regionProfileData,
  genesToPlot = plotData$geneId,
  matSource = regionMatType,
  matBins = regionMatDim,
  drawClusterAn = FALSE,
  profileColors = colorList,
  expressionColor = NULL,
  plotExpression = showExpressionHeatmap,
  column_title_gp = gpar(fontsize = 12),
  ylimFraction = ylimList
)

tssProfiles_peak <- multi_profile_plots(
  exptInfo = tssProfileData,
  genesToPlot = plotData$geneId,
  targetType = "point",
  targetName = "ATG",
  matSource = tssMatType,
  matBins = tssMatDim,
  clusters = dplyr::select(plotData, geneId, cluster = group),
  drawClusterAn = TRUE,
  profileColors = colorList,
  column_title_gp = gpar(fontsize = 12),
  ylimFraction = ylimList
)



## gene length annotation
anGl_peaks <- gene_length_heatmap_annotation(
  bedFile = file_genes, genes = plotData$geneId,
  # pointSize = unit(4, "mm"),
  axis_param = list(at = c(2000, 4000), labels = c("2kb", "> 4kb"))
)

peaks_htlist <- NULL
peaks_htlist <- peaks_htlist + anGl_peaks$an
# peaks_htlist <- peaks_htlist + regionProfiles_peak$heatmapList[, "clusterAn"]
peaks_htlist <- peaks_htlist + tssProfiles_peak$heatmapList
peaks_htlist <- peaks_htlist + regionProfiles_peak$heatmapList

## make sure that the order of genes in the heatmap list and in the dataframe is same
if(all(rownames(peaks_htlist@ht_list[[ tfData$profileName[1] ]]@matrix) == plotData$geneId)){
  
  rowOrd_peaks <- order(plotData[[ tfCols$peakDist[[1]] ]], decreasing = TRUE)
  
}

pdfWd <- 2 + 
  (length(peaks_htlist) * 2) +
  (length(polII_ids) * 0.25 * showExpressionHeatmap)

# wd <- 500 + (nrow(exptData) * 700) + (length(polII_ids) * 500 * showExpressionHeatmap)
title_peak= paste(comparisonName, ": KERS, TBP & TFIIB binding (n = ", nGenes, ")", sep = "")


pdf(file = paste(outPrefix, ".common_peaks.profiles.pdf", sep = ""), width = pdfWd, height = 13)

draw(
  peaks_htlist,
  main_heatmap = exptData$profileName[1],
  # annotation_legend_list = list(profile1$legend),
  column_title = title_peak,
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_sub_title_side = "left",
  heatmap_legend_side = "bottom",
  gap = unit(7, "mm"),
  # row_order = rowOrd_peaks,
  # split = rep(1, nrow(plotData)),
  padding = unit(rep(0.5, times = 4), "cm")
)


dev.off()




##################################################################################
## KERS+PIC mean profile for showing promoter of genes bound by KERS
meanProfile <- geneset_average_profile(
  exptInfo = tfData, genes = kersPeaksDf$geneId, profileMats = matList[tfData$sampleId]
)


lineData <- dplyr::left_join(
  x = meanProfile, y = dplyr::select(tfData, sampleId, TF, complex, color),
  by = "sampleId")

tfLineColors <- dplyr::select(tfData, TF, color) %>%
  tidyr::separate(col = color, sep = ",", into = c(NA, "color")) %>% 
  tibble::deframe()
tfLineColors["TBP"] <- "#ff7f00"
tfLineColors["TFIIB"] <- "black"

pt_meanLine <- ggplot2::ggplot(data = lineData) +
  geom_line(
    mapping = aes(x = bin, y = mean, color = TF, linetype = TF), 
    alpha = 0.7, size = 1.3
  ) +
  scale_color_manual(
    values = tfLineColors
  ) +
  scale_linetype_manual(
    values = setNames(c(1,1,1,1,2,2), names(tfLineColors))
  ) +
  scale_x_continuous(
    breaks = c(1, 100, 200, 300, 400, 500, 600),
    labels = c("-3kb", "-2kb", "-1kb", "ATG", "1kb", "2kb", "3kb"),
    limits = c(100, 500)
  ) +
  labs(
    title = "Mean binding signal of KERS, TBP & TFIIB",
    y = "Normalized coverage"
  )+
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 15), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 18),
    panel.grid = element_blank(),
    legend.position = c(0.9, 0.9),
    legend.key.width = unit(2,"cm"),
    legend.text = element_text(size = 15),
    # legend.key.size = unit(1, "cm"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.justification = c(1,1)
  )  


ggsave(
  plot = pt_meanLine, filename = paste(outPrefix, ".meanProfile.pdf", sep = ""),
  width = 10, height = 8
)


## median
pt_medianLine <- ggplot2::ggplot(data = lineData) +
  geom_line(
    mapping = aes(x = bin, y = median, color = TF, linetype = TF), 
    alpha = 0.7, size = 1.3
  ) +
  scale_color_manual(
    values = tfLineColors
  ) +
  scale_linetype_manual(
    values = setNames(c(1,1,1,1,2,2), names(tfLineColors))
  ) +
  scale_x_continuous(
    breaks = c(1, 100, 200, 300, 400, 500, 600),
    labels = c("-3kb", "-2kb", "-1kb", "ATG", "1kb", "2kb", "3kb"),
    limits = c(100, 500)
  ) +
  labs(
    title = "Mean binding signal of KERS, TBP & TFIIB",
    y = "Normalized coverage"
  )+
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 15), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 18),
    panel.grid = element_blank(),
    legend.position = c(0.9, 0.9),
    legend.key.width = unit(2,"cm"),
    legend.text = element_text(size = 15),
    # legend.key.size = unit(1, "cm"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.justification = c(1,1)
  )  


ggsave(
  plot = pt_medianLine, filename = paste(outPrefix, ".medianProfile.pdf", sep = ""),
  width = 10, height = 8
)


##################################################################################


hasPeakDf <- dplyr::select(expressionData, geneId, kersPeak, unname(tfCols$hasPeak[picIds])) %>% 
  dplyr::filter_at(.vars = vars(!geneId), .vars_predicate = any_vars(. == TRUE)) %>% 
  dplyr::group_by_at(.vars = vars(!geneId)) %>% 
  dplyr::mutate(group = cur_group_id()) %>% 
  dplyr::ungroup()


## topGO enrichment
goEnrich <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(!geneId)) %>%
  do(
    topGO_enrichment(
      genes = .$geneId, orgdb = orgDb, inKeytype = "GID",
      type = "BP", goNodeSize = 5, genenameKeytype = "GENE_NAME"
    )
  )


readr::write_tsv(x = goEnrich, file = paste(outPrefix, ".peakGroups.topGO.tab", sep = ""))

picSpecificGo <- dplyr::filter(.data = goEnrich, kersPeak == FALSE) %>% 
  dplyr::filter_at(.vars = vars(tfCols$hasPeak[picIds]), .vars_predicate = all_vars(. == TRUE)) %>% 
  dplyr::arrange(desc(log10_pval)) %>% 
  dplyr::slice_head(n = 10) %>% 
  dplyr::mutate(
    Term = stringr::str_replace(string = Term, pattern = "\\(.*\\)", replacement = ""),
    Term = forcats::as_factor(Term)
    )

pt_goBar <- ggplot2::ggplot(data = picSpecificGo) +
  geom_bar(
    mapping = aes(x = log10_pval, y = Term, fill = richness),
    stat = "identity"
  ) +
  scale_fill_viridis(name = "enrichment") +
  scale_x_continuous(expand = expansion(add = c(0, 1))) +
  labs(
    title = "Enriched GO BPs in genes bound by TBP+TFIIB and NOT by KERS",
    x = "log10(p-value)", y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(hjust = 1, size = 14, face = "bold")
  )

ggsave(
  filename = paste(outPrefix, ".TBP_TFIIB_specific_GO.pdf", sep = ""),
  plot = pt_goBar, width = 8, height = 5
)

## pathway enrichment
keggEnr <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(!geneId)) %>%
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















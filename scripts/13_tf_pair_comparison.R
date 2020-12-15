library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(scales)
library(ggplot2)
library(summarytools)


## 1) compares the binding of two TF pairs
## 2) perform GO enrichment, GO grouping, KEGG enrichment for different groups
## 3) generate TF binding change plot for TF pairs w.r.t. polII fold change


rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/GO_enrichment/topGO_functions.R")


##################################################################################
## main configuration
comparisonName <- "kdmB_48h_vs_20h"
outPrefix <- here::here("kdmB_analysis", comparisonName, comparisonName)

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

tf1 <- tfDiffPairs[[mainTfPair]]$samples[1]
tf2 <- tfDiffPairs[[mainTfPair]]$samples[2]

otherHist <- c("An_H3_20h_HIST_1", "An_H3_48h_HIST_1")


## polII signal fold change pairs
polIIDiffPairs <- list(
  p1 = list(
    name = "untagged_48h_vs_untagged_20h_polII",
    title = "polII log2(untagged_48h \n vs untagged_20h)",
    samples = c("An_untagged_20h_polII_1", "An_untagged_48h_polII_1")
  ),
  p2 = list(
    name = "kdmB_del_48h_vs_kdmB_del_20h_polII",
    title = "polII log2(kdmB_del_48h \n vs kdmB_del_20h)",
    samples =  c("An_kdmB_del_20h_polII_1", "An_kdmB_del_48h_polII_1")
  ),
  p3 = list(
    name = "kdmB_del_20h_vs_untagged_20h_polII",
    title = "polII log2(kdmB_del_20h \n vs untagged_20h_polII)",
    samples = c("An_untagged_20h_polII_1", "An_kdmB_del_20h_polII_1")
  ),
  p4 = list(
    name = "kdmB_del_48h_vs_untagged_48h_polII",
    title = "polII log2(kdmB_del_48h \n vs untagged_48h_polII)",
    samples = c("An_untagged_48h_polII_1", "An_kdmB_del_48h_polII_1")
  )
)

mainPolIIPair <- "p1"

polII1 <- polIIDiffPairs[[mainPolIIPair]]$samples[1]
polII2 <- polIIDiffPairs[[mainPolIIPair]]$samples[2]

# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
matrixType <- "normalizedmatrix_5kb"
matrixDim = c(500, 200, 100, 10)

clusterStorePath <- paste(outPrefix, "_profile.kmeans.clusters.txt", sep = "")

## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")

orgDb <- org.Anidulans.FGSCA4.eg.db

showExpressionHeatmap <- FALSE

outPrefix_tfPolII <- paste(outPrefix, "_tfPolII", sep = "")
outPrefix_tfSpecific <- paste(outPrefix, "_specific_binding", sep = "")
outPrefix_peakExp <- paste(outPrefix, "_pkExpGenes", sep = "")

anLables <- list()
anLables[["is_SM_gene"]] <- "SM gene"
anLables[["is_TF"]] <- "Transcription Factor"
anLables[["gene_length"]] <- "Gene Length"


##################################################################################

## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- suppressMessages(
  AnnotationDbi::select(x = orgDb, keys = geneSet$gene, columns = "DESCRIPTION", keytype = "GID")
)

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))

## gene information annotations: cluster and TF and polII expression values
geneInfo <- add_gene_info(file = file_geneInfo, clusterDf = geneSet)

head(geneInfo)

##################################################################################

## read the experiment sample details and select only those which are to be plotted

tfData <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = unique(purrr::map(tfDiffPairs, "samples") %>% unlist() %>% unname()),
  dataPath = TF_dataPath,
  matrixSource = matrixType)

polIIData <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = unique(purrr::map(polIIDiffPairs, "samples") %>% unlist() %>% unname()),
  dataPath = polII_dataPath,
  matrixSource = matrixType)

histData <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = otherHist,
  dataPath = hist_dataPath,
  matrixSource = matrixType)

exptData <- dplyr::bind_rows(tfData, polIIData, histData)
exptDataList <- purrr::transpose(exptData)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

polII_ids <- exptData$sampleId[which(exptData$IP_tag == "polII")]
tfIds <- exptData$sampleId[which(exptData$IP_tag %in% c("HA", "MYC", "TAP") & exptData$TF != "untagged")]


polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)


tfCols <- sapply(c("hasPeak", "pval", "peakType", "tesPeakType", "peakCoverage", "peakDist", "summitDist", "upstreamExpr", "peakExpr", "relativeDist"),
                 FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
                 simplify = F, USE.NAMES = T)


expressionData <- get_polII_expressions(exptInfo = exptData,
                                        genesDf = geneInfo)

## add fold change columns
for (i in names(polIIDiffPairs)) {
  expressionData <- get_fold_change(
    df = expressionData,
    nmt = polIIDiffPairs[[i]]$samples[2],
    dmt = polIIDiffPairs[[i]]$samples[1],
    newCol = polIIDiffPairs[[i]]$name,
    isExpressedCols = polIICols$is_expressed)
}


expressionData <- get_TF_binding_data(exptInfo = tfData,
                                      genesDf = expressionData)


expressionData %>% 
  dplyr::select(gene, starts_with("hasPeak")) %>% 
  dplyr::group_by_at(.vars = vars(starts_with("hasPeak"))) %>% 
  dplyr::summarise(n = n())


expressionData$group <- dplyr::group_by_at(expressionData, unname(tfCols$hasPeak[c(tf1, tf2)])) %>%
  dplyr::group_indices()

rownames(expressionData) = expressionData$gene

fwrite(x = expressionData, file = paste(outPrefix, "_data.tab", sep = ""),
       sep = "\t", col.names = T, quote = F, na = "NA")

view(dfSummary(expressionData))

## peak data
hasPeakDf <- dplyr::filter_at(.tbl = expressionData,
                              .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)])),
                              .vars_predicate = any_vars(. == "TRUE"))

rownames(hasPeakDf) <- hasPeakDf$gene


dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)]), group)) %>%
  dplyr::summarise(n = n())

## check the MA plot
plot_MA_gg(df = expressionData, s1 = polII2, s2 = polII1, title = "MA plot", colorCol = "group")
# plot_MA_gg(df = expressionData, s1 = polII1, s2 = polII2, title = "MA plot", colorCol = "group")
plot_MA_gg(df = expressionData, s1 = polIIDiffPairs$p2$samples[1],
           s2 = polIIDiffPairs$p2$samples[2], title = "MA plot", colorCol = "group")
# plot_MA_gg(df = expressionData, title = "MA plot", colorCol = "group",
#            s1 = polIIDiffPairs$p2$samples[1], s2 = polIIDiffPairs$p2$samples[2])


plot_scatter(df = expressionData, s1 = polII1, s2 = polII2, title = "Scatter plot", colorCol = "group")
plot_scatter(df = expressionData, s1 = polII2, s2 = polII1, title = "Scatter plot", colorCol = "group")

peakCovDf <- tidyr::gather(data = hasPeakDf, key = "sample", value = "coverage", starts_with("peakCoverage.")) %>% 
  dplyr::select(gene, sample, coverage) %>% 
  dplyr::mutate(
    sample = gsub(pattern = "peakCoverage.", replacement = "", x = sample, fixed = T),
    tf = gsub(pattern = "An_(\\w+)_(20h|48h)_HA_1", replacement = "\\1", x = sample, perl = T),
    time = gsub(pattern = "An_(\\w+)_(20h|48h)_HA_1", replacement = "\\2", x = sample, perl = T)) %>% 
  as.tibble()

ggplot(data = peakCovDf,
       mapping = aes(x = coverage, group = sample, color = tf, linetype = time)) +
  # geom_density(size = 1) +
  geom_line(stat = "density", size = 1) +
  coord_cartesian(xlim = c(0, 100)) +
  facet_wrap(tf ~ ., nrow = 2, ncol = 2)


##################################################################################

## join the profile matrices and do clustering
# mergedKm = merged_profile_matrix_cluster(name = comparisonName,
#                                          exptInfo = exptData,
#                                          genes = geneSet$gene,
#                                          clusterStorePath = clusterStorePath,
#                                          k = 12)

# clusterInfo <- fread(file = clusterStorePath, sep = "\t", header = T, stringsAsFactors = F)
# 
# expressionData <- dplyr::left_join(x = expressionData, y = clusterInfo, by = c("gene" = "gene"))




##################################################################################
## topGO enrichment
goEnrich <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)]), group)) %>%
  do(topGO_enrichment(goMapFile = file_topGoMap, genesOfInterest = .$gene, goNodeSize = 5))


fwrite(x = goEnrich,
       file = paste(outPrefix, "_GO.tab", sep = ""), sep = "\t", col.names = T, quote = F)


## clusterProfiler groupGO assignment
grpGo <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)]), group)) %>% 
  do(clusterProfiler_groupGO(genes = .$gene, org = orgDb, goLevel = 3, type = "BP", keyType ="GID"))

fwrite(x = grpGo,
       file = paste(outPrefix, "_groupGO.tab", sep = ""), sep = "\t", col.names = T, quote = F)


## pathway enrichment
keggEnr <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(unname(tfCols$hasPeak[c(tf1, tf2)]), group)) %>% 
  do(keggprofile_enrichment(genes = .$gene, orgdb = orgDb, keytype = "GID", keggOrg = "ani", pvalCut = 0.05))


fwrite(x = keggEnr,
       file = paste(outPrefix, "_KEGG.tab", sep = ""), sep = "\t", col.names = T, quote = F)


##################################################################################
# ## GO enrichment and KEGG testing
# gl <- hasPeakDf$gene[ which( !hasPeakDf[[tfCols$hasPeak[tf1]]] & hasPeakDf[[tfCols$hasPeak[tf2]]] ) ]
# 
# 
# ## KEGG pathway enrichment using enrichKEGG
# keggIds <- AnnotationDbi::select(x = orgDb, keys = gl, columns = c("GID", "KEGG_ID"), keytype = "GID") %>% 
#   dplyr::filter(!is.na(KEGG_ID))
# 
# kk <- clusterProfiler::enrichKEGG(gene = keggIds$KEGG_ID, organism = 'ani', keyType = "kegg",
#                                   pvalueCutoff = 0.05, pAdjustMethod = "none", qvalueCutoff = 1,
#                                   minGSSize = 1)
# 
# 
# kkRes <- dplyr::filter(kk@result, pvalue <= 0.05)
# 
# ## using KEGGprofile package
# kp <- KEGGprofile::find_enriched_pathway(gene = keggIds$KEGG_ID, species = 'ani',
#                                          returned_pvalue = 0.05, returned_adjpvalue = 1,
#                                          returned_genenumber = 1, download_latest = TRUE)
# 
# kp$stastic

##################################################################################
## clusterProfiler enrichment: testing

# ## groupGO
# grpGo <- clusterProfiler::groupGO(gene = gl,
#                                   OrgDb = orgDb,
#                                   keyType = "GID", ont = "BP", level = 3)
# 
# grpGoDf <- dplyr::filter(grpGo@result, Count > 0)
# 
# grpGoDf[, 1:4]
# 
# ## enricher
# t2g <- AnnotationDbi::select(x = orgDb, keys = keys(orgDb), columns = c("GO", "GID","ONTOLOGY"), keytype = "GID") %>% 
#   dplyr::filter(ONTOLOGY == "BP") %>% 
#   dplyr::select(GO, GID)
# 
# t2n <- AnnotationDbi::select(x = GO.db, keys = unique(t2g$GO), columns = c("GOID", "TERM"), keytype = "GOID")
# 
# enrGo <- clusterProfiler::enricher(gene = gl, TERM2GENE = t2g, TERM2NAME = t2n)
# enrGoDf <- dplyr::filter(enrGo@result, pvalue <= 0.05)
# 
# enrGoDf[, 1:7]
# 
# 
# ## enricher
# enrGo2 <- clusterProfiler::enrichGO(gene = gl, OrgDb = orgDb, keyType = "GID", ont = "BP",
#                                     pvalueCutoff = 0.05)
# 
# enrGo2Df <- dplyr::filter(enrGo2@result, pvalue <= 0.05)
# enrGo2Df[, 1:7]
# 
# write.table(enrGo2Df, "clipboard", sep = "\t", quote = F)
# 
# 
# ## formula based enrichment
# keytypes(orgDb)
# reformulate(termlabels = c(tfCols$hasPeak, "group"), response = "gene")
# 
# formulaRes <- clusterProfiler::compareCluster(
#   geneClusters = gene ~ group,
#   fun = "enrichGO",
#   data = hasPeakDf,
#   OrgDb = orgDb,
#   keyType = "GID",
#   ont = "BP", minGSSize = 2)

# formSimple <- clusterProfiler::simplify(formulaRes, cutoff = 0.7)
# 
# plt <- dotplot(formSimple, showCategory = NULL) +
#   scale_y_discrete(labels = wrap_format(60)) +
#   ggtitle(paste(titleName, ": GO enrichment for DEGs")) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.8, size = 16, face = "bold"),
#         axis.text.x = element_text(size = 14),
#         axis.text.y = element_text(size = 14),
#         axis.title = element_text(face = "bold"),
#         panel.grid = element_blank(),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14, face = "bold"))
# 

##################################################################################

## polII signal matrix
polIIMat <- data.matrix(log2(expressionData[, polII_ids] + 1))

quantile(polIIMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.997, 0.999, 0.9999, 1), na.rm = T)


polII_color <- colorRamp2(breaks = c(0, quantile(polIIMat, c(0.5, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99, 0.995, 0.999))),
                          colors = c("white", RColorBrewer::brewer.pal(n = 9, name = "RdPu")))


## polII signal fold change matrix
lfcMat <- as.matrix(expressionData[, purrr::map_chr(polIIDiffPairs, "name"), drop = FALSE])

lfc_color <- colorRamp2(breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2),
                        colors = RColorBrewer::brewer.pal(n = 7, name = "PuOr"))

##################################################################################

# ## peak pval heatmap
# quantile(hasPeakDf[, tfCols$pval], c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
# 
# peakPvalCol <- colorRamp2(breaks = c(0, 1, quantile(hasPeakDf[, tfCols$pval], c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 0.95), na.rm = T)),
#                           colors = c("grey", "white", RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")))
# 
# 
# macs2Ht <- signal_heatmap(log2_matrix = hasPeakDf[, tfCols$pval],
#                           htName = "macs2_peak",
#                           col_title = tfCols$pval,
#                           legend_title = "log10(macs2_peak)",
#                           color = peakPvalCol,
#                           htSplit = hasPeakDf["group"],
#                           cluster_rows = TRUE, cluster_columns = FALSE)
# 
# 
# polIIMat_peaks <- polIIMat[hasPeakDf$gene , ]
# 
# polIIht_peaks <- signal_heatmap(log2_matrix = polIIMat,
#                                 htName = "polII_exp",
#                                 col_title = polII_ids,
#                                 legend_title = "log2(polII_singal)",
#                                 color = polII_color)

##################################################################################
## colors for profile matrix
matList <- import_profiles(exptInfo = exptData,
                               geneList = geneInfo$gene,
                               source = matrixType,
                               up = matrixDim[1], target = matrixDim[2], down = matrixDim[3])



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
    return(colorRamp2(breaks = quantile(tfMeanProfile, c(0.50, 0.995), na.rm = T),
                      colors = c("white", exptDataList[[x]]$color)))
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
      colorRamp2(quantile(matList[[x]], c(0.50, 0.99), na.rm = T), c("white", exptDataList[[x]]$color))
    )
  }
)


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
    return(colorRamp2(breaks = quantile(histMeanProfile, c(0.20, 0.995), na.rm = T),
                      colors = c("black", exptDataList[[x]]$color)))
  }
)

## LFC(TF2/TF1) matrix
## adding Tf1's median to TF2 matrix and vice a versa. Doing this will remove the background skewness of fold change
## matrix and the median will be around 0
tfLfcMat <- log2((matList[[tf2]] + quantile(matList[[tf1]], 0.5)) / (matList[[tf1]] + quantile(matList[[tf2]], 0.5)))

quantile(log2((matList[[tf2]] + quantile(matList[[tf2]], 0.5)) / (matList[[tf1]] + quantile(matList[[tf1]], 0.5))),
         c(0, 0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 1))

plot(density(log2((matList[[tf2]] + quantile(matList[[tf1]], 0.5)) / (matList[[tf1]] + quantile(matList[[tf2]], 0.5))),
             na.rm = TRUE)
)

lfcProfileCol <- colorRamp2(breaks = c(-2, -1.5, -1, -0.75, 0, 0.75, 1, 1.5, 2),
                            colors = RColorBrewer::brewer.pal(n = 9, name = "PuOr"))

## scalled LFC matrix
tfScalledLfcMat <- scale(x = tfLfcMat, center = TRUE, scale = TRUE)

quantile(tfScalledLfcMat,
         c(0, 0.005, 0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 0.995, 1))

## TF1 scalled matrix
tf1ScalledMat <- scale(x = matList[[tf1]], center = TRUE, scale = TRUE)
quantile(tf1ScalledMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
tf1ScalledColor <- colorRamp2(quantile(tf1ScalledMat, c(0.50, 0.99), na.rm = T), c("white", exptDataList[[tf1]]$color))

## TF2 scalled matrix
tf2ScalledMat <- scale(x = matList[[tf2]], center = TRUE, scale = TRUE)
quantile(tf2ScalledMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
tf2ScalledColor <- colorRamp2(quantile(tf2ScalledMat, c(0.50, 0.99), na.rm = T), c("white", exptDataList[[tf2]]$color))

## Difference between TF2 and TF1 scalled matrix
scalledTfDiffMat <- tf2ScalledMat - tf1ScalledMat
plot(density(scalledTfDiffMat))
quantile(scalledTfDiffMat, c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)


scalledTfDiffColor <- colorRamp2(breaks = c(-3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3),
                                 colors = RColorBrewer::brewer.pal(n = 11, name = "PuOr"))


##################################################################################
## genes which have peak in both TFs and polII signal in either one or both TFs

peakExpDf <- dplyr::filter_at(.tbl = expressionData,
                              .vars = unname(tfCols$hasPeak[c(tf1, tf2)]),
                              .vars_predicate = all_vars(. == TRUE)) %>% 
  dplyr::filter_at(.vars = unname(polIICols$is_expressed[c(polII1, polII2)]),
                   .vars_predicate = any_vars(. == TRUE)) 

## group genes into polII fold change bins
peakExpDf <- dplyr::mutate(peakExpDf,
                           group = case_when(
                             !! as.name(polIIDiffPairs$p1$name) >= 2 ~ "G1: LFC >= 2",
                             !! as.name(polIIDiffPairs$p1$name) >= 1 ~ "G2: 1 <= LFC < 2",
                             !! as.name(polIIDiffPairs$p1$name) >= 0.5 ~ "G3: 0.5 <= LFC < 1",
                             !! as.name(polIIDiffPairs$p1$name) >= 0 ~ "G4: 0 <= LFC < 0.5",
                             !! as.name(polIIDiffPairs$p1$name) <= -2 ~ "G8: -2 >= LFC",
                             !! as.name(polIIDiffPairs$p1$name) <= -1 ~ "G7: -2 < LFC <= -1",
                             !! as.name(polIIDiffPairs$p1$name) <= -0.5 ~ "G6: -1 < LFC <= -0.5",
                             !! as.name(polIIDiffPairs$p1$name) < 0 ~ "G5: -0.5 < LFC < 0",
                             TRUE ~ "0"
                           ))

plot_scatter(df = peakExpDf, s1 = polII1, s2 = tfCols$peakCoverage[tf2], title = "Scatter plot", colorCol = "group")
plot_MA_gg(df = peakExpDf, s1 = tfCols$peakCoverage[tf2], s2 = tfCols$peakCoverage[tf1], colorCol = "group")
plot_scatter(df = peakExpDf, s1 = tfCols$peakCoverage[tf1], s2 = tfCols$peakCoverage[tf2], title = "Scatter plot", colorCol = "group")

dplyr::group_by(peakExpDf, group) %>% 
  dplyr::summarise(n = n()) %>% 
  readr::write_tsv(path = paste(outPrefix_peakExp, "_stats.tab", sep = ""))


peakExp_clusters <- dplyr::select(peakExpDf, gene, group) %>% 
  dplyr::rename(cluster = group)


## profile heatmap: no need to include input
multiProf_peakExp <- multi_profile_plots(
  exptInfo = tfData[tfData$sampleId %in% tfIds, ],
  genesToPlot = peakExpDf$gene,
  matSource = matrixType,
  matBins = matrixDim,
  clusters = peakExp_clusters,
  profileColors = tfWiseColors,
  # ylimFraction = ylimList,
  column_title_gp = gpar(fontsize = 12))


## Scalled TF diff profile
scalledTfDiffProf_peakExp <- profile_heatmap(
  profileMat = scalledTfDiffMat[peakExpDf$gene, ],
  signalName = comparisonName,
  profileColor = scalledTfDiffColor,
  column_title_gp = gpar(fontsize = 12),
  geneGroups = peakExp_clusters,
  ylimFraction = c(-2.3, 1.7))


## polII signal heatmap
polIIMat_peakExp <- polIIMat[peakExpDf$gene, ]

polIIht_peakExp <- signal_heatmap(
  log2_matrix = polIIMat_peakExp,
  htName = "polII_exp",
  col_title = polII_ids,
  legend_title = "log2(polII_singal)",
  color = polII_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE)

## polII fold change heatmap
lfc_peakExp <- lfcMat[peakExpDf$gene, ]

lfcHt_peakExp <- signal_heatmap(
  log2_matrix = lfc_peakExp,
  htName = "polII_lfc",
  col_title = purrr::map_chr(polIIDiffPairs, "title"),
  legend_title = "log2(fold change)",
  color = lfc_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE)

## gene length annotation
gl_peakExp <- gene_length_heatmap_annotation(bedFile = file_genes, genes = peakExpDf$gene)

htlist_peakExp <- gl_peakExp$an +
  lfcHt_peakExp + 
  polIIht_peakExp +
  multiProf_peakExp$heatmapList +
  scalledTfDiffProf_peakExp$heatmap



if( all(rownames(multiProf_peakExp$profileHeatmaps[[1]]$heatmap@matrix) == peakExpDf$gene) ){
  ## row order by polII LFC
  rowOrd_peakExp <- order(peakExpDf[[ polIIDiffPairs$p1$name ]], decreasing = TRUE)
}


title_peakExp <- "TF binding vs polII transcription: genes bound both at 20h and 48h and atleast one timepoint showing polII signal "

pdfWd <- 4 + (length(multiProf_peakExp$profileHeatmaps) * 3) + 3 +
  (length(polII_ids) * 0.5 * showExpressionHeatmap) + 
  (length(polIIDiffPairs) * 0.5)

# draw Heatmap and add the annotation name decoration
pdf(file = paste(outPrefix_peakExp, "_profile.pdf", sep = ""), width = pdfWd, height = 12)

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

## decorate the annotations
row_annotation_axis(an = "gene_length",
                    at = c(0, 2000, 4000),
                    labels = c("0kb", "2kb", ">4kb"),
                    slice = length(unique(peakExp_clusters$cluster)))

dev.off()


## box plots
boxDt <- data.table::melt(
  data = as.data.table(peakExpDf),
  measure.vars = list(tfCols$hasPeak, tfCols$peakCoverage),
  variable.name = "sample",
  value.name = c("hasPeak", "peakCoverage")
)

boxDt$sample <- forcats::fct_recode(
  .f = boxDt$sample,
  structure(.Data = levels(boxDt$sample), names = names(tfCols$hasPeak))
)

boxDt <- dplyr::left_join(boxDt, dplyr::select(tfData, sampleId, TF, timepoint), by = c("sample" = "sampleId")) %>% 
  dplyr::group_by(sample) %>% 
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


pdf(file = paste(outPrefix_peakExp, "_box.pdf", sep = ""), width = 8, height = 12)
pt
dev.off()


##################################################################################
## profile plot with genes showing peak specific and common between the conditions
# tfSpecificDf <- dplyr::filter(hasPeakDf,
#                               !! as.name(tfCols$hasPeak[tf1]) == FALSE | !! as.name(tfCols$hasPeak[tf2]) == FALSE )

tfSpecificDf <- hasPeakDf

rownames(tfSpecificDf) <- tfSpecificDf$gene

clusters_tfSpecific <- dplyr::select(tfSpecificDf, gene, group) %>% 
  dplyr::rename(cluster = group)

ylimList <- sapply(X = tfData$sampleId, FUN = function(x){c(0, 30)}, USE.NAMES = T, simplify = F)

## TF profile plot
multiProf_tfSpecific <- multi_profile_plots(
  exptInfo = tfData[tfData$sampleId %in% tfIds, ],
  genesToPlot = tfSpecificDf$gene,
  matSource = matrixType,
  matBins = matrixDim,
  profileColors = tfMeanColorList,
  clusters = clusters_tfSpecific,
  column_title_gp = gpar(fontsize = 12),
  ylimFraction = ylimList)


## Scalled TF diff profile
scalledTfDiffProf_tfSpecific <- profile_heatmap(
  profileMat = scalledTfDiffMat[tfSpecificDf$gene, ],
  signalName = comparisonName,
  profileColor = scalledTfDiffColor,
  column_title_gp = gpar(fontsize = 12),
  geneGroups = clusters_tfSpecific)

## histone profile plots
histProfiles <- multi_profile_plots(
  exptInfo = histData,
  genesToPlot = tfSpecificDf$gene,
  matSource = matrixType,
  matBins = matrixDim,
  profileColors = histColorList,
  clusters = clusters_tfSpecific,
  column_title_gp = gpar(fontsize = 12),
  drawClusterAn = FALSE)

## polII signal heatmap
polIIMat_tfSpecific <- polIIMat[tfSpecificDf$gene, ]

polIIHt_tfSpecific <- signal_heatmap(
  log2_matrix = polIIMat_tfSpecific,
  htName = "polII_exp",
  col_title = polII_ids,
  legend_title = "log2(polII_singal)",
  color = polII_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE)

## polII signal fold change heatmap
lfc_tfSpecific <- lfcMat[tfSpecificDf$gene, ]

lfc_heatmap <- signal_heatmap(
  log2_matrix = lfc_tfSpecific,
  htName = "polII_lfc",
  col_title = purrr::map_chr(polIIDiffPairs, "title"),
  legend_title = "log2(fold change)",
  color = lfc_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE)


gl_tfSpecific <- gene_length_heatmap_annotation(bedFile = file_genes, genes = tfSpecificDf$gene)

htlist_tfSpecific <- gl_tfSpecific$an + 
  multiProf_tfSpecific$heatmapList + 
  scalledTfDiffProf_tfSpecific$heatmap +
  polIIHt_tfSpecific + 
  lfc_heatmap +
  histProfiles$heatmapList



if( all(rownames(multiProf_tfSpecific$profileHeatmaps[[1]]$heatmap@matrix) == tfSpecificDf$gene) ){
  rowOrd_tfSpecific <- order(tfSpecificDf[[ tfCols$peakDist[tf1] ]], tfSpecificDf[[ tfCols$peakDist[tf2] ]],
                             decreasing = TRUE)
}


title_tfSpecific <- "Differential binding of kdmB at 20h and 48h"

pdfWd <- 2 + (length(multiProf_tfSpecific$profileHeatmaps) * 2) + 3 +
  (length(polII_ids) * 0.5 * showExpressionHeatmap) + 
  (length(polIIDiffPairs) * 0.5) + 
  (length(histProfiles$profileHeatmaps) * 2) + 3

# draw Heatmap and add the annotation name decoration
pdf(file = paste0(outPrefix_tfSpecific, "_profile.pdf", collapse = ""), width = pdfWd, height = 12)
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

## decorate the annotations
row_annotation_axis(an = "gene_length",
                    at = c(0, 2000, 4000),
                    labels = c("0kb", "2kb", ">4kb"),
                    slice = length(unique(tfSpecificDf$group)))

dev.off()


##################################################################################
## genes with either peak detected or with top 10% polII signal 

expressionData %>% 
  dplyr::group_by_at(.vars = vars(unname(c(tfCols$hasPeak[c(tf1, tf2)], polIICols$is_expressed[c(polII1, polII2)])))) %>% 
  dplyr::summarise(n = n())

tfPolIIDf <- dplyr::filter_at(.tbl = expressionData,
                              .vars = unname(c(tfCols$hasPeak[c(tf1, tf2)], polIICols$is_expressed[c(polII1, polII2)])),
                              .vars_predicate = any_vars(. == TRUE))

tfPolIIDf$group <- tfPolIIDf %>% 
  dplyr::group_by_at(.vars = vars(unname(c(tfCols$hasPeak[c(tf1, tf2)], polIICols$is_expressed[c(polII1, polII2)])))) %>% 
  dplyr::group_indices()

tfpolII_clusters <- dplyr::mutate(tfPolIIDf, group = sprintf(fmt = "%02d", group)) %>% 
  dplyr::select(gene, group) %>% 
  dplyr::rename(cluster = group)

tfPolIIDf %>% 
  dplyr::group_by_at(
    .vars = vars(unname(c(tfCols$hasPeak[c(tf1, tf2)], polIICols$is_expressed[c(polII1, polII2)])), "group")) %>% 
  dplyr::summarise(n = n())

## profile heatmap
multiProf_tfPolII <- multi_profile_plots(
  exptInfo = tfData[tfData$sampleId %in% tfIds, ],
  genesToPlot = tfPolIIDf$gene,
  matSource = matrixType,
  matBins = matrixDim,
  profileColors = tfMeanColorList,
  column_title_gp = gpar(fontsize = 12),
  clusters = tfpolII_clusters)


## Scalled TF diff profile
scalledTfDiffProf_tfPolII <- profile_heatmap(
  profileMat = tfLfcMat[tfPolIIDf$gene, ],
  signalName = comparisonName,
  profileColor = scalledTfDiffColor,
  geneGroups = tfpolII_clusters)

## polII signal heatmap
polIIMat_tfPolII <- polIIMat[tfPolIIDf$gene, ]

polIIht_tfPolII <- signal_heatmap(
  log2_matrix = polIIMat_tfPolII,
  htName = "polII_exp",
  col_title = polII_ids,
  legend_title = "log2(polII_singal)",
  color = polII_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE)

## polII fold change heatmap
lfc_tfPolII <- lfcMat[tfPolIIDf$gene, ]

lfcHt_tfPolII <- signal_heatmap(
  log2_matrix = lfc_tfPolII,
  htName = "polII_lfc",
  col_title = purrr::map_chr(polIIDiffPairs, "title"),
  legend_title = "log2(fold change)",
  color = lfc_color,
  column_title_gp = gpar(fontsize = 12),
  cluster_columns = FALSE)

## gene length annotation
gl_tfPolII <- gene_length_heatmap_annotation(bedFile = file_genes, genes = tfPolIIDf$gene)

htlist_tfPolII <- gl_tfPolII$an +
  multiProf_tfPolII$heatmapList +
  scalledTfDiffProf_tfPolII$heatmap +
  polIIht_tfPolII + lfcHt_tfPolII



if( all(rownames(multiProf_tfPolII$profileHeatmaps[[1]]$heatmap@matrix) == tfPolIIDf$gene) ){
  ## order by peak distance
  rowOrd_tfPolII <- order(tfPolIIDf[[ tfCols$peakDist[tf1] ]], tfPolIIDf[[ polIIDiffPairs$p1$name ]],
                          decreasing = TRUE)
}


title_tfPolII <- "Differential binding of TF at 20h and 48h: genes with macs2 peaks or top 10% polII signal"


pdfWd <- 2 + (length(multiProf_tfPolII$profileHeatmaps) * 3.5) + 3 +
  (length(polII_ids) * 0.5 * showExpressionHeatmap) + 
  (length(polIIDiffPairs) * 0.5)

# draw Heatmap and add the annotation name decoration
pdf(file = paste0(outPrefix_tfPolII, "_profile.pdf", collapse = ""), width = pdfWd, height = 12)

draw(htlist_tfPolII,
     main_heatmap = tfData$profileName[1],
     # annotation_legend_list = list(profile1$legend),
     column_title = title_tfPolII,
     column_title_gp = gpar(fontsize = 14, fontface = "bold"),
     row_sub_title_side = "left",
     gap = unit(6, "mm"),
     row_order = rowOrd_tfPolII,
     padding = unit(rep(0.5, times = 4), "cm")
)

## decorate the annotations
row_annotation_axis(an = "gene_length",
                    at = c(0, 2000, 4000),
                    labels = c("0kb", "2kb", ">4kb"),
                    slice = length(unique(tfpolII_clusters$cluster)))

dev.off()


##################################################################################



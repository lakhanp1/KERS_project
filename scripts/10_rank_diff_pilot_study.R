library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(scales)
library(ggplot2)
library(corrplot)

## 1) compares the binding of two TF replicate pairs
## 2) find the correlation between the replicate peaks
## 3) a) peak coverage correlation for master peak set formed by merging replicate peaks
##    b) peak coverage correlation using original replicate peakset regions, wherever 
##       peaks are not found, master region coverage is used
##    c) correlation of peak coverage (from b) quartile ranks formed by grouping into 100 quartiles 


rm(list = ls())

source(file = "E:/Chris_UM/Codes/GO_enrichment/topGO_functions.R")


##################################################################################
## main configuration
comparisonName <- "replicate_cor"
outPrefix <- here::here("kdmB_analysis", "rank_diff_pilot_study", comparisonName)

file_tfReps <- here::here("kdmB_analysis", "rank_diff_pilot_study", "tf_replicates.tab")

tfReps <- readr::read_tsv(file = file_tfReps)


# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
matrixType <- "normalizedmatrix_5kb"
matrixDim = c(500, 200, 100, 10)

# clusterStorePath <- paste(outPrefix, "_profile.kmeans.clusters.txt", sep = "")

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

# ## genes to read
# geneSet <- data.table::fread(file = file_genes, header = F,
#                              col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
#   dplyr::mutate(length = end - start)
# 
# geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, columns = "DESCRIPTION", keytype = "GID")
# 
# geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))
# 
# ## gene information annotations: cluster and TF and polII expression values
# geneInfo <- add_gene_info(file = file_geneInfo, clusterDf = geneSet)
# 
# head(geneInfo)

##################################################################################


tfData <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = unique(c(tfReps$rep1, tfReps$rep2)),
                                 dataPath = TF_dataPath,
                                 matrixSource = matrixType)

polIIData <- NULL
histData <- NULL

exptData <- dplyr::bind_rows(tfData, polIIData, histData)

tfIds <- exptData$sampleId[which(exptData$IP_tag %in% c("HA", "MYC", "TAP") & exptData$TF != "untagged")]


tfCols <- purrr::map(
  .x = tfIds,
  .f = function(x){
    list(
      id = x,
      cov = paste("cov.", x, sep = ""),
      qtBin = paste("qtbin.", x, sep = ""),
      op = paste("op.", x, sep = "")
    )
  }) %>%
  purrr::set_names(tfIds)


exptDataList <- purrr::transpose(exptData) %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))


skipChr <- "mito_A_nidulans_FGSC_A4"

corList <- list()
i <- 5

for(i in 1:nrow(tfReps)){
  tf1 <- tfReps$rep1[i]
  tf2 <- tfReps$rep2[i]
  sampleId <- tfReps$sample[i]
  
  ## import narrowPeak files
  peaks1 <- rtracklayer::import(con = exptDataList[[tf1]]$narrowpeakFile, format = "narrowPeak")
  peaks2 <- rtracklayer::import(con = exptDataList[[tf2]]$narrowpeakFile, format = "narrowPeak")
  
  
  if(! is.null(skipChr)){
    peaks1 <- peaks1[! seqnames(peaks1) %in% skipChr, ]
    peaks2 <- peaks2[! seqnames(peaks2) %in% skipChr, ]
  }
  
  ## peak region coverage for each peakset
  peaks1 <- region_coverage(regions = peaks1, bwFile = exptDataList[[tf1]]$bwFile, name = "cov")
  peaks2 <- region_coverage(regions = peaks2, bwFile = exptDataList[[tf2]]$bwFile, name = "cov")
  
  
  ## create master peak regions
  masterSet <- GenomicRanges::union(peaks1, peaks2)
  
  ## find coverage of each sample for master peak set
  masterSet <- region_coverage(regions = masterSet, bwFile = exptDataList[[tf1]]$bwFile, name = tf1)
  masterSet <- region_coverage(regions = masterSet, bwFile = exptDataList[[tf2]]$bwFile, name = tf2)
  
  plot_MA_gg(df = as.data.frame(mcols(masterSet)), s1 = tf1, s2 = tf2, colorCol = "peaks")
  plot_scatter(df = as.data.frame(mcols(masterSet)), s1 = tf1, s2 = tf2, colorCol = "peaks") 
  
  
  ## find overlap of each peak set with master peak set
  op1 <- GenomicRanges::findOverlaps(query = masterSet, subject = peaks1)
  op2 <- GenomicRanges::findOverlaps(query = masterSet, subject = peaks2)
  
  
  # cor(x = mcols(masterSet)[[tf1]][queryHits(op1)], y = peaks1$cov[subjectHits(op1)])
  # cor(x = mcols(masterSet)[[tf2]][queryHits(op2)], y = peaks2$cov[subjectHits(op2)])
  # cor(x = mcols(masterSet)[[tf1]][queryHits(op1)], y = mcols(masterSet)[[tf2]][queryHits(op1)])
  # cor(x = mcols(masterSet)[[tf1]][queryHits(op2)], y = mcols(masterSet)[[tf2]][queryHits(op2)])
  # cor(x = mcols(masterSet)[[tf1]], y = mcols(masterSet)[[tf2]])
  
  ## set the original peak region coverage instead of master region coverage if peak is found for sample 
  mcols(masterSet)[[ tfCols[[tf1]]$cov ]] <- NA
  mcols(masterSet)[[ tfCols[[tf1]]$cov ]][queryHits(op1)] <- peaks1$cov[subjectHits(op1)]
  
  mcols(masterSet)[[ tfCols[[tf2]]$cov ]] <- NA
  mcols(masterSet)[[ tfCols[[tf2]]$cov ]][queryHits(op2)] <- peaks2$cov[subjectHits(op2)]
  
  ## replace NA values with master region coverage values
  ## and generate quartile ranks
  df <- as.data.frame(mcols(masterSet)) %>% 
    dplyr::mutate(
      !! as.name(tfCols[[tf1]]$cov) := if_else(condition = is.na(!! as.name(tfCols[[tf1]]$cov)),
                                               true = !! as.name(tfCols[[tf1]]$id),
                                               false = !! as.name(tfCols[[tf1]]$cov)),
      !! as.name(tfCols[[tf2]]$cov) := if_else(condition = is.na(!! as.name(tfCols[[tf2]]$cov)),
                                               true = !! as.name(tfCols[[tf2]]$id),
                                               false = !! as.name(tfCols[[tf2]]$cov))
    ) %>% 
    dplyr::mutate(
      !! tfCols[[tf1]]$qtBin := ntile(!! as.name(tfCols[[tf1]]$cov), 100),
      !! tfCols[[tf2]]$qtBin := ntile(!! as.name(tfCols[[tf2]]$cov), 100)
    )
  
  
  
  # plot_scatter(df = df, s1 = tfCols[[tf1]]$qtBin, s2 = tfCols[[tf2]]$qtBin,
  #              colorCol = "peaks", transformation = "identity", pseudoCount = 0)
  
  cr <- cor(x = df[, c(1, 3, 5)], y = df[, c(2, 4, 6)], method = "spearman")
  cr <- cor(x = df, method = "spearman")
  # corrplot.mixed(cr, lower.col = "black", number.cex = .7, number.digits = 3, tl.pos = "lt")
  
  corList[[sampleId]] <- c(
    corMaster = cor(x = df[[tfCols[[tf1]]$id]], y = df[[tfCols[[tf2]]$id]], method = "spearman"),
    corCov = cor(x = df[[tfCols[[tf1]]$cov]], y = df[[tfCols[[tf2]]$cov]], method = "spearman"),
    corRankBin = cor(x = df[[tfCols[[tf1]]$qtBin]], y = df[[tfCols[[tf2]]$qtBin]], method = "spearman")
  )
  
}


corDf <- as.data.frame(do.call(rbind, corList)) %>% 
  tibble::rownames_to_column(var = "sample")

readr::write_tsv(x = corDf, path = paste(outPrefix, "_data.tab", sep = ""))




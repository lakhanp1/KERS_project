suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))


## plot 

rm(list = ls())

##################################################################################

## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
comparisonName <- "kdmB_polII_corr"

tf_sample <- "An_kdmB_20h_HA_1"
polII_sample <- "An_untagged_20h_polII_1"

outDir <- here::here("analysis", "11_KERS_polII_corr", comparisonName)
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


## colors
colList <- list()


##################################################################################
## read the experiment sample details and select only those which are to be plotted
tfData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = tf_sample,
  dataPath = TF_dataPath
)


polIIData <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = polII_sample,
  dataPath = polII_dataPath
)


exptData <- dplyr::bind_rows(tfData, polIIData)

tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", tf_sample, sep = ""), names = tf_sample) },
  simplify = F, USE.NAMES = T)

polIICols <- list(
  exp = structure(polII_sample, names = polII_sample),
  is_expressed = structure(paste("is_expressed", ".", polII_sample, sep = ""), names = polII_sample)
)


## genes to read
geneSet <- suppressMessages(
  readr::read_tsv(
    file = file_genes,
    col_names = c("chr", "start", "end", "geneId", "score", "strand")
  )
)

genesGr <- rtracklayer::import.bed(con = file_genes)
tssUp <- GenomicFeatures::promoters(x = genesGr, upstream = 500, downstream = 0)
tssDown <- GenomicFeatures::promoters(x = genesGr, upstream = 0, downstream = 500)


expressionData <- get_TF_binding_data(
  exptInfo = tfData, genesDf = geneSet
)

expressionData <- get_polII_expressions(
  exptInfo = polIIData, genesDf = expressionData
)


bindingCov <- chipmine::region_coverage(regions = tssUp, bwFile = tfData$bwFile)
polIICov <- chipmine::region_coverage(regions = tssDown, bwFile = polIIData$bwFile)

if(all(bindingCov$name == polIICov$name)){
  covDf <- tibble::tibble(
    geneId = bindingCov$name,
    binding = bindingCov$coverage,
    polII = polIICov$coverage
  )
} else{
  stop("ERROR")
}

expressionData <- dplyr::left_join(
  x = expressionData, y = covDf, by = "geneId"
) %>% 
  tidyr::unite(
    col = "peakPolII", remove = F,
    !!!unname(c(tfCols$hasPeak, polIICols$is_expressed))
  )

hasPeakDf <- filter_at(
  .tbl = expressionData,
  .vars = unname(tfCols$hasPeak),
  .vars_predicate = any_vars(. == "TRUE")
)

peakExpDf <- filter_at(
  .tbl = expressionData,
  .vars = unname(c(tfCols$hasPeak, polIICols$is_expressed)),
  .vars_predicate = any_vars(. == "TRUE")
)


ggplot2::ggplot(
  data = expressionData,
  mapping = aes(x = log2(polII+1), y = log2(binding+1),
                color = peakPolII)
) +
  geom_point() +
  stat_summary_bin(
    fun = "mean", bins = 100, color='black', size=2, geom='point',
    orientation = "x"
  ) +
  # geom_bin2d(bins = 50) +
  # geom_smooth(method=glm, se = FALSE, formula = y ~ x, color = "red") +
  # geom_smooth(se = FALSE, color = "red") +
  # ggpubr::stat_cor(method = "spearman", size = 6) +
  theme_bw()








library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(TxDb.Anidulans.AspGD.GFF)
library(here)


rm(list = ls())

##################################################################################

analysisName <- "creA_target_coverage"
outDir <- here::here("kdmB_analysis", "89_creA_target_comparison")
outPrefix <- paste(outDir, "/", analysisName, sep = "")

file_creA_targets <- here::here("kdmB_analysis", "89_creA_target_comparison", "creA_targets.txt")
file_samples <- here::here("kdmB_analysis", "89_creA_target_comparison", "sample_ids.txt")

## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")
other_dataPath <- here::here("data", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.AspGD.GFF
##################################################################################

sampleList <- suppressMessages(readr::read_tsv(file = file_samples, col_names = T, comment = "#"))
creA_targets <- suppressMessages(
  readr::read_tsv(file = file_creA_targets, col_names = T, comment = "#"))

## genes to read
geneSet <- suppressMessages(
  readr::read_tsv(file = file_genes, col_names = c("chr", "start", "end", "geneId", "score", "strand"))
) %>% 
  dplyr::select(geneId)

geneDesc <- suppressMessages(
  AnnotationDbi::select(x = orgDb, keys = geneSet$geneId,
                        columns = c("DESCRIPTION"), keytype = "GID")) %>% 
  dplyr::rename(geneId = GID)

creA_targets <- dplyr::left_join(x = creA_targets, y = geneDesc, by = "geneId")

##################################################################################
## prepare sample information dataframe

## read the experiment sample details and select only those which are to be plotted
tempSInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = sampleList$sampleId,
                                    dataPath = TF_dataPath)

# polII_ids <- tempSInfo$sampleId[which(tempSInfo$data_type == "polII")]
tfIds <- tempSInfo$sampleId[which(tempSInfo$data_type == "TF" & tempSInfo$TF != "untagged")]
inputIds <- tempSInfo$sampleId[which(tempSInfo$data_type == "TF" & tempSInfo$TF == "untagged")]
# histIds <- tempSInfo$sampleId[which(tempSInfo$data_type == "HIST")]
otherIds <- tempSInfo$sampleId[which(tempSInfo$data_type == "other")]


## read the experiment sample details and select only those which are to be plotted
tfData <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = tfIds,
                                 dataPath = TF_dataPath)


inputData <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = inputIds,
                                    dataPath = TF_dataPath)

otherData <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = otherIds,
                                    dataPath = other_dataPath)

exptData <- dplyr::bind_rows(tfData, inputData, otherData)

exptDataList <- purrr::transpose(exptData) %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

##################################################################################
## get coverage matrix

genesGr <- GenomicFeatures::genes(x = txDb, columns = "gene_id",
                                  filter = list(gene_id = creA_targets$geneId))

upstreamGr <- GenomicRanges::promoters(x = genesGr, upstream = 500, downstream = 100,
                                       use.names = FALSE)

# mcols(upstreamGr)$name <- mcols(upstreamGr)$gene_id
# mcols(upstreamGr)$gene_id <- NULL
mcols(upstreamGr)$region <- paste(
  seqnames(upstreamGr), ":", start(upstreamGr), "-", end(upstreamGr), sep = "")

upstreamCov <- region_coverage_matrix(regions = upstreamGr, exptInfo = exptData)

creA_targets <- dplyr::left_join(x = creA_targets, y = upstreamCov, by = c("geneId" = "gene_id"))

readr::write_tsv(x = creA_targets, path = paste(outPrefix, ".tab", sep = ""))


##################################################################################
## heatmap plots

covMat <- as.matrix(dplyr::select(creA_targets, exptData$sampleId))
rownames(covMat) <- creA_targets$geneId

quantile(covMat, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)

cov_color <- colorRamp2(
  breaks = c(0, quantile(covMat, c(0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99))),
  colors = c("white", RColorBrewer::brewer.pal(n = 9, name = "YlOrRd"))
)


ht <- Heatmap(
  matrix = covMat,
  col = cov_color,
  name = "CreA_targets",
  column_title = "CreA target gene's KERS and PIC binding profile",
  split = creA_targets$group,
  cluster_columns = F,
  show_row_names = F,
  show_row_dend = F,
  row_title_rot = 0,
  border = TRUE
)


png(filename = paste(outPrefix, ".heatmap.png", sep = ""), width = 2000, height = 2500, res = 200)

draw(
  ht,
  row_title = "CreA targets"
)

dev.off()














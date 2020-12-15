library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(TxDb.Anidulans.FGSCA4.AspGD.GFF)
library(here)

## KERS 20h+48h samples target gene matrix 

rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/04_GO_enrichment/topGO_functions.R")


## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
analysisName <- "KERS_targets"
workDir <- here::here("kdmB_analysis", "03_KERS_complex", analysisName)
outPrefix <- paste(workDir, "/", analysisName, sep = "")

file_plotSamples <- paste(workDir, "/", "samples.txt", sep = "")
orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

matrixType <- "normalizedMatrix_5kb"
up <- 5000
body <- 2000
down <- 1000
binSize <- 10
matrixDim = c(c(up, body, down)/binSize, binSize)

showExpressionHeatmap = FALSE

## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")

anLables <- list()

##################################################################################

sampleList <- suppressMessages(readr::read_tsv(file = file_plotSamples, col_names = T, comment = "#"))

## genes to read
geneSet <- suppressMessages(
  readr::read_tsv(file = file_genes, col_names = c("chr", "start", "end", "geneId", "score", "strand"))
) %>% 
  dplyr::mutate(length = end - start) %>% 
  dplyr::select(geneId)

geneDesc <- suppressMessages(
  AnnotationDbi::select(x = orgDb, keys = geneSet$geneId,
                        columns = c("DESCRIPTION"), keytype = "GID"))

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

head(geneInfo)


##################################################################################

## read the experiment sample details and select only those which are to be plotted
tempSInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = sampleList$sampleId,
                                    dataPath = TF_dataPath,
                                    profileMatrixSuffix = matrixType)

polII_ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "polII")]
tfIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag %in% c("HA", "MYC", "TAP") & tempSInfo$TF != "untagged")]
inputIds <- tempSInfo$sampleId[which(! tempSInfo$IP_tag %in% c("polII", "HIST") & tempSInfo$TF == "untagged")]
histIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST")]


## read the experiment sample details and select only those which are to be plotted
tfData <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = tfIds,
                                 dataPath = TF_dataPath,
                                 profileMatrixSuffix = matrixType)


inputData <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = inputIds,
                                    dataPath = TF_dataPath,
                                    profileMatrixSuffix = matrixType)

polIIData <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = polII_ids,
                                    dataPath = polII_dataPath,
                                    profileMatrixSuffix = matrixType)

histData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = histIds,
                                   dataPath = hist_dataPath,
                                   profileMatrixSuffix = matrixType)

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

##################################################################################

bindingMat = combinatorial_binding_matrix(sampleInfo = tfData, summitRegion = 50,
                                          peakCols = c("peakId"))

bindingMat$region <- paste(bindingMat$seqnames, ":", bindingMat$start, "-", bindingMat$end, sep = "")

hist(bindingMat$end - bindingMat$start)

##################################################################################
## binding stats

i <- 1
sampleInfo = tfData

regionMat <- dplyr::select(bindingMat, name)

for (i in 1:nrow(sampleInfo)) {
  
  peakAn <- import_peak_annotation(
    sampleId = tfData$sampleId[i], peakAnnoFile = tfData$peakAnno[i],
    columns = c("peakId", "geneId", "peakPval"))
  
  joinby <- unname(tfCols$peakId[i])
  peakRegionAn <- dplyr::left_join(x = bindingMat, y = peakAn, by = joinby) %>% 
    dplyr::select(name, geneId, colnames(peakAn)) %>% 
    dplyr::filter_at(.vars = vars(starts_with("peakId.")),
                     .vars_predicate = any_vars(!is.na(.)))
  
  if(i == 1){
    regionMat <- dplyr::left_join(x = regionMat, y = peakRegionAn, by = c("name"))
  } else{
    regionMat <- dplyr::full_join(x = regionMat, y = peakRegionAn, by = c("name", "geneId"))
  }
  
}

targetMatrix <- dplyr::left_join(
  x = dplyr::select(bindingMat, name, region),
  y = regionMat, by = "name") %>% 
  # dplyr::arrange(seqnames, start) %>%
  dplyr::add_count(name, name = "n") %>% 
  dplyr::select(name, region, n, everything()) %>% 
  dplyr::filter_at(.vars = vars(starts_with("peakId.")),
                   .vars_predicate = any_vars((!is.na(.) | n == 1))) %>%
  dplyr::add_count(name, name = "n") %>%
  dplyr::distinct()

readr::write_tsv(x = targetMatrix, path = paste(outPrefix, ".summit_based.tab", sep = ""))

##################################################################################

peakTargetMat <- peak_target_matrix(sampleInfo = tfData, position = "both")

readr::write_tsv(x = peakTargetMat, path = paste(outPrefix, ".summit_based.genes.tab", sep = ""))




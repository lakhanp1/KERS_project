suppressPackageStartupMessages(library(chipmine))
# suppressPackageStartupMessages(library(markPeaks))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(here))

## generate profile matrix for regions of interest from bigWig files

rm(list = ls())

cl <- makeCluster(4) #not to overload your computer
registerDoParallel(cl)

##################################################################################

file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
txDb <- TxDb.Anidulans.FGSCA4.AspGD.GFF

file_tf_macs2 <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")
file_tf <- paste(TF_dataPath, "/", "sample_tf.list", sep = "")

matrixType <- "normalizedmatrix"
up <- 5000
down <- 1000
body <- 2000
bin <- 10
matrixDim = c(c(up, body, down)/bin, bin)

# geneSet <- suppressMessages(
#   readr::read_tsv(
#     file = file_genes,col_names = c("chr", "start", "end", "geneId", "score", "strand")
#   )) %>%
#   dplyr::select(geneId)
# 
# geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$geneId, columns = "DESCRIPTION", keytype = "GID")
# 
# geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID"))

##################################################################################
sampleList <- tibble::tibble(
  sampleId = c(
    # "An_kdmB_20h_HA_1", "An_sntB_20h_HA_1", "An_H3_20h_HIST_1", "An_H3K4Me3_20h_HIST_1",
    # "An_H3K9Ac_20h_HIST_1", "An_H3K9Me3_20h_HIST_1", "An_H3K36Me3_20h_HIST_1", "An_H3_48h_HIST_1", 
    # "An_H3K9Ac_48h_HIST_1", "An_rpdA_20h_HA_1", "An_sntB_48h_HA_1",  "An_kdmB_48h_HA_1",
    # "An_H3K36Me3_48h_HIST_1", "An_H3K4Me3_48h_HIST_1", "An_H3K9Me3_48h_HIST_1",
    # "An_ecoA_20h_HA_1", "An_ecoA_20h_HA_2", "An_ecoA_48h_HA_1", "An_ecoA_48h_HA_2
    # An_kdmB_20h_HA_2", "An_kdmB_48h_HA_2", "An_rpdA_20h_HA_2", "An_rpdA_48h_HA_1",
    # "An_rpdA_48h_HA_2", "An_sntB_20h_HA_2", "An_sntB_48h_HA_2", "An_untagged_20h_HA_1",
    # "An_untagged_20h_HA_2", "An_untagged_48h_HA_1", "An_untagged_48h_HA_2",
    "H3_ANM_NH4_CW360_411", "H3_NH4_CW466_525", "H3K27ac_NH4_CW466_525", 
    "H3K27ac_ANM_NH4_CW360_411", "H3K4me3_NH4_CW466_525", "H3K4me3_ANM_NH4_CW360_411",
    "H3K9ac_NH4_CW466_525", "H3K9ac_ANM_NH4_CW360_411", "H3ac_NH4_CW466_525",
    "veA1_Rpb1_15h_mix22_1", "veA1_Rpb3_15h_mix22_1", "veA1_TFIIB_15h_mix22_1",
    "veA1_TBP_15h_mix22_1", "veA_wt_Rpb1_15h_mix22_1", "veA_wt_Rpb1_15h_mix22_2",
    "veA_wt_Rpb3_15h_mix22_1", "veA_wt_Rpb3_15h_mix22_2", "veA_wt_TBP_15h_mix22_1",
    "veA_wt_TBP_15h_mix22_2", "veA_wt_TFIIB_15h_mix22_2", "veA_wt_TFIIB_15h_mix22_1",
    "WT_H3K27ac_20h_1", "WT_H3K27ac_20h_2", "WT_H3K27ac_20h_3", "WT_H3K27ac_48h_1",
    "WT_H3K27ac_48h_2" 
  )
)


tempSInfo <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = sampleList$sampleId,
  dataPath = TF_dataPath, profileMatrixSuffix = matrixType
)

tfIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag %in% c("HA", "MYC", "TAP") & tempSInfo$TF != "untagged")]
inputIds <- tempSInfo$sampleId[which(! tempSInfo$IP_tag %in% c("polII", "HIST") & tempSInfo$TF == "untagged")]
polII_ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "polII")]
histIds <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST" & tempSInfo$TF != "H3")]
histH3Ids <- tempSInfo$sampleId[which(tempSInfo$IP_tag == "HIST" & tempSInfo$TF == "H3")]

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

h3Data <- get_sample_information(
  exptInfoFile = file_exptInfo, samples = histH3Ids,
  dataPath = hist_dataPath, profileMatrixSuffix = matrixType
)



exptData <- dplyr::bind_rows(tfData, inputData, histData, h3Data, polIIData)

exptDataList <- purrr::transpose(exptData)  %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))


genesGr <- rtracklayer::import.bed(con = file_genes)
tssGr <- GenomicFeatures::promoters(x = genesGr, upstream = 0, downstream = 1)
tesGr <- chipmine::get_TES(gr = genesGr, up = 0, down = 1)
##################################################################################

rowId <- 1

foreach::foreach(
  rowId = 1:nrow(exptData),
  .packages = c("chipmine")
  ) %dopar% {
# for(rowId in 1:nrow(exptData)){
  
  # ## 2kb - 2kb - 1kb matrix
  # bwMat <- chipmine::bigwig_profile_matrix(
  #   bwFile = exptData$bwFile[rowId],
  #   regions = file_genes,
  #   signalName = exptData$sampleId[rowId],
  #   storeLocal = TRUE,
  #   localPath = exptData$matFile[rowId]
  # )

  # mat5Kb <- gsub(pattern = ".tab.gz", replacement = "_5kb.tab.gz", x = exptData$matFile[rowId])
  # 
  # ## 5kb - 2kb - 1kb matrix
  # bwMat <- chipmine::bigwig_profile_matrix(
  #   bwFile = exptData$bwFile[rowId],
  #   regions = file_genes,
  #   signalName = exptData$sampleId[rowId],
  #   storeLocal = TRUE,
  #   localPath = mat5Kb,
  #   extend = c(5000, 1000),
  #   target_ratio = 0.25
  # )
  
  ## -3kb - ATG - 3kb
  matTss <- gsub(pattern = ".tab.gz", replacement = "_3kbATG3kb.tab.gz", x = exptData$matFile[rowId])
  bwMat <- chipmine::bigwig_profile_matrix(
    bwFile = exptData$bwFile[rowId],
    regions = tssGr,
    signalName = exptData$sampleId[rowId],
    extend = c(3000, 3000),
    storeLocal = TRUE,
    localPath = matTss,
    targetName = "tss"
  )
  
  # ## -2kb - TES - 2kb
  # matTes <- gsub(pattern = ".tab.gz", replacement = "_2kbTES2kb.tab.gz", x = exptData$matFile[rowId])
  # bwMat <- chipmine::bigwig_profile_matrix(
  #   bwFile = exptData$bwFile[rowId],
  #   regions = tesGr,
  #   signalName = exptData$sampleId[rowId],
  #   extend = c(2000, 2000),
  #   storeLocal = TRUE,
  #   localPath = matTss,
  #   targetName = "tes"
  # )
  #
  # matTss <- gsub(pattern = ".tab.gz", replacement = "_3kbTSS3kb.tab.gz", x = exptData$matFile[rowId])
  # ## -5kb - TSS - 1kb
  # bwMat <- chipmine::bigwig_profile_matrix(
  #   bwFile = exptData$bwFile[rowId],
  #   regions = tssGr,
  #   signalName = exptData$sampleId[rowId],
  #   storeLocal = TRUE,
  #   localPath = matTss,
  #   extend = c(3000, 3000),
  #   targetName = "tss"
  # )
  #
  #
  # matTes <- gsub(pattern = ".tab.gz", replacement = "_3kbTES3kb.tab.gz", x = exptData$matFile[rowId])
  # ## -1kb - TES - 5kb
  # bwMat <- chipmine::bigwig_profile_matrix(
  #   bwFile = exptData$bwFile[rowId],
  #   regions = tesGr,
  #   signalName = exptData$sampleId[rowId],
  #   storeLocal = TRUE,
  #   localPath = matTes,
  #   extend = c(3000, 3000),
  #   targetName = "tes"
  # )
  # 
  
  # print(exptData$sampleId[rowId])
  
  
}


parallel::stopCluster(cl = cl)
##################################################################################

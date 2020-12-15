suppressPackageStartupMessages(library(chipmine))
# suppressPackageStartupMessages(library(markPeaks))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(TxDb.Anidulans.FGSCA4.AspGD.GFF))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(here))

## 1) annotate peaks
## 2) create gene level peak annotation data

rm(list = ls())

# cl <- makeCluster(4) #not to overload your computer
# registerDoParallel(cl)

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

geneSet <- suppressMessages(
  readr::read_tsv(
    file = file_genes,col_names = c("chr", "start", "end", "geneId", "score", "strand")
  )) %>%
  dplyr::select(geneId)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$geneId, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID"))

##################################################################################

tfSampleList <- suppressMessages(
  readr::read_tsv(file = file_tf_macs2, col_names = c("id"),  comment = "#")
)

txInfo <- suppressMessages(
  AnnotationDbi::select(
    x = txDb, keys = AnnotationDbi::keys(x = txDb, keytype = "TXID"),
    columns = c("GENEID", "TXNAME", "TXTYPE"), keytype = "TXID")) %>%
  dplyr::mutate(TXID = as.character(TXID)) %>%
  dplyr::rename(geneId = GENEID, txName = TXNAME, txType = TXTYPE)

txInfo <- dplyr::filter(txInfo, !txType %in% c("tRNA", "rRNA", "snRNA", "snoRNA")) %>% 
  dplyr::filter(!grepl(pattern = "uORF", x = geneId))

# tfSampleList <- data.frame(id = c("An_untagged_48h_input_1", "An_untagged_20h_HA_1", "An_untagged_20h_HA_2", "An_untagged_48h_HA_1", "An_untagged_48h_HA_2", "An_untagged_20h_MYC_1", "An_untagged_20h_MYC_2", "An_untagged_48h_MYC_1", "An_untagged_48h_MYC_2"),
#                           stringsAsFactors = F)

# tfSampleList <- data.frame(id = c("An_kdmB_20h_HA_1", "An_kdmB_laeA_del_20h_HA_2", "An_sudA_20h_HA_1", "An_rpdA_20h_HA_1", "An_sntB_20h_HA_1", "An_ecoA_20h_HA_1", "An_kdmB_48h_HA_1", "An_rpdA_48h_HA_1", "An_sntB_48h_HA_1", "An_sudA_48h_HA_1", "An_ecoA_48h_HA_1"),
#                            stringsAsFactors = F)

tfInfo <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$id,
  dataPath = TF_dataPath,
  profileMatrixSuffix = matrixType)


i <- 1

for(i in 1:nrow(tfInfo)){
  
  ## annotate peaks and prepare gene level annotation file
  peakType <- dplyr::case_when(
    tfInfo$peakType[i] == "narrow" ~ "narrowPeak",
    tfInfo$peakType[i] == "broad" ~ "broadPeak"
  )
  
  cat("Annotating", tfInfo$sampleId[i],"\n")
  
  peakAn <- narrowPeak_annotate(
    peakFile = tfInfo$peakFile[i],
    txdb = txDb,
    summitRegion = 1,
    txIds = txInfo$TXID,
    fileFormat = peakType,
    promoterLength = 700,
    upstreamLimit = 1500,
    bidirectionalDistance = 500,
    bidirectionalSkew = 0.2,
    includeFractionCut = 0.7,
    bindingInGene = FALSE,
    insideSkewToEndCut = 0.7,
    output = tfInfo$peakAnno[i],
    removePseudo = TRUE)
  
  ## based on the updated package markPeaks
  # peakAn <- markPeaks::annotate_peaks(
  #   peakFile = tfInfo$peakFile[i],
  #   txdb = txDb,
  #   txIds = txInfo$TXID,
  #   summitRegion = 1,
  #   fileFormat = peakType,
  #   promoterLength = 700,
  #   upstreamLimit = 1500,
  #   bidirectionalDistance = 500,
  #   includeFractionCut = 0.7,
  #   bindingInGene = FALSE,
  #   insideSkewToEndCut = 0.7,
  #   removePseudo = TRUE,
  #   output = tfInfo$peakAnno[i]
  # )
  
  # if( !is.null(peakAn) ){
  #   tfDf <- gene_level_peak_annotation(
  #     sampleId = tfInfo$sampleId[i],
  #     peakAnnotation = tfInfo$peakAnno[i],
  #     genesDf = geneSet,
  #     peakFile = tfInfo$peakFile[i],
  #     bwFile = tfInfo$bwFile[i],
  #     outFile = tfInfo$peakTargetFile[i])
  # }
  
  
  
  # ## 2kb - 2kb - 1kb matrix
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                          bedFile = file_genes,
  #                                          signalName = tfInfo$sampleId[i],
  #                                          genes = geneSet$gene,
  #                                          readLocal = FALSE,
  #                                          storeLocal = TRUE,
  #                                          localPath = tfInfo$matFile[i])
  
  # mat5Kb <- gsub(pattern = ".tab.gz", replacement = "_5kb.tab.gz", x = tfInfo$matFile[i])
  # 
  # ## 5kb - 2kb - 1kb matrix
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                          bedFile = file_genes,
  #                                          signalName = tfInfo$sampleId[i],
  #                                          genes = geneSet$gene,
  #                                          readLocal = FALSE,
  #                                          storeLocal = TRUE,
  #                                          localPath = mat5Kb,
  #                                          extend = c(5000, 1000),
  #                                          target_ratio = 0.25)
  
  # matTss <- gsub(pattern = ".tab.gz", replacement = "_4kbTSS2kb.tab.gz", x = tfInfo$matFile[i])
  # ## -4kb - TSS - 2kb
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                          bedFile = file_genes,
  #                                          signalName = tfInfo$sampleId[i],
  #                                          genes = geneSet$gene,
  #                                          readLocal = FALSE,
  #                                          storeLocal = TRUE,
  #                                          localPath = matTss,
  #                                          extend = c(4000, 2000),
  #                                          target = "tss")
  #
  # matTes <- gsub(pattern = ".tab.gz", replacement = "_2kbTES4kb.tab.gz", x = tfInfo$matFile[i])
  # ## -2kb - TES - 4kb
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                          bedFile = file_genes,
  #                                          signalName = tfInfo$sampleId[i],
  #                                          genes = geneSet$gene,
  #                                          readLocal = FALSE,
  #                                          storeLocal = TRUE,
  #                                          localPath = matTes,
  #                                          extend = c(2000, 4000),
  #                                          target = "tes")
  #
  # matTss <- gsub(pattern = ".tab.gz", replacement = "_3kbTSS3kb.tab.gz", x = tfInfo$matFile[i])
  # ## -5kb - TSS - 1kb
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                          bedFile = file_genes,
  #                                          signalName = tfInfo$sampleId[i],
  #                                          genes = geneSet$gene,
  #                                          readLocal = FALSE,
  #                                          storeLocal = TRUE,
  #                                          localPath = matTss,
  #                                          extend = c(3000, 3000),
  #                                          target = "tss")
  #
  #
  # matTes <- gsub(pattern = ".tab.gz", replacement = "_3kbTES3kb.tab.gz", x = tfInfo$matFile[i])
  # ## -1kb - TES - 5kb
  # bwMat <- chipmine::bigwig_profile_matrix(bwFile = tfInfo$bwFile[i],
  #                                          bedFile = file_genes,
  #                                          signalName = tfInfo$sampleId[i],
  #                                          genes = geneSet$gene,
  #                                          readLocal = FALSE,
  #                                          storeLocal = TRUE,
  #                                          localPath = matTes,
  #                                          extend = c(3000, 3000),
  #                                          target = "tes")
  #
  
  # print(tfInfo$sampleId[i])
  
  
}



# ##################################################################################
#
#
# ## specific processing for samples where binding is seen in gene body
# # tfSampleList = c("An_laeA_20h_HA", "An_laeA_48h_HA", "An_kdmB_20h_HA", "An_kdmB_48h_HA")

samplesWithBindingInGene <- c("An_cclA_20h_HA_1", "An_cclA_20h_HA_2", "An_cclA_48h_HA_1", "An_cclA_48h_HA_2", "An_cclA_kdmA_del_20h_HA_1", "An_cclA_kdmA_del_20h_HA_2", "An_cclA_kdmA_del_48h_HA_1", "An_cclA_kdmA_del_48h_HA_2")


tfInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = samplesWithBindingInGene,
                                 dataPath = TF_dataPath,
                                 profileMatrixSuffix = "normalizedmatrix")

i <- 1

for (i in 1:nrow(tfInfo)) {
  
  peakType <- dplyr::case_when(
    tfInfo$peakType[i] == "narrow" ~ "narrowPeak",
    tfInfo$peakType[i] == "broad" ~ "broadPeak"
  )
  
  peakAn <- narrowPeak_annotate(
    # peakFile = tfInfo$peakFile[i],
    # txdb = txDb,
    # txIds = txInfo$TXID,
    # fileFormat = peakType,
    # includeFractionCut = 0.7,
    # bindingInGene = TRUE,
    # promoterLength = 700,
    # upstreamLimit = 1500,
    # insideSkewToEndCut = 0.7,
    # removePseudo = FALSE,
    # output = tfInfo$peakAnno[i]
    peakFile = tfInfo$peakFile[i],
    txdb = txDb,
    summitRegion = 1,
    txIds = txInfo$TXID,
    fileFormat = peakType,
    promoterLength = 700,
    upstreamLimit = 1500,
    bidirectionalDistance = 500,
    bidirectionalSkew = 0.2,
    includeFractionCut = 0.7,
    bindingInGene = TRUE,
    insideSkewToEndCut = 0.7,
    output = tfInfo$peakAnno[i],
    removePseudo = TRUE
  )
  
  if( !is.null(peakAn) ){
    tfDf <- gene_level_peak_annotation(
      sampleId = tfInfo$sampleId[i],
      peakAnnotation = tfInfo$peakAnno[i],
      preference = c("nearStart", "peakInFeature", "featureInPeak",
                     "nearEnd", "upstreamTss", "intergenic"),
      genesDf = geneSet,
      peakFile = tfInfo$peakFile[i],
      bwFile = tfInfo$bwFile[i],
      outFile = tfInfo$peakTargetFile[i])
    
  }
  
}








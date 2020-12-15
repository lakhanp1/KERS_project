library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(BSgenome.Anidulans.FGSCA4.AspGD)
library(here)


rm(list = ls())


## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
analysisName <- "motifs_KERS_20h"
workDir <- here::here("kdmB_analysis", "03_KERS_complex", "motif_analysis", analysisName)
outPrefix <- paste(workDir, "/", analysisName, sep = "")

file_plotSamples <- paste(workDir, "/", "samples.txt", sep = "")
orgDb <- org.Anidulans.FGSCA4.eg.db
bsGnom <- BSgenome.Anidulans.FGSCA4.AspGD


## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/annotation_resources/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")

##################################################################################

sampleList <- suppressMessages(readr::read_tsv(file = file_plotSamples, col_names = T, comment = "#"))

## genes to read
geneSet <- suppressMessages(
  readr::read_tsv(file = file_genes, col_names = c("chr", "start", "end", "geneId", "score", "strand"))
) %>% 
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
tfData <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = sampleList$sampleId,
                                 dataPath = TF_dataPath)

tfIds <- tfData$sampleId[which(tfData$IP_tag %in% c("HA", "MYC", "TAP") & tfData$TF != "untagged")]

exptData <- tfData

exptDataList <- purrr::transpose(exptData) %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))


tfCols <- sapply(
  c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
    "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
    "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
  simplify = F, USE.NAMES = T)



# dplyr::group_by_at(hasPeakDf, .vars = vars(starts_with("hasPeak."))) %>% 
#   dplyr::summarise(n = n_distinct(geneId))
# 
# ## group label data mean profile facet plot
# groupLabelDf <- dplyr::group_by_at(.tbl = hasPeakDf, .vars = vars(starts_with("hasPeak."), group)) %>%
#   dplyr::summarise(n = n()) %>% 
#   dplyr::mutate(groupLabels = paste(group, ": ", n, " genes", sep = ""))
# 
# groupLabels <- structure(groupLabelDf$groupLabels, names = groupLabelDf$group)

##################################################################################
## binding stats
bindingMat = combinatorial_binding_matrix(
  sampleInfo = tfData, summitRegion = 50,
  peakCols = c("peakId", "peakPval"),
  genome = bsGnom, summitSeqLen = 100
)


commonBinding <- dplyr::filter_at(
  .tbl = bindingMat, .vars = vars(starts_with("overlap.")),
  .vars_predicate = all_vars(. == TRUE))

readr::write_tsv(x = commonBinding, path = paste(outPrefix, ".data.tab", sep = ""))


seqLen <- 500

i <- 1

for (i in 1:nrow(tfData)) {
  
  summitSeq <- get_peak_summit_seq(
    file = tfData$peakFile[i], sampleId = tfData$sampleId[i],
    genome = bsGnom, length = seqLen, column_name_prefix = FALSE)
  
  peaksDf <- import_peaks_as_df(file = tfData$peakFile[i], sampleId = tfData$sampleId[i],
                                rename = F)
  
  summitSeq <- dplyr::left_join(x = summitSeq, y = peaksDf, by = "peakId") %>% 
    dplyr::filter(peakPval >= 20)
  
  seqSet <- DNAStringSet(x = summitSeq$summitSeq)
  names(seqSet) <- summitSeq$peakId
  
  file_fasta <- paste(workDir, "/", tfData$sampleId[i], ".summitSeq.",seqLen, "bp.fasta", sep = "")
  
  writeXStringSet(x = seqSet, filepath = file_fasta)
}










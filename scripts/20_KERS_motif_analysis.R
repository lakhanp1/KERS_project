suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(BSgenome.Anidulans.FGSCA4.AspGD))
suppressPackageStartupMessages(library(here))


rm(list = ls())

##################################################################################
## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
analysisName <- "motifs_KERS_48h"
workDir <- here::here("analysis", "03_KERS_complex", "motif_analysis", analysisName)
outPrefix <- paste(workDir, "/", analysisName, sep = "")

seqLen <- 500

file_plotSamples <- paste(workDir, "/", "samples.txt", sep = "")
orgDb <- org.Anidulans.FGSCA4.eg.db
bsGnom <- BSgenome.Anidulans.FGSCA4.AspGD

## genes to read
file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

##################################################################################

sampleList <- suppressMessages(readr::read_tsv(file = file_plotSamples, comment = "#"))

## read the experiment sample details and select only those which are to be plotted

tfData <- get_sample_information(
  exptInfoFile = file_exptInfo, dataPath = TF_dataPath, samples = sampleList$sampleId
)

tfIds <- tfData$sampleId[which(tfData$IP_tag %in% c("HA", "MYC", "TAP") & tfData$TF != "untagged")]

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

readr::write_tsv(x = commonBinding, file = paste(outPrefix, ".data.tab", sep = ""))



rowId <- 1

for (rowId in 1:nrow(tfData)) {
  
  summitSeq <- get_peak_summit_seq(
    file = tfData$peakFile[rowId], sampleId = tfData$sampleId[rowId],
    genome = bsGnom, length = seqLen, column_name_prefix = FALSE)
  
  peaksDf <- import_peaks_as_df(file = tfData$peakFile[rowId], sampleId = tfData$sampleId[rowId],
                                rename = F)
  
  summitSeq <- dplyr::left_join(x = summitSeq, y = peaksDf, by = c("name"="peakId")) %>% 
    dplyr::filter(peakPval >= 20)
  
  seqSet <- DNAStringSet(x = summitSeq$summitSeq)
  names(seqSet) <- summitSeq$peakId
  
  file_fasta <- paste(workDir, "/", tfData$sampleId[rowId], ".summitSeq.",seqLen, "bp.fasta", sep = "")
  
  writeXStringSet(x = seqSet, filepath = file_fasta)
}










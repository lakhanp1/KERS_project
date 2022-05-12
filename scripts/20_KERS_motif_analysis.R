suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(BSgenome.Anidulans.FGSCA4.AspGD))
suppressPackageStartupMessages(library(here))


rm(list = ls())

##################################################################################
## IMP: the first sampleID will be treated primary and clustering will be done/used for/of this sample
analysisName <- "KERS_20h"
workDir <- here::here("analysis", "03_KERS_complex", "motif_analysis", analysisName)
outPrefix <- paste(workDir, "/", analysisName, sep = "")

seqLen <- 500
cutoff_pval <- 10

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


tfCols <- sapply(
  X = c("hasPeak", "peakId", "peakEnrichment", "peakPval", "peakQval", "peakSummit", "peakDist", "summitDist",
        "peakType", "bidirectional", "featureCovFrac", "relativeSummitPos", "peakRegion", "peakPosition",
        "peakCoverage", "pvalFiltered", "summitSeq"),
  FUN = function(x){ structure(paste(x, ".", tfIds, sep = ""), names = tfIds) },
  simplify = F, USE.NAMES = T)

##################################################################################
## extract common KERS sites and summit sequences
bindingMat = combinatorial_binding_matrix(
  sampleInfo = tfData, summitRegion = 50,
  peakCols = c("peakId", "peakPval"),
  genome = bsGnom, summitSeqLen = seqLen
)

readr::write_tsv(
  file = paste(outPrefix, ".binding_matrix.tab", sep = ""), x = bindingMat
)

commonBinding <- dplyr::filter_at(
  .tbl = bindingMat, .vars = vars(starts_with("overlap.")),
  .vars_predicate = all_vars(. == TRUE))

readr::write_tsv(
  x = commonBinding,
  file = paste(outPrefix, ".common_sites_summitSeq", seqLen, ".data.tab", sep = "")
)

seqSet <- NULL
seqSet <- DNAStringSet(x = commonBinding[[tfCols$summitSeq[1]]])
names(seqSet) <- commonBinding[[tfCols$peakId[1]]]
writeXStringSet(
  x = seqSet, filepath = paste(outPrefix, ".common_sites_summitSeq", seqLen, ".fasta", sep = "")
)

##################################################################################
## extract summit sequence for individual member
rowId <- 1

for (rowId in 1:nrow(tfData)) {
  
  summitSeq <- get_peak_summit_seq(
    file = tfData$peakFile[rowId], sampleId = tfData$sampleId[rowId],
    genome = bsGnom, length = seqLen, column_name_prefix = FALSE)
  
  peaksDf <- import_peaks_as_df(file = tfData$peakFile[rowId], sampleId = tfData$sampleId[rowId],
                                rename = F)
  
  summitSeq <- dplyr::left_join(x = summitSeq, y = peaksDf, by = c("name"="peakId")) %>%
    dplyr::filter(peakPval >= cutoff_pval)
  
  seqSet <- NULL
  seqSet <- DNAStringSet(x = summitSeq$summitSeq)
  names(seqSet) <- summitSeq$name
  
  file_fasta <- 
    
    writeXStringSet(
      x = seqSet,
      filepath = paste(
        workDir, "/", tfData$sampleId[rowId], ".summitSeq.pval_", cutoff_pval,
        ".", seqLen, "bp.fasta", sep = ""
      )
    )
}


##################################################################################









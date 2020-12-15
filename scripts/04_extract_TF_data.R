library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)



rm(list = ls())

path <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/kdmB_analysis"
setwd(path)


dataPath <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data"


file_exptInfo <- paste(dataPath, "referenceData/sampleInfo.txt", sep = "/")
TF_dataPath <- paste(dataPath, "TF_data", sep = "/")
polII_dataPath <- paste(dataPath, "polII_data", sep = "/")
file_genes <-  paste(dataPath, "referenceData/AN_genesForPolII.bed", sep = "/")


geneInfoFile <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"


sampleList <- c("An_kdmB_20h_HA_1", "An_kdmB_20h_HA_2", "An_kdmB_48h_HA_1", "An_kdmB_48h_HA_2")

orgDb <- org.Anidulans.FGSCA4.eg.db

##################################################################################
## read the experiment sample details and select only those which are to be plotted
exptData <- get_sample_information(exptInfoFile = file_exptInfo,
                                   samples = sampleList,
                                   dataPath = TF_dataPath,
                                   matrixSource = "deeptools")

## prepare the gene set
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::select(-score) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- select(x = orgDb, keys = geneSet$gene, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))


colsToSel <- lapply(X = exptData$sampleId,
                    FUN = function(x) {
                      paste(c("hasPeak", "peakType", "hasTesPeak", "tesPeakType"), ".", x, sep = "") 
                    }
)


tfData <- get_TF_binding_data(exptInfo = exptData,
                              genesDf = geneSet,
                              allColumns = TRUE) %>% 
  dplyr::select(gene, DESCRIPTION, !!! unlist(colsToSel))


fwrite(x = tfData, file = "kdmB_peaks.tab", sep = "\t", col.names = T, quote = F, na = "NA")


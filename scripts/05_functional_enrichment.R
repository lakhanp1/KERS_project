suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(org.Anidulans.FGSCA4.eg.db))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(foreach))


## 1) perform the GO enrichment using topGO for the genes which show binding signal
## 2) perform the pathway enrichment using keggprofile for the genes which show binding signal

rm(list = ls())

source(file = "https://raw.githubusercontent.com/lakhanp1/omics_utils/master/04_GO_enrichment/s01_enrichment_functions.R")

# path <- "E:/Chris_UM/Analysis/21_CL2017_ChIPmix_ULAS_MIX/ULAS_AN/data/"
# setwd(path)


## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"

TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")

orgDb <- org.Anidulans.FGSCA4.eg.db
keggOrgCode <- "ani"

##################################################################################

## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))



##################################################################################
## TF samples functional enrichment

tfSampleFile <- paste(TF_dataPath, "/", "tf_macs2_samples.list", sep = "")
# tfSampleFile <- paste(TF_dataPath, "/", "tf_samples.list", sep = "")

tfSampleList <- readr::read_tsv(file = tfSampleFile, col_names = c("id"),  comment = "#") %>% 
  as.data.frame()

# 
# tfSampleList <- data.frame(id = c("An_kdmA_20h_HA_1", "An_cclA_20h_HA_1", "An_mcmA_20h_MYC_1", "An_ecmB_20h_HA_1", "An_rstB_20h_HA_1", "An_kdmA_48h_HA_1", "An_cclA_48h_HA_1", "An_mcmA_48h_MYC_1", "An_ecmB_48h_HA_1", "An_rstB_48h_HA_1"),
#                            stringsAsFactors = F)




tfInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = tfSampleList$id,
                                 dataPath = TF_dataPath,
                                 matrixSource = "normalizedmatrix")

i <- 66

foreach(i = 1:nrow(tfInfo), .packages= c("chipmine", "XLConnect", "org.Anidulans.FGSCA4.eg.db")) %do% {
  
  tfFuncDir <- paste(TF_dataPath, "/", tfInfo$sampleId[i], "/", "functional_analysis", sep = "")
  if(! dir.exists(tfFuncDir)){
    dir.create(tfFuncDir)
  }
  
  tfOutPrefix <- paste(tfFuncDir, "/", tfInfo$sampleId[i], sep = "")
  tfExcelOut <- paste(tfOutPrefix, "_functional_analysis.xlsx", sep = "")
  
  ## select the genes which has peak near ATG
  peakData <- get_TF_binding_data(genesDf = geneSet, exptInfo = tfInfo[i,]) %>% 
    dplyr::filter_at(.vars = vars(starts_with("hasPeak.")),
                     .vars_predicate = all_vars(. == TRUE))
  
  ## store data in excel
  unlink(tfExcelOut, recursive = FALSE, force = FALSE)
  exc <- loadWorkbook(tfExcelOut , create = TRUE)
  xlcFreeMemory()
  
  if(nrow(peakData) > 0){
    ## topGO enrichment for genes which has macs2 peak near ATG
    goEnrich <- topGO_enrichment(goMapFile = file_topGoMap, genesOfInterest = peakData$gene, goNodeSize = 5)

    readr::write_tsv(x = goEnrich, path = paste(tfOutPrefix, "_macs2_topGO.tab", sep = ""), col_names = T)
    wrkSheet = "macs2_topGO"
    createSheet(exc, name = wrkSheet)
    createFreezePane(exc, sheet = wrkSheet, 2, 2)
    setMissingValue(object = exc, value = "NA")
    writeWorksheet(object = exc, data = goEnrich, sheet = wrkSheet, header = T)
    setAutoFilter(object = exc, sheet = wrkSheet, reference = aref(topLeft = "A1", dimension = dim(goEnrich)))
    
    
    ## pathway enrichment for genes which has macs2 peak near ATG
    keggEnr <- keggprofile_enrichment(genes = peakData$gene,
                                      orgdb = orgDb, keytype = "GID",
                                      keggOrg = keggOrgCode, pvalCut = 0.05)
    
    readr::write_tsv(x = keggEnr, path = paste(tfOutPrefix, "_macs2_keggprofile.tab", sep = ""), col_names = T)
    wrkSheet = "macs2_kegg"
    createSheet(exc, name = wrkSheet)
    createFreezePane(exc, sheet = wrkSheet, 2, 2)
    setMissingValue(object = exc, value = "NA")
    writeWorksheet(object = exc, data = keggEnr, sheet = wrkSheet, header = T)
    setAutoFilter(object = exc, sheet = wrkSheet, reference = aref(topLeft = "A1", dimension = dim(keggEnr)))
    
    
    
  }
  
  xlcFreeMemory()
  saveWorkbook(exc)
  
  
}








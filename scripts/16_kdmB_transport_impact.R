library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(scales)
library(ggplot2)


## As per Ozgur's input, kdmB-del strain use less nutrients
## this script checks the transport BP gene's fold change in kdmB_del_20h / 20h_polII


rm(list = ls())

source(file = "E:/Chris_UM/GitHub/omics_util/GO_enrichment/topGO_functions.R")

##################################################################################

analysisName <- "kdmb_del_transport"
outPrefix <- here::here("kdmB_analysis", "kdmB_del_20h", analysisName)

polII1 <- "An_untagged_20h_polII_1"
polII2 <- "An_untagged_48h_polII_1"
otherPolII <- c("An_kdmB_del_20h_polII_1", "An_kdmB_del_48h_polII_1")

# "deeptools", "miao", "normalizedmatrix", "normalizedmatrix_5kb"
matrixType <- "normalizedmatrix_5kb"
matrixDim = c(500, 200, 100, 10)

## polII signal fold change pairs
polIIDiffPairs <- list(
  p1 = list(
    name = "untagged_48h_vs_untagged_20h_polII",
    title = "polII log2(untagged_48h \n vs untagged_20h)",
    samples = c(polII2, polII1)
  ),
  p2 = list(
    name = "kdmB_del_48h_vs_kdmB_del_20h_polII",
    title = "polII log2(kdmB_del_48h \n vs kdmB_del_20h)",
    samples = otherPolII[c(2, 1)]
  ),
  p3 = list(
    name = "kdmB_del_20h_vs_untagged_20h_polII",
    title = "polII log2(kdmB_del_20h \n vs untagged_20h_polII)",
    samples = c(otherPolII[1], polII1)
  ),
  p4 = list(
    name = "kdmB_del_48h_vs_untagged_48h_polII",
    title = "polII log2(kdmB_del_48h \n vs untagged_48h_polII)",
    samples = c(otherPolII[2], polII2)
  )
)


sampleList <- c(polII1, polII2, otherPolII)


## genes to read
file_exptInfo <- here::here("data", "referenceData/sampleInfo.txt")
file_genes <- here::here("data", "referenceData/AN_genesForPolII.bed")
file_topGoMap <- "E:/Chris_UM/Database/A_Nidulans/ANidulans_OrgDb/geneid2go.ANidulans.topGO.map"
file_geneInfo <- "E:/Chris_UM/Database/A_Nidulans/A_nidulans_FGSC_A4_geneClasses.txt"


TF_dataPath <- here::here("data", "TF_data")
polII_dataPath <- here::here("data", "polII_data")
hist_dataPath <- here::here("data", "histone_data")


orgDb <- org.Anidulans.FGSCA4.eg.db

##################################################################################

## genes to read
geneSet <- data.table::fread(file = file_genes, header = F,
                             col.names = c("chr", "start", "end", "geneId", "score", "strand")) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$geneId, columns = "DESCRIPTION", keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("geneId" = "GID"))

## gene information annotations: cluster and TF and polII expression values
geneInfo <- add_gene_info(file = file_geneInfo, clusterDf = geneSet)

head(geneInfo)

##################################################################################

polIIData <- get_sample_information(exptInfoFile = file_exptInfo,
                                    samples = c(polII1, polII2, otherPolII),
                                    dataPath = polII_dataPath)


exptData <- dplyr::bind_rows(polIIData)

polII_ids <- exptData$sampleId[which(exptData$IP_tag == "polII")]


polIICols <- list(
  exp = structure(polII_ids, names = polII_ids),
  is_expressed = structure(paste("is_expressed", ".", polII_ids, sep = ""), names = polII_ids)
)


expressionData <- get_polII_expressions(exptInfo = exptData,
                                        genesDf = geneInfo)

## add fold change columns
for (i in names(polIIDiffPairs)) {
  expressionData <- get_fold_change(df = expressionData,
                                    nmt = polIIDiffPairs[[i]]$samples[1],
                                    dmt = polIIDiffPairs[[i]]$samples[2],
                                    newCol = polIIDiffPairs[[i]]$name,
                                    isExpressedCols = polIICols$is_expressed)
}


kdmB20h_diff <- dplyr::mutate(
  expressionData,
  group = if_else(
    condition = !! as.name(polIIDiffPairs$p3$name) > 0.5,
    true = "up",
    false = if_else(
      condition = !! as.name(polIIDiffPairs$p3$name) < -0.5,
      true = "down",
      false = "noDEG"
    ))) %>% 
  dplyr::select(geneId, !!polIIDiffPairs$p3$name, !!!polIIDiffPairs$p3$samples, group) %>% 
  dplyr::filter(abs(!! as.name(polIIDiffPairs$p3$name)) > 1)


# keytypes(x = orgDb)
# AnnotationDbi::select(x = orgDb, keys = "GO:0007005", columns = c("ONTOLOGYALL", "GID", "DESCRIPTION"), keytype = "GOALL")


## topGO enrichment
goEnrich <- dplyr::group_by(kdmB20h_diff, group) %>%
  do(topGO_enrichment(goMapFile = file_topGoMap, genes = .$geneId, goNodeSize = 5))

gt <- c("GO:0046323", "GO:0006810", "GO:0008645", "GO:0019740")


gtMap <- dplyr::group_by(kdmB20h_diff, group) %>%
  do(GO_map(genes = .$geneId, goTerms = gt, org = orgDb)) %>% 
  dplyr::ungroup()

transportGenes <- gtMap %>% 
  dplyr::mutate(geneId = strsplit(x = genes, split = ";", fixed = T)) %>% 
  tidyr::unnest() %>% 
  dplyr::select(geneId, GOID, TERM) %>% 
  dplyr::left_join(y = geneDesc, by = c("geneId" = "GID")) %>% 
  dplyr::left_join(y = kdmB20h_diff, by = "geneId")

readr::write_tsv(x = transportGenes, path = paste(outPrefix, ".genes.tab", sep = ""))






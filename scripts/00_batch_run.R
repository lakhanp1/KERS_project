###########################################################################
#######################
### DESeq2 pairwise ###
#######################

## RNAseq DESeq2 differential gene expression batches
file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_conf.txt")
diffDataPath <- here::here("analysis", "02_polII_DEGs")
script_deseq2 <- here::here("scripts", "05_DESeq2_diff.R")

runConfig <- suppressMessages(readr::read_tsv(file_RNAseq_info))

rowId <- 1

for (rowId in 1:nrow(runConfig)) {
  cat("Processing DESeq2:", runConfig$comparison[rowId], "...\n")

  system2(
    command = "Rscript",
    args = c(script_deseq2, "--config", file_RNAseq_info, "--deg", runConfig$comparison[rowId]),
    stdout = paste("logs/stdouterr.DESeq2.", runConfig$comparison[rowId],".log", sep = ""),
    stderr = paste("logs/stdouterr.DESeq2.", runConfig$comparison[rowId],".log", sep = "")
  )

}

###########################################################################
####################################
### DESeq2 Functional enrichment ###
####################################

# file_productionData <- here::here("data", "reference_data", "production_data.summary.tab")
# diffDataPath <- here::here("analysis", "06_polII_diff")
# script_deseq2_GO <- here::here("scripts", "a_explore_polII", "a09_polII_DEGs_functional_enrichment.R")
# file_RNAseq_info <- here::here("data", "reference_data", "polII_DESeq2_DEG_info.txt")
# 
# 
# productionData <- suppressMessages(readr::read_tsv(file = file_productionData)) %>% 
#   dplyr::filter(has_polII_ChIP == "has_data")
# 
# 
# rowId <- 1
# 
# for (rowId in 1:nrow(productionData)) {
#   
#   system2(
#     command = "Rscript",
#     # args = c(script_deseq2_GO, "--help"),
#     args = c(script_deseq2_GO, "--config", file_RNAseq_info, "--deg", productionData$degId[rowId]),
#     stdout = paste("logs/stdouterr.functional_enrichment.", productionData$degId[rowId],".log", sep = ""),
#     stderr = paste("logs/stdouterr.functional_enrichment.", productionData$degId[rowId],".log", sep = "")
#   )
#   
#   print(productionData$degId[rowId])
# }


###########################################################################



# 
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_hepA_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_hepA_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_hepA_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_hepA_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_laeA_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_laeA_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_laeA_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_laeA_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_laeA_hepA_GFP_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_laeA_hepA_GFP_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_laeA_hepA_GFP_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_laeA_hepA_GFP_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_laeA_kdmB_del_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_laeA_kdmB_del_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_cclA_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_cclA_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_cclA_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_cclA_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_cclA_kdmA_del_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_cclA_kdmA_del_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_cclA_kdmA_del_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_cclA_kdmA_del_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecmB_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecmB_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecmB_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecmB_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecmB_kdmA_del_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecmB_kdmA_del_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecmB_kdmA_del_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecmB_kdmA_del_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmA_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmA_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmA_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmA_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmA_short_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmA_short_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmA_short_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmA_short_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_mcmA_20h_MYC_1  An_untagged_20h_MYC_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_mcmA_20h_MYC_2  An_untagged_20h_MYC_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_mcmA_48h_MYC_1  An_untagged_48h_MYC_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_mcmA_48h_MYC_2  An_untagged_48h_MYC_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_mcmA_kdmA_del_20h_MYC_1  An_untagged_20h_MYC_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_mcmA_kdmA_del_20h_MYC_2  An_untagged_20h_MYC_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_mcmA_kdmA_del_48h_MYC_1  An_untagged_48h_MYC_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_mcmA_kdmA_del_48h_MYC_2  An_untagged_48h_MYC_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rstB_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rstB_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rstB_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rstB_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rstB_kdmA_del_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rstB_kdmA_del_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rstB_kdmA_del_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rstB_kdmA_del_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecoA_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecoA_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecoA_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecoA_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecoA_kdmB_del_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecoA_kdmB_del_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecoA_kdmB_del_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecoA_kdmB_del_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecoA_sntB_del_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecoA_sntB_del_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecoA_sntB_del_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_ecoA_sntB_del_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmB_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmB_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmB_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmB_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmB_laeA_del_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmB_laeA_del_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmB_laeA_del_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmB_laeA_del_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmB_sntB_del_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmB_sntB_del_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmB_sntB_del_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_kdmB_sntB_del_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_laeA_kdmB_del_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_laeA_kdmB_del_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rpdA_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rpdA_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rpdA_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rpdA_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rpdA_kdmB_del_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rpdA_kdmB_del_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rpdA_kdmB_del_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rpdA_kdmB_del_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rpdA_sntB_del_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rpdA_sntB_del_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rpdA_sntB_del_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_rpdA_sntB_del_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sntB_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sntB_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sntB_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sntB_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sntB_kdmB_del_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sntB_kdmB_del_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sntB_kdmB_del_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sntB_kdmB_del_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sudA_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sudA_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sudA_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sudA_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sudA_ecoA_teton_20h_HA_1  An_untagged_20h_HA_1  An_untagged_dox_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sudA_ecoA_teton_20h_HA_2  An_untagged_20h_HA_2  An_untagged_dox_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sudA_ecoA_teton_48h_HA_1  An_untagged_48h_HA_1  An_untagged_dox_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sudA_ecoA_teton_48h_HA_2  An_untagged_48h_HA_2  An_untagged_dox_48h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sudA_kdmB_del_20h_HA_1  An_untagged_20h_HA_1  An_untagged_20h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sudA_kdmB_del_20h_HA_2  An_untagged_20h_HA_2  An_untagged_20h_polII_2")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sudA_kdmB_del_48h_HA_1  An_untagged_48h_HA_1  An_untagged_48h_polII_1")
# system("Rscript E:/Chris_UM/Codes/Mevlut_data_analysis/TF_polII_args.R An_sudA_kdmB_del_48h_HA_2  An_untagged_48h_HA_2  An_untagged_48h_polII_2")

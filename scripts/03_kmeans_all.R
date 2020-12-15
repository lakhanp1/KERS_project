suppressPackageStartupMessages(library(chipmine))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))

## this script performs k-means clustering on (-2kb)-TSS-TES-(+1kb) profile matrix
## generated for all genes in each TF ChIPseq data

rm(list = ls())

##################################################################################

file_exptInfo <- here::here("data", "reference_data", "sampleInfo.txt")
file_genes <- here::here("data", "reference_data", "AN_genesForPolII.bed")

TF_dataPath <- here::here("..", "data", "A_nidulans", "TF_data")
polII_dataPath <- here::here("..", "data", "A_nidulans", "polII_data")
hist_dataPath <- here::here("..", "data", "A_nidulans", "histone_data")
other_dataPath <- here::here("..", "data", "A_nidulans", "other_data")

file_tf_macs2 <- paste(TF_dataPath, "/", "sample_tf_macs2.list", sep = "")
file_tf <- paste(TF_dataPath, "/", "sample_tf.list", sep = "")


set.seed(20)
doClustering <- TRUE
clusters <- 7
tfYlim <- 0.996              ##0.999

geneFilter <- c("AN5245", "AN3245")

cl <- makeCluster(4) #not to overload your computer
registerDoParallel(cl)


##################################################################################

# geneSet <- data.table::fread(file = file_genes, header = F,
#                              col.names = c("chr", "start", "end", "name", "score", "strand")) %>% 
#   dplyr::mutate(length = end - start) %>% 
#   dplyr::filter(! name %in% geneFilter)

geneSet <- suppressMessages(
  readr::read_tsv(
    file = file_genes,col_names = c("chr", "start", "end", "geneId", "score", "strand")
  )) %>%
  dplyr::mutate(length = end - start) %>% 
  dplyr::filter(! geneId %in% geneFilter)

tfSampleList <- suppressMessages(
  readr::read_tsv(file = file_tf_macs2, col_names = c("id"),  comment = "#")
)

# tfSampleList <- data.frame(id = c("An_ecoA_20h_HA_1", "An_ecoA_48h_HA_1", "An_kdmB_20h_HA_1", "An_kdmB_48h_HA_1", "An_rpdA_20h_HA_1", "An_rpdA_48h_HA_1", "An_sntB_20h_HA_1", "An_sntB_48h_HA_1", "An_kdmB_20h_HA_2", "An_kdmB_48h_HA_2", "An_rpdA_20h_HA_2", "An_rpdA_48h_HA_2", "An_sntB_20h_HA_2", "An_sntB_48h_HA_2", "An_ecoA_kdmB_del_20h_HA_1", "An_ecoA_kdmB_del_48h_HA_1", "An_rpdA_kdmB_del_20h_HA_1", "An_rpdA_kdmB_del_48h_HA_1", "An_sntB_kdmB_del_20h_HA_1", "An_sntB_kdmB_del_48h_HA_1", "An_ecoA_20h_HA_2", "An_ecoA_48h_HA_2", "An_ecoA_kdmB_del_20h_HA_2", "An_ecoA_kdmB_del_48h_HA_2", "An_rpdA_kdmB_del_20h_HA_2", "An_rpdA_kdmB_del_48h_HA_2", "An_sntB_kdmB_del_20h_HA_2", "An_sntB_kdmB_del_48h_HA_2", "An_ecoA_sntB_del_20h_HA_2", "An_ecoA_sntB_del_48h_HA_2", "An_kdmB_laeA_del_20h_HA_1", "An_kdmB_laeA_del_48h_HA_1", "An_laeA_kdmB_del_20h_HA_1", "An_laeA_kdmB_del_48h_HA_1", "An_sudA_kdmB_del_20h_HA_1", "An_sudA_kdmB_del_48h_HA_1"))

tf_info <- get_sample_information(
  exptInfoFile = file_exptInfo,
  samples = tfSampleList$id,
  dataPath = TF_dataPath,
  profileMatrixSuffix = "normalized_profile"
)

# i <- 1

foreach(i = 1:nrow(tf_info),
        .packages = c("chipmine")) %dopar% {
          
          ## read the profile matrix
          mat1 <- chipmine::import_profile_from_file(
            file = tf_info$matFile[i],
            source = "deeptools",
            signalName = tf_info$sampleId[i],
            selectGenes = geneSet$geneId)
          
          ## check the distribution in data
          quantile(mat1, c(seq(0, 0.9, by = 0.1), 0.95, 0.99, 0.992, 0.995, 0.999, 0.9999, 1), na.rm = T)
          # col_fun <- colorRamp2(quantile(mat1, c(0.50, 0.995), na.rm = T), c("white", "red"))
          
          km <- chipmine::profile_matrix_kmeans(
            mat = mat1,
            km = clusters,
            clustFile = "temp.kmeans.tab",
            # clustFile = tf_info$clusterFile[i],
            name = tf_info$sampleId[i])
          
          
          cat(as.character(Sys.time()), "Done...", tf_info$sampleId[i], "\n\n")
          
        }


stopCluster(cl)








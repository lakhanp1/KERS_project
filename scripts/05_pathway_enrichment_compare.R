library(chipmine)
library(org.Anidulans.FGSCA4.eg.db)
library(here)

## 1) compare the keggprofile KEGG enrichment result between pairs of samples
## 2) generate comparison bar chart

rm(list = ls())

##################################################################################
analysisName <- "kdmA_complex_KEGG_cmp"
outPrefix <- here::here("kdmA_analysis/functional_analysis", analysisName)


tfs20h <- c("An_kdmA_20h_HA_1", "An_cclA_20h_HA_1", "An_mcmA_20h_MYC_1", "An_ecmB_20h_HA_1", "An_rstB_20h_HA_1")
tfs48h <- c("An_kdmA_48h_HA_1", "An_cclA_48h_HA_1", "An_mcmA_48h_MYC_1", "An_ecmB_48h_HA_1", "An_rstB_48h_HA_1")
tfNames <- c("kdmA", "cclA", "mcmA", "ecmB", "rstB")

tfPairs <- purrr::pmap(.l = list(tfs20h, tfs48h, tfNames),
                       .f = function(s1, s2, n){
                         pr <- list(s1 = s1, s2 = s2, name = n)
                         return(pr)
                       }) %>% 
  purrr::set_names(paste("p", 1:length(tfNames), sep = ""))



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
                             col.names = c("chr", "start", "end", "gene", "score", "strand")) %>% 
  dplyr::mutate(length = end - start)

geneDesc <- AnnotationDbi::select(x = orgDb, keys = geneSet$gene, columns = c("DESCRIPTION"), keytype = "GID")

geneSet <- dplyr::left_join(x = geneSet, y = geneDesc, by = c("gene" = "GID"))


##################################################################################

tfInfo <- get_sample_information(exptInfoFile = file_exptInfo,
                                 samples = c(tfs20h, tfs48h),
                                 dataPath = TF_dataPath,
                                 matrixSource = "normalizedmatrix")



s1 <- tfPairs$p5$s1
s2 <- tfPairs$p5$s2
cmpName <- tfPairs$p5$name

## read KEGG enrichment file and store
f1 <- paste(TF_dataPath, "/", s1, "/functional_analysis/", s1, "_macs2_keggprofile.tab", sep = "")
f2 <- paste(TF_dataPath, "/", s2, "/functional_analysis/", s2, "_macs2_keggprofile.tab", sep = "")

s1Rename <- structure(c("Gene_Found", "pvalue", "geneIds"),
                      names = paste("s1", c("Gene_Found", "pvalue", "geneIds"), sep = "."))

s2Rename <- structure(c("Gene_Found", "pvalue", "geneIds"),
                      names = paste("s2", c("Gene_Found", "pvalue", "geneIds"), sep = "."))

newColNames <- list()
newColNames[[s1]] <- structure(names(s1Rename), names = unname(s1Rename))
newColNames[[s2]] <- structure(names(s2Rename), names = unname(s2Rename))


## prepare kegg1 list
kegg1 <- suppressMessages(readr::read_tsv(file = f1, col_names = TRUE)) %>% 
  dplyr::filter(Gene_Pathway >= 5) %>% 
  dplyr::select(-Percentage, -pvalueAdj) %>% 
  dplyr::rename(!!! s1Rename) %>% 
  dplyr::mutate(sample1 = s1)


## prepare kegg2 list
kegg2 <- suppressMessages(readr::read_tsv(file = f2, col_names = TRUE)) %>% 
  dplyr::filter(Gene_Pathway >= 5) %>% 
  dplyr::select(-Percentage, -pvalueAdj) %>% 
  dplyr::rename(!!! s2Rename) %>% 
  dplyr::mutate(sample2 = s2)


  # dplyr::filter(Gene_Pathway >= 5) %>% 
  # purrr::transpose() %>% 
  # purrr::set_names(nm = map(., "pathway_id")) %>% 
  # purrr::map(.f = function(a){
  #   a[["genes"]] <- unlist(strsplit(x = a[["geneIds"]], split = "/", fixed = TRUE))
  #   return(a)
  # })

## merged pair data
keggData <- dplyr::full_join(x = kegg1, y = kegg2, by = c("pathway_id", "Pathway_Name", "Gene_Pathway")) %>% 
  dplyr::mutate_at(.vars = vars(ends_with(".Gene_Found")),
                   .funs = funs(if_else(condition = is.na(.), true = 0, false = as.numeric(.)))) %>% 
  tidyr::replace_na(list(sample1 = s1, sample2 = s2))

## get the overlap statistics
overlapDf <- purrr::transpose(keggData) %>% 
  purrr::set_names(nm = purrr::map(., "pathway_id")) %>% 
  purrr::map_dfr(.f = function(x){
    dt <- list()
    
    dt <- x[c("pathway_id", "Pathway_Name", "Gene_Pathway", "sample1", "sample2")]
    # dt[["pathway_id"]] <- x[["pathway_id"]]
    # dt[["Pathway_Name"]] <- x[["Pathway_Name"]]
    # dt[["Gene_Pathway"]] <- x[["Gene_Pathway"]]
    
    s1Genes <- NA
    s2Genes <- NA
    
    if(! is.na(x[["s1.geneIds"]])){
      s1Genes <- unlist(strsplit(x = x[["s1.geneIds"]], split = "/", fixed = TRUE))
    }
    
    if(! is.na(x[["s2.geneIds"]])){
      s2Genes <- unlist(strsplit(x = x[["s2.geneIds"]], split = "/", fixed = TRUE))
    }
    
    
    dt[["ns1"]] <- length(na.exclude(s1Genes))
    dt[["ns2"]] <- length(na.exclude(s2Genes))
    dt[["total"]] <- length(na.exclude(union(s1Genes, s2Genes)))
    
    ## s1-s2 common
    if(dt[["total"]] == 0){
      dt[["nCommon"]] <- NA
      dt[["s1Specific"]] <- NA
      dt[["s2Specific"]] <- NA
    } else{
      dt[["nCommon"]] <-  length(na.exclude(intersect(s1Genes, s2Genes)))
      dt[["s1Specific"]] <- length(na.exclude(setdiff(s1Genes, s2Genes)))
      dt[["s2Specific"]] <- length(na.exclude(setdiff(s2Genes, s1Genes)))
    }
    
    
    return(dt)
  })



pltDf <- tidyr::gather(overlapDf, key = "group", value = "count",
                       s1Specific, nCommon, s2Specific,
                       factor_key = FALSE)


pltDf$group <- factor(pltDf$group, levels = c("s2Specific", "nCommon", "s1Specific"))

title <- paste("KEGG enrichment:", s1, "vs", s2)

pt <- ggplot(data = pltDf, mapping = aes(x = Pathway_Name, y = count)) +
  geom_bar(mapping = aes(fill = group),
           stat = 'identity', position = "fill", color = "black") +
  geom_text(mapping = aes(label = count),
            stat='identity', position=position_fill(vjust=0.5)) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(
    name = "",
    breaks = c("s1Specific", "nCommon", "s2Specific"),
    labels = c(s1, "common", s2),
    values = c("#ff4d4d", "#e6ccff", "#0080ff")
  ) +
  ggtitle(str_wrap(title, 80)) +
  ylab("Relative %") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(1, 4, 1, 1), "cm"))

# png(filename = paste(outPrefix, "_", cmpName, ".png", sep = ""), width = 4000, height = 4000, res = 400)
pdf(file = paste(outPrefix, "_", cmpName, ".pdf", sep = ""), width = 10, height = 10)
pt
dev.off()



# KdmB-EcoA-RpdA-SntB (KERS) chromatin complex ChIPseq data analysis

KERS is a tetrameric nucleosome remodelling complex in the filamentous fungus Aspergillus nidulans that is composed of the H3K4 histone demethylase KdmB, a cohesin acetyltransferase (EcoA), a histone deacetylase (RpdA) and a reader protein (SntB). This repository has the scripts related to KERS ChIPseq data analysis.

### Required tools for data analysis
#### *A. nidulans* AnnotationHub packages like `org.db`, `TxDb` and `BSgenome`
```R
## org.db
devtools::install_github(
  "lakhanp1/fungal_resources/A_nidulans/org.Anidulans.FGSCA4.eg.db",
  upgrade = "never")
  
## TxDb
devtools::install_github(
  "lakhanp1/fungal_resources/A_nidulans/TxDb.Anidulans.FGSCA4.AspGD.GFF",
  upgrade = "never")
  
## BSgenome
devtools::install_github(
  "lakhanp1/fungal_resources/A_nidulans/BSgenome.Anidulans.FGSCA4.AspGD",
  upgrade = "never")

```

#### R package `chipmine` with utility functions for ChIPseq data management and processing
`chipmine` is under active development. For the analysis involving current manuscript, `chimpine` from branch `KERS_paper` was used.
``` R
devtools::install_github(repo = "lakhanp1/chipmine", ref = "KERS_paper")
```

### Data processing
All the ChIPseq data was processed as per steps mentioned in [fungal ChIPseq data processing workflow of Koon Ho Wong lab](https://github.com/lakhanp1/bioinformatics_notes/blob/master/data/ChIPseq/01_CL_ChIPseq_pipeline.md).
<br>

### Data Availability
**NCBI SRA data:** ChIPseq raw data is deposited to NCBI SRA Bioproject [PRJNA591094](https://dataview.ncbi.nlm.nih.gov/object/PRJNA591094?reviewer=2vodlkj4hlsqn2bd1p37bp3723).

**UCSC genome browser tracks:** All the normalized bigWig tracks for ChIPseq data are made available on UCSC genome browser and can be accessed using link <https://genome-asia.ucsc.edu/s/lakhanp1/A.%20nidulans%20KERS>

<br><br><br>

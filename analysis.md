## bowtie2 mapping
```bash
bowtie2 -p 6 --local  -x <bt2-idx> -1 <> -2 <> | samtools view -bS - | samtools sort  -O bam -o <>.bam
```

## Index and alignment stats
```bash
for i in `cat sample_ANidulans.list`
do
cd $i
bash /home/lakhanp/scripts/ChIPseq_scripts/template_ChIPseq_config.sh -c $REF_DIR/reference_genomes.yaml -o A_nidulans --polII >> generalJob.sh
sed "s/SAMPLE_ID/$i/g" /home/lakhanp/scripts/ChIPseq_scripts/template_ChIPseq_process.sh >> generalJob.sh
cd ..
done
```


## Calculate GC bias
```bash
for i in `cat sample_ANidulans.list`
do
cd $i
printf "##Compute GC bias
computeGCBias -b %s_bt2.bam --effectiveGenomeSize 29850950 -g $REF_DIR/A_nidulans_FGSC_A4/reference/A_nidulans_FGSC_A4_version_s10-m04-r03_chromosomes.2bit -l 60 --GCbiasFrequenciesFile %s_GCbiasFreq.txt --biasPlot %s_GCbias.png 
error_exit \$? \n\n" $i $i $i >> generalJob.sh
cd ..
done
```

##Pol-II raw read count
```bash
for i in `cat sample_polII.list`
do
cd $i
printf "##Pol-II raw read count
bedtools multicov -bams %s_bt2.bam -bed /home/lakhanp/database/A_nidulans_FGSC_A4/annotation/A_nidulans_FGSC_A4_version_s10-m04-r03_CDS_Unique.bed > %s.CDS.readCount.tab
error_exit \$? \n\n" $i $i >> generalJob.sh
printf "##Pol-II raw read count for intergenic region
bedtools multicov -bams %s_bt2.bam -bed /home/lakhanp/database/A_nidulans_FGSC_A4/annotation/intergenic_bins.bed > %s.intergenicBin.readCount.tab
error_exit \$? \n\n" $i $i >> generalJob.sh
cd ..
done
```


## Print control information and peak type (narrow/broad) information
```bash
## IMP
## use the printf information from the excel file
for i in `cat sample_tf_macs2.list`
do cd $i
printf "peakType=\'narrow\'\n\n" >> generalJob.sh
cd ..
done
```


## MACS2 peak calling: use template
```bash
for i in `cat sample_tf_macs2.list`
do cd $i
sed "s/SAMPLE_ID/$i/g" /home/lakhanp/scripts/ChIPseq_scripts/template_ChIPseq_macs2.sh >> generalJob.sh
cd ..
done
```



## Generate profile matrix
```bash
for i in `cat sample_ANidulans.list`
do
cd $i
printf "##generate profile matris (-2kb == normalized(geneBody) to 2kb == +1kb)\n"
computeMatrix scale-regions -S %s_normalized.bw  -R $REF_DIR/A_nidulans_FGSC_A4/annotation/A_nidulans_FGSC_A4_version_s10-m04-r03_CDS_Unique.bed -m 2000 -b 2000 -a 1000 --numberOfProcessors 2  --outFileName %s_normalized_profile.tab.gz
error_exit $?
cd ..
done
```


##   Motif enrichment analysis  ####

### meme-chip analysis: with control
```bash
while IFS=$'\t' read -r name fasta neg
do
	printf "## meme-chip: ${name}\n"
	printf "## meme-chip: ${name}
	meme-chip -order 1 -meme-minw 5 -meme-maxw 30 -meme-nmotifs 5 -meme-mod zoops \
	-db /home/lakhanp/tools/meme_motif_databases_12.19/JASPAR/JASPAR2018_CORE_fungi_non-redundant.meme \
	-desc ${name} -oc ${name} -neg ${neg} ${fasta} < /dev/null
	"
	printf "## done...\n\n"
done < memechip_de_conf.tab
```

### meme-chip analysis: without control
```bash
while IFS=$'\t' read -r name fasta
do
printf "## meme-chip: ${name}\n"
dm meme-chip -order 1 -meme-nmotifs 3 -meme-mod anr -meme-p 8 \
-db JASPAR2018_CORE_fungi_non-redundant.meme \
-desc ${name} -oc ${name} ${fasta}
#"
printf "## done...\n\n"
done < memechip_conf.tab


## run meme-chip in parallel using GNU parallel
parallel --keep-order --joblog $PWD/parallel_job.logs --halt now,fail=1 -j 6 -a memechip_conf.tab --colsep '\t' \
meme-chip -order 1 -meme-nmotifs 3 -meme-mod anr -meme-p 8 \
-db $TOOLS_PATH/motif_databases/JASPAR/JASPAR2020_CORE_fungi_non-redundant.meme \
-desc {1} -oc {1} {2}
```


## copy data to local
```bash
for i in `cat sample_ANidulans.list`
do 
cd $i
cp ${i}_normalized.bw ${i}_normalized_profile.tab.gz ${i}_polii_expr.tab.rel.mat ../localCopy/
cd ..
done

for i in `cat tf_AN.list`; do cd $i; 
cp macs2_*/${i}*{.narrowPeak,.broadPeak,.tab} ../localCopy/
cd ..
done
```
## concat all raw fastqs from different lanes:
cd /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/
mkdir concat_raw_fastq

cd raw_data/
for I in *L001*; do 
 base=`cut -d'_' -f1 <<< $I`;
 printf "$I\t$base\n";
 cat $base* > ../concat_raw_fastq/$base.fastq.gz
done

### run pipeline
### run_AQUAS.pl script take the sample file, identifies the correct sample & input and runs the Kundaje pipeline
module load java/1.8.0_latest
perl run_AQUAS.pl
perl run_AQUAS.no_input.pl

### couple samples failed the pipeline; try to rerun again manually:
## NN161
python /well/got2d/agata/bin/TF_chipseq_pipeline/chipseq.py TF idr -species hg19 -peak_caller macs2 -fastq1 concat_raw_fastq/CHE324.fastq.gz -ctl_fastq1 concat_raw_fastq/CHE323.fastq.gz --out-dir aquas/NN161_rerun --title NN161

## NN078 - try this one without the pseudo-reps, maybe this is the problem?
python /well/got2d/agata/bin/TF_chipseq_pipeline/chipseq.py TF idr -species hg19 -peak_caller macs2 -fastq1 concat_raw_fastq/CHO232.fastq.gz -ctl_fastq1 concat_raw_fastq/CHO231.fastq.gz --out-dir aquas/NN078_rerun --title NN078 -no_pseudo_rep

## NN078 - try this one without the pseudo-reps, maybe this is the problem?
#### try changing the IDR threshold ?
python /well/got2d/agata/bin/TF_chipseq_pipeline/chipseq.py TF idr -species hg19 -peak_caller macs2 -fastq1 concat_raw_fastq/CHO232.fastq.gz -ctl_fastq1 concat_raw_fastq/CHO231.fastq.gz --out-dir aquas/NN078 --title NN078 -idr_thresh 0.1




#######  for some reason the AQUAS pipeline is still not happy to produce the blacklist filtered peak lists for these two samples, I will filter out the blacklist regions with bedtools for all samples 

#/apps/well/bedtools/2.24.0/intersectBed -v -a <(zcat -f NN078/peak/macs2/rep1/CHO232.nodup.tagAlign_x_CHO231.nodup.tagAlign.narrowPeak.gz) -b <(zcat /well/mccarthy/production/chip-seq/resources/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz ) | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' | grep -P 'chr[\dXY]+[ \t]' | gzip -nc > NN078.filt.narrowPeak.gz

mkdir filt_narrowPeak
for I in NN*; do
	printf "$I\n";
	/apps/well/bedtools/2.24.0/intersectBed -v -a <(zcat -f $I/peak/macs2/rep1/*.nodup.tagAlign.narrowPeak.gz) -b <(zcat /well/mccarthy/production/chip-seq/resources/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz ) | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' | grep -P 'chr[\dXY]+[ \t]' | gzip -nc > filt_narrowPeak/$I.filt.narrowPeak.gz
done

#### calculate ... values using SPP
/apps/well/R/3.2.2/bin/Rscript /well/mccarthy/production/chip-seq/code/run_spp.R -c=NN064/align/rep1/CHO244.nodup.tagAlign.gz -odir=NN064_spp

mkdir spp_stats
for I in NN*; do
	printf "$I\n"
	tag=`ls $I/align/rep1/*.nodup.tagAlign.gz`
	printf "$tag\n"
	echo "/apps/well/R/3.2.2/bin/Rscript /well/mccarthy/production/chip-seq/code/run_spp.R -c=$tag > spp_stats/$I.spp.txt" > spp_stats/spp_stats.$I.sh
	qsub -V -N $I -cwd -q short.qc -P mccarthy.prjc spp_stats/spp_stats.$I.sh
done

# Minimum cross-correlation value 0.3069481
# Minimum cross-correlation shift 1500
# Top 3 cross-correlation values 0.310868950406642,0.309511464628273,0.308871078129451
# Top 3 estimates for fragment length 110,175,215
# Window half size 280
# Phantom peak location 80
# Phantom peak Correlation 0.3103518
# Normalized Strand cross-correlation coefficient (NSC) 1.012774
# Relative Strand cross-correlation Coefficient (RSC) 1.151932
# Phantom Peak Quality Tag 1

#### try calling peaks with spp - for NN064:
/apps/well/R/3.2.2/bin/Rscript /well/mccarthy/production/chip-seq/code/run_spp.R -c=NN064/align/rep1/CHO244.nodup.tagAlign.gz -i=NN064/align/ctl1/CHO243.nodup.tagAlign.gz -savp=NN064.pdf -savn=NN064.spp.narrowPeak


#### set up manually the IDR for replicated samples:
## HNF4A: D14 pancreatic prog 
python /well/got2d/agata/bin/TF_chipseq_pipeline/chipseq.py TF idr -species hg19 -peak_caller macs2 -fastq1 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHO238.fastq.gz -fastq2 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHC1366.fastq.gz -ctl_fastq1 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHO237.fastq.gz -ctl_fastq2 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHC1365.fastq.gz --out-dir /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/aquas_IDR/HNF4A_D14_PP --title HNF4A_D14_PP

## HNF4A: D20 endocrine prog
python /well/got2d/agata/bin/TF_chipseq_pipeline/chipseq.py TF idr -species hg19 -peak_caller macs2 -fastq1 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHO240.fastq.gz -fastq2 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHC1370.fastq.gz -ctl_fastq1 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHO239.fastq.gz -ctl_fastq2 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHC1369.fastq.gz --out-dir /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/aquas_IDR/HNF4A_D20_EP --title HNF4A_D20_EP

## D20 - 1 of the IDR calculations failed - RERUN!!!

## HNF4A: D35 beta-like cells
python /well/got2d/agata/bin/TF_chipseq_pipeline/chipseq.py TF idr -species hg19 -peak_caller macs2 -fastq1 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHO242.fastq.gz -fastq2 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHC1372.fastq.gz -ctl_fastq1 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHO241.fastq.gz -ctl_fastq2 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHC1371.fastq.gz --out-dir /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/aquas_IDR/HNF4A_D35_BLC --title HNF4A_D35_BLC

## D35 - IDR failed -- RERUN!

## HNF4A: D8 hepatoblasts
python /well/got2d/agata/bin/TF_chipseq_pipeline/chipseq.py TF idr -species hg19 -peak_caller macs2 -fastq1 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/SHH002.fastq.gz -fastq2 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHC1368.fastq.gz -ctl_fastq1 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/SHH001.fastq.gz -ctl_fastq2 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHC1367.fastq.gz --out-dir /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/aquas_IDR/HNF4A_D8_HEP --title HNF4A_D8_HEP

## HNF1A: EndoC-bH1
python /well/got2d/agata/bin/TF_chipseq_pipeline/chipseq.py TF idr -species hg19 -peak_caller macs2 -fastq1 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHO234.fastq.gz -fastq2 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHE320.fastq.gz -ctl_fastq1 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHO233.fastq.gz -ctl_fastq2 /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/concat_raw_fastq/CHE319.fastq.gz --out-dir /well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore/aquas_IDR/HNF1A_EndoC_bH1 --title HNF1A_EndoC_bH1


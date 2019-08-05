#### for samples with replicates - see overlap of peaks with bedtools intersect

mkdir /well/mccarthy/users/agata/StemBANCC/HNF4_1_ChIP_Singapore/biological_reps
cd /well/mccarthy/users/agata/StemBANCC/HNF4_1_ChIP_Singapore/biological_reps/

### EndoC-bH1_HNF1A: NN064 and NN092
/apps/well/bedtools/2.24.0/intersectBed -a <(zcat -f ../aquas_flt_narrowPeak/NN064.filt.narrowPeak.gz) -b <(zcat -f ../aquas_flt_narrowPeak/NN092.filt.narrowPeak.gz) > EndoC-bH1_HNF1A.intersection.bed

### HB_D8_HNF4A: NN117 and NN159_D8
/apps/well/bedtools/2.24.0/intersectBed -a <(zcat -f ../aquas_flt_narrowPeak/NN117.filt.narrowPeak.gz) -b <(zcat -f ../aquas_flt_narrowPeak/NN159_D8.filt.narrowPeak.gz) > HB_D8_HNF4A.intersection.bed

### PP_D14_HNF4A: NN158_D14 and NN119
/apps/well/bedtools/2.24.0/intersectBed -a <(zcat -f ../aquas_flt_narrowPeak/NN158_D14.filt.narrowPeak.gz) -b <(zcat -f ../aquas_flt_narrowPeak/NN119.filt.narrowPeak.gz) > PP_D14_HNF4A.intersection.bed

wc -l *
#    11 EndoC-bH1_HNF1A.intersection.bed
#  2715 HB_D8_HNF4A.intersection.bed
#   818 PP_D14_HNF4A.intersection.bed

#### for samples where Natasha shared results - check overlap in called peaks
cd /well/mccarthy/users/agata/StemBANCC/HNF4_1_ChIP_Singapore/annotated_peaks_Singapore/
for I in *peaks.txt; do
 awk '{if (NR>1){print $2 "\t" $3 "\t" $4}}' $I > $I.bed
done

## NN161
/apps/well/bedtools/2.24.0/intersectBed -a D35_NN161_HNF1A.annotated_63peaks.txt.bed -b <(zcat -f ../aquas_flt_narrowPeak/NN161.filt.narrowPeak.gz) > NN161.intersection.bed
## NN092
/apps/well/bedtools/2.24.0/intersectBed -a EndoC_NN092_HNF1A.annotated_101peaks.txt.bed -b <(zcat -f ../aquas_flt_narrowPeak/NN092.filt.narrowPeak.gz) > NN092.intersection.bed
## NN177
/apps/well/bedtools/2.24.0/intersectBed -a HepG2_NN177_HNF1A.annotated_5687peaks.txt.bed -b <(zcat -f ../aquas_flt_narrowPeak/NN177.filt.narrowPeak.gz) > NN177.intersection.bed

wc -l *intersection.bed
#     9 NN092.intersection.bed
#    47 NN161.intersection.bed
#  5544 NN177.intersection.bed


#### run HOMER to find motifs enriched in peaks
mkdir HOMER_motifs
for I in aquas_flt_narrowPeak/NN0*; do
base=`basename $I .filt.narrowPeak.gz`
findMotifsGenome.pl <(zcat -f $I) hg19 HOMER_motifs/$base -size given &
done

## also on the biological reps
cd /well/mccarthy/users/agata/StemBANCC/HNF4_1_ChIP_Singapore/biological_reps/
mkdir HOMER_motifs
for I in *.bed; do
base=`basename $I .intersection.bed`
findMotifsGenome.pl $I hg19 HOMER_motifs/$base -size given &
done

### and on the IDR peaks
cd /well/mccarthy/users/agata/StemBANCC/HNF4_1_ChIP_Singapore/aquas_IDR_narrowPeak/
mkdir HOMER_motifs
for I in NN*; do
base=`basename $I .rep1_pr.IDR0.05.filt.narrowPeak.gz`
findMotifsGenome.pl <(zcat -f $I) hg19 HOMER_motifs/$base -size given &
done


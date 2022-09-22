#!/bin/bash

echo "MACS2 Peak Calling"
mkdir MACS2_outputs
mkdir MACS2_outputs/broad
mkdir MACS2_outputs/narrow
macs2 callpeak -t BAM_files/final/H3K9ac-APOE3.final.bam -c BAM_files/final/input-mus-uglia-APOE3.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K9ac-APOE3 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K9ac-APOE3_macs2.log
macs2 callpeak -t BAM_files/final/H3K9ac-APOE4.final.bam -c BAM_files/final/input-mus-uglia-APOE4.final.bam -f BAMPE --mfold 5 50 -p 0.001 -g mm -n H3K9ac-APOE4 --outdir MACS2_outputs/narrow | tee MACS2_outputs/narrow/H3K9ac-APOE4_macs2.log


#Making Big WIG Files
mkdir BIGWIG_files
bamCoverage --bam BAM_files/final/H3K9ac-APOE3.final.bam -o BIGWIG_files/H3K9ac-APOE3.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --ignoreForNormalization chrX --extendReads --numberOfProcessors 2 --skipNonCoveredRegions
bamCoverage --bam BAM_files/final/H3K9ac-APOE4.final.bam -o BIGWIG_files/H3K9ac-APOE4.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --ignoreForNormalization chrX --extendReads --numberOfProcessors 2 --skipNonCoveredRegions


#Intersecting Bed Files
multiIntersectBed -i  Peak_files/narrow/H3K9ac-APOE3_peaks.narrowPeak Peak_files/narrow/H3K9ac-APOE4_peaks.narrowPeak > master_peak_H3K9ac.bed
bedtools merge -i master_peak_H3K9ac.bed -d 100 -c 4 -o count > master_peak_merged100_H3K9ac.bed


#Extracting Counts from Peak-file
bedtools multicov -bams BAM_files/final/H3K9ac-APOE3.final.bam BAM_files/final/H3K9ac-APOE4.final.bam BAM_files/final/H3K27ac-APOE3.final.bam BAM_files/final/H3K27ac-APOE4.final.bam -bed master_peak_merged100_all.bed > extracted_data


#Annotation of Peaks 
annotatePeaks.pl Annotation_file_H3K9ac.txt mm10   > annotated_DPs_H3K9ac

echo 'Detected following samples'
find ./fastq -name "*.fastq.gz" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d'_' -f1 | sort -u > sample_list.txt
cat sample_list.txt

echo 'Unzipping files'
gunzip fastq/*.fastq.gz

echo 'Checking quality of untrimmed reads'
mkdir 'QC/untrimmed'
fastqc fastq/*.fastq -t 4 -o 'QC/untrimmed'
multiqc QC/untrimmed -o QC/untrimmed

echo 'Trimming reads - stand by'
mkdir fastq/trimmed_reads
cat sample_list.txt | while read sample; do
	echo $sample
	cutadapt --cores 0 --minimum-length 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o fastq/trimmed_reads/${sample}.1.trimmed.fastq.gz -p fastq/trimmed_reads/${sample}.2.trimmed.fastq.gz fastq/${sample}_R1_001.fastq fastq/${sample}_R2_001.fastq
done

echo 'Zipping untrimmed files'
gzip fastq/*.fastq

echo 'Checking quality of trimmed_reads'
mkdir 'QC/trimmed'
fastqc fastq/trimmed_reads/*.fastq.gz -t 4 -o 'QC/trimmed'
multiqc QC/trimmed -o QC/trimmed 

echo 'Performing Bowtie2 alignment to SAM files'
mkdir SAM_files
cat sample_list.txt | while read sample; do
	bowtie2 -p 2 -q --local -x ~/Desktop/mm10/mm10 -1 fastq/trimmed_reads/${sample}.1.trimmed.fastq.gz -2 fastq/trimmed_reads/${sample}.2.trimmed.fastq.gz -S SAM_files/${sample}.unsorted.sam
done

echo "Creating BAM files from SAM files"
mkdir BAM_files
mkdir BAM_files/unsorted
cat sample_list.txt | while read sample; do
	samtools view -h -S -b -o BAM_files/unsorted/${sample}.unsorted.bam SAM_files/${sample}.unsorted.sam
done

find ./unsorted -name "*.unsorted.bam" -maxdepth 1 -type f -exec basename "{}" \; |  cut -d'_' -f1 | sort -u > sample_list.txt
cat sample_list.txt

echo "Sorting BAM files"
mkdir BAM_files/sorted
cat sample_list.txt | while read sample; do
	sambamba sort -t 4 -o BAM_files/sorted/${sample}.sorted.bam BAM_files/unsorted/${sample}
done

echo "Filtering sorted BAM files"
mkdir BAM_files/final
cat sample_list.txt | while read sample; do
	sambamba view -h -t 4 -f bam -F "[XS] == null and not unmapped and not duplicate"  BAM_files/sorted/${sample}.sorted.bam > BAM_files/final/${sample}.final.bam
done

echo "Indexing final BAM files"
cat sample_list.txt | while read sample; do
	samtools index BAM_files/final/${sample}.final.bam
done

echo "Alignment statistics"
cat sample_list.txt | while read sample; do
	samtools flagstat BAM_files/sorted/${sample}.sorted.bam
done

echo "Making BIGWIG_files"
mkdir BIGWIG_files
bamCoverage --bam BAM_files/final/'SAMPLE_NAME'.final.bam -o BIGWIG_files/'SAMPLE_NAME'.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --ignoreForNormalization chrX --extendReads --numberOfProcessors 10 --skipNonCoveredRegions


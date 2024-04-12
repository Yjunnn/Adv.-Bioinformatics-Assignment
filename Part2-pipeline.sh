
#!/bin/bash

# Data Preparation

mkdir pipeline # make directory for this pipeline

cd ~/pipeline # change directory
mkdir data meta results logs scripts # make subdirectory

cd ~/pipeline/data
mkdir trimmed_fastq untrimmed_fastq

# download the data and rename it
wget -O NGS0001.R1.fastq.gz https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget -O NGS0001.R2.fastq.gz https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget -O annotation.bed https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

# move the data to untrimmed file
mv *fastq.gz ~/pipeline/data/untrimmed_fastq

# download the reference
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

# use fastqc to assess quality on raw data and move the results to results file
cd ~/pipeline/data/untrimmed_fastq
fastqc -t 4 *.fastq.gz
mkdir ~/pipeline/results/fastqc_untrimmed_reads
mv ~/pipeline/data/untrimmed_fastq/*fastqc* ~/pipeline/results/fastqc_untrimmed_reads/


# Trimming
trimmomatic PE  \
  -threads 4 \
  -phred33 \
  /home/ubuntu/pipeline/data/untrimmed_fastq/NGS0001.R1.fastq.gz /home/ubuntu/pipeline/data/untrimmed_fastq/NGS0001.R2.fastq.gz \
  -baseout ~/pipeline/data/trimmed_fastq/trimmed_data \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50

# fastqc on paired trimmed fastq data and save the results to results file
fastqc -t 4 /home/ubuntu/pipeline/data/trimmed_fastq/trimmed_data_1P \
	/home/ubuntu/pipeline/data/trimmed_fastq/trimmed_data_2P
mkdir ~/pipeline/results/fastqc_trimmed_reads
mv ~/pipeline/data/trimmed_fastq/*fastqc* ~/pipeline/results/fastqc_trimmed_reads/


# Alignment

# index the reference genome first
mkdir -p ~/pipeline/data/reference
mv ~/pipeline/data/hg19.fa.gz ~/pipeline/data/reference/
bwa index ~/pipeline/data/reference/hg19.fa.gz

# alignment with BWA on trimmed data
mkdir ~/pipeline/data/aligned_data

# use group info in this dataset
bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1.111.D1375ACXX\tSM:NGS0001\tPL:ILLUMINA' -I 250,50  ~/pipeline/data/reference/hg19.fa.gz ~/pipeline/data/trimmed_fastq/trimmed_data_1P ~/pipeline/data/trimmed_fastq/trimmed_data_2P > ~/pipeline/data/aligned_data/NGS0001.sam

# convert sam to bam and sort it
cd ~/pipeline/data/aligned_data
samtools view -h -b NGS0001.sam > NGS0001.bam  
samtools sort NGS0001.bam > NGS0001_sorted.bam 
samtools index NGS0001_sorted.bam #This will generate a .bai index file


# Basic Alignment post processing

# mark duplicates
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
samtools index NGS0001_sorted_marked.bam

# filter BAM file
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam

# index filtered bam_
samtools index NGS0001_sorted_filtered.bam

# generate standard alignment statistics (i.e. flagstats, idxstats, depth of coverage, insert size) 
samtools flagstat NGS0001_sorted_filtered.bam
samtools idxstats NGS0001_sorted_filtered.bam
bedtools coverage -a annotation.bed -b NGS0001_sorted_filtered.bam > coverage.txt
picard CollectInsertSizeMetrics I=NGS0001_sorted_filtered.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf


# Basic Variant calling

# index the reference
zcat ~/pipeline/data/reference/hg19.fa.gz > ~/pipeline/data/reference/hg19.fa
samtools faidx ~/pipeline/data/reference/hg19.fa

# call variants with freebayes
freebayes --bam ~/pipeline/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/pipeline/data/reference/hg19.fa --vcf ~/pipeline/results/NGS0001.vcf --bed annotation.bed

bgzip ~/pipeline/results/NGS0001.vcf
tabix -p vcf ~/pipeline/results/NGS0001.vcf.gz

# filter the VCF
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
	~/pipeline/results/NGS0001.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf

# using bedtools to filter the vcf file for the annotation regions:
bedtools intersect -header -wa -a ~/pipeline/results/NGS0001_filtered.vcf -b ../annotation.bed \
 	> ~/pipeline/results/NGS0001_filtered_annotation.vcf
bgzip ~/pipeline/results/NGS0001_filtered_annotation.vcf
tabix -p vcf ~/pipeline/results/NGS0001_filtered_annotation.vcf.gz


# Variant annotation

# convert VCF to Annovar input format
 ./convert2annovar.pl -format vcf4 ~/pipeline/results/NGS0001_filtered_annotation.vcf.gz > ~/pipeline/results/NGS0001_filtered_annotation.avinput

# run Annovar table function and get the csv output
 ./table_annovar.pl ~/pipeline/results/NGS0001_filtered_annotation.avinput humandb/ -buildver hg19  \
   -out ~/pipeline/results/NGS0001_filtered_annotation -remove   \
      -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout






























































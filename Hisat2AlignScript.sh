#!/bin/bash

## This script is to do the read alignments using hisat2

## Inputs: processed data directory

## STEP 1: Read the directory

echo "Script starts from $PWD";
dir=$1
filelist=$(ls $dir | grep fastq.gz);
echo "Directory read: $filelist"

## STEP 2: create a file with the names of each pair of samples

cd $dir
filenames=$(ls $filelist | sed 's/_q30l50_R1.fastq.gz//' | sed 's/_q30l50_R2.fastq.gz//' | uniq)
echo "filenames: $filenames"
cd ..
pwd 

mkdir Results/mapping_bam
mkdir Results/mapping_bam/mappingsummaries

## STEP 3: To do the alignments

for file in $filenames
	do
	echo "$file about to be processed"
	hisat2 -p 16 --dta -x Genome/GCA_008087245.1_ASM808724v1_genomic -1 ./00_processed_data/$file"_q30l50_R1_fastq.gz" -2 ./00_processed_data/$file"_q30l50_R2_fastq.gz" --summary-file $file".hisat2.mappingsummary.txt" | samtools view -F 4 -Sb -o $file".hisat2.bam" -
	mv $file".hisat2.mappingsummary.txt" Results/mapping_bam/mappingsummaries/ 
	mv $file".hisat2.bam" Results/mapping_bam
	echo "alignment of $file done"
	done

echo "all alignments done"

## STEP 4: Get the number of mapped reads with samtools stats

mkdir Results/mapping_bam/samstats

for file in $filenames
	do
	echo "Doing samstats of $file"
	samtools stats Results/mapping_bam/$file".hisat2.bam" > $file".hisat2.samstats.txt"
	mv $file".hisat2.samstats.txt" Results/mapping_bam/samstats
	echo "Samstats of $file done"
	done

cd Results/mapping_bam/samstats
rep "reads mapped:" *.hisat2.samstats.txt >> number_reads_mapped.txt
echo "file number_reads_mapped.txt created"

## STEP 5: Sort bam files

cd Results
for file in $filenames
	do
	samtools sort -@ 8 -o mapping_bam/$file".hisat2.bam" unsorted_mapping_bam/$file".hisat2.bam"
	done


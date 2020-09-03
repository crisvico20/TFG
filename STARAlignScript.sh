#!/bin/bash

## This script is to do the STAR Alignment. The input directory is the reads files directory

## STEP 1: GENERATE GENOME INDEXES

## 1.1 Read the directory

echo "Script starts from $PWD";
dir=$1
filelist=$(ls $dir | grep fastq.gz);
echo "Directory read: $filelist"

## 1.2 Create a file with the names of each pair of samples

cd $dir
filenames=$(ls $filelist | sed 's/_q30l50_R1_fastq.gz//' | sed 's/_q30l50_R2_fastq.gz//' | uniq)
echo "filenames: $filenames"
cd ..

## 1.3 Create the genomeDir

mkdir Results/STAR/index/

## 1.4 Run STAR in the genomeGenerate mode

pwd
echo "Start to generate genome indexes"
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir Results/STAR/index/ --genomeFastaFiles MyNewAvocado_Reference.fasta --sjdbGTFfile Genome/MyNewAvocado_Reference.gff --sjdbGTFtagExonParentTranscript Parent
echo "Generation of genome indexes done"


## STEP 2: MAP READS TO THE GENOME

## 2.1 Create a results directory for STAR

mkdir Results/STAR/Files

## 2.2 Run STAR in the alignReads mode
pwd
for file in $filenames; 
	do
	echo "Starting alignment of $file"
	STAR --runThreadN 8 --runMode alignReads --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --genomeDir Results/STAR/index/ --readFilesIn 00_processed_data/$file"_q30l50_R1_fastq.gz" 00_processed_data/$file"_q30l50_R2_fastq.gz" --outFileNamePrefix Results/STAR/Files/$file
	echo "Alignment of $file done"
	done

echo "STAR Alignment finished"

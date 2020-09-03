#!/bin/bash

## This script is  to process the files and count the reads and the bp per base of each one
## Arguments: directory with raw files

## STEP 1: Fastq-stats of each file of raw data
## 1.1 Read the files in the directory

echo "Script starts from $PWD";
dir=$1
filelist=$(ls $dir | grep fastq.gz);
echo "Directory read: $filelist";

echo "Changing to dir $dir";
cd $dir

## 1.2 Create the directory for the stats
mkdir Raw_stats

## 1.3 Make a for loop that goes to each file and does the stats, saving the stats in the raw files stats directory
echo "Loop starts at the dir $PWD";
for file in $filelist
	do
	echo "Interacting file $file";
	fastq-stats $file > $file.stats.txt
	mv $file.stats.txt Raw_stats/
	echo "File stats of $file done";
	done

cd Raw_stats
pwd
ls
grep "reads" *.stats.txt >> numberreads_raw.txt
echo "file numberreads_raw.txt created"
cd ..

## STEP 2: Process reads removing adapters and the low-quality (<30) nucleotides in the extreme of the reads
## 2.1 Create a list with the names of each pair of samples

filenames=$(ls $filelist | sed 's/_R1.fastq.gz//' | sed 's/_R2.fastq.gz//' | uniq)
echo "uniq file name of samples"
echo $filenames

## 2.2 Create a directory for processed files

cd ..
pwd
mkdir processed_data


## 2.3 Make a for loop that goes to each file and execute fastq-mcf on it

for file in $filenames
	do
	echo "$file about to be processed"
	fastq-mcf -q 30 -l 50 -o $file"_q30l50_R1_fastq.gz" -o $file"_q30l50_R2_fastq.gz" 00_sources/IlluminaAdapters_V2.fasta ./$dir/$file"_R1.fastq.gz" ./$dir/$file"_R2.fastq.gz"
	mv $file"_q30l50_R1_fastq.gz" processed_data
	mv $file"_q30l50_R2_fastq.gz" processed_data
	echo "$file processed and moved"
	done 

## STEP 3: Fastq-stats of processed files

## 3.1 Create the directory for the stats

mkdir processed_data/Processed_stats

## 3.2 Make a for loop that goes to each file and does the stats, saving the stats in the stats directory for raw files

Processed_files=$(ls processed_data | grep fastq.gz);
echo "Directory read: $Processed_files";
cd processed_data

for file in $Processed_files
        do
        fastq-stats $file > $file.stats.txt
        mv $file.stats.txt Processed_stats
        echo "stats of $file done"
        done

cd Processed_stats
grep "reads" *.stats.txt >> numberreads_processed.txt
echo "file numberreads_processed.txt created"
cd ..


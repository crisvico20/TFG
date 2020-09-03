#!/bin/bash

## This script is to ensamble transcripts with stringtie

## Input: sorted bam files

## STEP 1:  read the input directory

echo "Script starts from $PWD";
dir=$1
filelist=$(ls $dir | grep .out.bam);
echo "Directory read: $filelist"

cd $dir
filenames=$(ls $filelist | sed 's/Aligned.sortedByCoord.out.bam//' | uniq)
echo "filenames: $filenames"

mkdir ../gtf_files


## STEP 2: Stringtie

for file in $filenames
	do
	stringtie -p 8 -G ../../../Genome/Persea_americana_annotation.gtf -o ../gtf_files/$file".gtf" -l $file $file"Aligned.sortedByCoord.out.bam"
	echo "Stringtie of $file done"
	done

cd ..

## STEP 3: Manually, create mergelist.txt and do stringtie --merge to obtain a gtf file with all the transcripts

## STEP 4: Estimate transcript abundances and create table counts for Ballgown

mkdir Ballgown

for file in $filenames
	do
	mkdir Ballgown/$file
	stringtie -e -B -p 8 -G gtf_files/stringtie_merged.gtf -o Ballgown/$file/$file".gtf" Files/$file"Aligned.sortedByCoord.out.bam"
	echo "stringtie of $file done"
	done

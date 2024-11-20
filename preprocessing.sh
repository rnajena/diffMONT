#!/bin/bash

## created by Daria Meyer @11/2024
## Bash script to created sorted bedmethyl file out of multiple bedmethyl files as input

## Input: bedmethyl files, 11th column has to contain methylation values (either 0-100 or 0-1)
## Output: merged and sorted bedmethyl file

## Loop through the given bedmethyl files and merge them
i=1;
for bedfile in "$@" 
do
    echo "Filename - $i: $bedfile";
	new_name=${bedfile##*/}
	new_name="${new_name%%_*}_$i"
	echo "New Filename - $i: $new_name"

	head $bedfile | awk '{if ($11 != "nan"){print $0}}' | sed "s/$/\t$new_name/" >> mergedMethylation.bed
	# awk '{if ($11 != "nan"){print $0}}' $bedfile | sed "s/$/\t$new_name/" 

	
	## increase pointer
    i=$((i + 1));
done

## Sort the bedmethyl file
sort -k 1,1 -V -s -k2,2 -o mergedMethylation.bed mergedMethylation.bed

#!/usr/bin/env python

"""
author: Daria Meyer
created: 12.04.2023

Contact:
daria.meyer@uni-jena.de

This script extracts regions which are differentially methylated between two conditions.

Input: 
* merged and sorted bedmethyl samples with data for all samples (merged by chromosome and position)

The following steps are made:
* For each CpG dinucleotide: 
*	- remove CpG dinucleotide if methylated > 10% in any control (Filter1)
* 	- calculate boxplot stats for remaining CpG dinucleotides
* 	- remove CpG dinucleotide if not median tumor sample methylation above upper whisker of control samples
* Find MSP primers with 18-24 nt length
* Find MSP region = combination of MSP primers, spanning a region of length 60-400 nt (with a sliding window approach)


Example:
	python diffONT.py --bedmethylFile $data/bedmethyl/sorted.bed --tumors t0044c t0085c t0126c --controls t0099n t0025n t0045n --outfolder $data/results
	python diffONT.py --bedmethylFile $data/bedmethyl/sorted.bed --tumors t0044c t0085c t0126c --controls t0099n t0025n t0045n --reference hg38.fa --annotation hg38.gtf --outfolder $data/results

"""

import numpy as np
import pandas as pd
import os
import argparse
import math
from collections import Counter

# Instantiate the parser
parser = argparse.ArgumentParser(description='Get differentially methylated regions')

## Define input
parser.add_argument('--bedmethylFile', required=True, 
					help=' merged and sorted bedmethyl file')
parser.add_argument('--controls', required=True, nargs='+', 
					help='list of control sample names')
parser.add_argument('--tumors', required=True, nargs='+', 
					help='list of tumor sample names')
## Define output
parser.add_argument('--outfolder', required=True,  
                    help='outfolder in which the results will be stored')

## Define maximum control methylation per sample per position
parser.add_argument('--maxMethControl', default = 10,
                    help='maximum average methylation allowed for the control sample')
## Define minimum amount of differently methylated cytosines  in primer
parser.add_argument('--minCpGs', default = 3,
                    help='minimum amount of differentially methylated cytosines in primer')
## Define primer length
parser.add_argument('--minPrimerLength', default = 18,
                    help='minimum length required for primers')
parser.add_argument('--maxPrimerLength', default = 24,
                    help='maximum length required for primers')
## Define amplicon length
parser.add_argument('--minAmpliconLength', default = 60,
                    help='minimum length required for MSP regions')
parser.add_argument('--maxAmpliconLength', default = 400,
                    help='maximum length required for MSP regions')

## Additional optional arguments
parser.add_argument("--annotation", help='Ensembl annotation file -- will be screened for genes overlapping MSP regions')
parser.add_argument("--reference", help='Genome reference file -- will be used to calculate genomic stats like GC-content')
parser.add_argument("--boxplotData", help='path to boxplotData.txt -- if given, this will be used instead of newly calculating boxplot data')
parser.add_argument("--primerData", help='path to primerData.txt -- if given, this will be used instead of newly calculating primer data')


## FUNCTIONS ##
###############################################################################
def importFasta(genomeFile):
	""" Read the given reference genome into a dictionary """
	genome = {}
	sequence = ""
	with open(genomeFile, 'r') as f:
		for line in f.readlines():
			if (line.startswith(">")):
				if sequence:
					genome[header] = sequence
					# break ## only for chromosome 1
				header = line[1:].split()[0].strip()
				sequence = ""
			else:
				sequence = sequence + line.strip()
		genome[header] = sequence
		return genome


def getFasta(genome, chro, strand, start, end):
	""" Extract the DNA sequence for a specific region from the given reference genome """
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	## extract the sequence of interest
	if(chro in genome.keys()):
		## Check if index out of bounds, e.g. when extracting PCR product +200nt
		if (int(end) >= len(genome[chro])):
			print(f"{end} out of bounds for {chro}. Taking last base of chro instead. This is ok, if extracting PCR product DNA (+-200nt).")
			end = len(genome[chro])-1
		seq = genome[chro][int(start):int(end)]
		if (strand == "-"):
			seq = "".join(complement.get(base, base) for base in reversed(seq))
	else:
		seq = ""
	return seq


def getRegions(bedMethyFile, maxMethControl, controls, tumors):
	"""
	Function to remove uninteresting CpG dinucleotides from a list of bedmethyl file entries.
	remove:
		* CpGs which are methylated in control samples above maxMethControl
		* CpGs which median methylation of tumor samples is below upper whisker methylation of control samples

	Args:
		bedMethyFile (str): 	Path to the merged and sorted bedmethyl file.
		maxMethControl (int): 	Maximum percentage of methylation allowed per control sample.
		controls (list): 		List of control sample names
		tumors (list): 			List of tumor sample names

	Returns:
		list: boxplot stats for each CpG dinucleotide which fulfills requirements
		list: updated bedmethyl file (prob. written into new file)
	"""
	resList = []
	tmp = []
	boxplotData = [['chro', 'pos', 'strand', 'median_ctrl', 'upper_quartile_ctrl', 'lower_quartile_ctrl', 'upper_whisker_ctrl', 'lower_whisker_ctrl', 'median_tmr', 'upper_quartile_tmr', 'lower_quartile_tmr', 'upper_whisker_tmr', 'lower_whisker_tmr', 'control_coverage', 'control_methylation', 'control_samples', 'tumor_coverage', 'tumor_methylation', 'tumor_samples']]
	keep = True

	## Write the header to boxplotData.tsv -- overwrite existing files
	with open(os.path.join(outFolder, 'boxplotData.tsv'), 'w') as outstream:
		outstream.write(f"chro\tpos\tstrand\tmedian_ctrl\tupper_quartile_ctrl\tlower_quartile_ctrl\tupper_whisker_ctrl\tlower_whisker_ctrl\tmedian_tmr\tupper_quartile_tmr\tlower_quartile_tmr\tupper_whisker_tmr\tlower_whisker_tmr\tcontrol_coverage\tcontrol_methylation\tcontrol_samples\ttumor_coverage\ttumor_methylation\ttumor_samples\n")

	## Start output stream
	with open(os.path.join(outFolder, 'boxplotData.tsv'), 'a') as outstream:

		## Initialization for the first line. Is there a better way?
		with open(bedMethyFile) as file:
			line = file.readline()
			old_chro = line.split()[0]
			old_pos = int(line.split()[1])


		## read in the merged bedmethyl line by line
		with open(bedMethyFile) as file:
			for line in file:

				chro = line.split()[0]
				pos = int(line.split()[1])
				strand = line.split()[5]
				coverage = int(line.split()[9])
				methylation = float(line.split()[10])
				sample = line.split()[len(line.split())-1]   ## There might be some additional columns or not. Sample name is always the last column

				if ((chro == old_chro) and (pos == old_pos)):
					## Same position, next sample. Add data from all samples into one list and process together
					tmp.append([chro, pos, strand, coverage, methylation, sample])
					## Set keep to false if methylated in any of the controls above threshold
					if sample in controls:
						if methylation > maxMethControl:
							keep = False
				else: 
					## If not methylated the control samples
					if keep:
						## Calculate stats for previous position (= data from tmp)
						data = pd.DataFrame(tmp, columns = ['chro', 'pos', 'strand', 'coverage', 'methylation', 'sample'])
						# print(data.info())

						tmp_boxplot = calcBoxplotStats(data, controls, tumors)
						tumor_median = tmp_boxplot[8]
						control_upper_whisker = tmp_boxplot[6]

						## Check if median from tumor sample is higher than max from control
						if (tumor_median > control_upper_whisker): ## False if one of the values is nan.
							## If all criteria are fulfilled, add data to final lists boxplotData and resList
							boxplotData.append(tmp_boxplot) 
							resList.extend(tmp)
							for value in tmp_boxplot:
								outstream.write(f"{value}\t")
							outstream.write(f"\n")
							outstream.flush()

					## Always initialize tmp, old_chro, old_pos, keep 
					tmp = [[chro, pos, strand, coverage, methylation, sample]]
					old_chro = chro
					old_pos = pos
					keep = True
					## Set keep to false if methylated in any of the controls above threshold
					if sample in controls:
						if methylation > maxMethControl:
							keep = False

		## Write the CpG dinucleotides of interest into a new file
		with open(os.path.join(outFolder, 'interestingCpGs.tsv'), 'w') as f:
			for line in resList:
				for value in line:
					f.write(f"{value}\t")
				f.write(f"\n")

		## Make pandas dataframe out of list
		# dtype_dict = {'chro': 'str', 'pos': 'int64', 'strand': 'str', 'median_ctrl' : 'float64', 'upper_quartile_ctrl': 'float64', 'lower_quartile_ctrl': 'float64', 'upper_whisker_ctrl': 'float64', 'lower_whisker_ctrl': 'float64', 'median_tmr': 'float64', 'upper_quartile_tmr': 'float64', 'lower_quartile_tmr': 'float64', 'upper_whisker_tmr': 'float64', 'lower_whisker_tmr': 'float64', 'control_coverage' : 'str', 'control_methylation' : 'str', 'control_samples' : 'str', 'tumor_coverage' : 'str', 'tumor_methylation' : 'str', 'tumor_samples' : 'str'}
		cols = ['chro', 'pos', 'strand', 'median_ctrl', 'upper_quartile_ctrl', 'lower_quartile_ctrl', 'upper_whisker_ctrl', 'lower_whisker_ctrl', 'median_tmr', 'upper_quartile_tmr', 'lower_quartile_tmr', 'upper_whisker_tmr', 'lower_whisker_tmr', 'control_coverage', 'control_methylation', 'control_samples', 'tumor_coverage', 'tumor_methylation', 'tumor_samples']
		boxplotData_df = pd.DataFrame(boxplotData, columns = cols)

	return(boxplotData_df)


def calcBoxplotStats(data, controls, tumors):
	""" Function to calculate boxplot statistics which are:
		* upper and lower whisker (= 0th and 100th percentile)
		* upper and lower quartile (= 25th and 75th quartile)
		* median (= 50th quartile)
	"""
	
	## Extract only control samples into pandas dataframe
	df_control = data.loc[data['sample'].isin(controls)]
	# print(df_control)
	control_coverage = ';'.join(str(x) for x in df_control['coverage'].tolist())
	# print("coverage:")
	# print(control_coverage)
	control_methylation = ';'.join(str(x) for x in df_control['methylation'].tolist())
	# print("methylation:")
	# print(control_methylation)
	control_samples = ';'.join(str(x) for x in df_control['sample'].tolist())

	## Calculate boxplot statistics for control samples
	control_median = np.nanmedian(df_control['methylation'])
	control_upper_quartile = np.nanpercentile(df_control['methylation'], 75)
	control_lower_quartile = np.nanpercentile(df_control['methylation'], 25)
	control_iqr = control_upper_quartile - control_lower_quartile
	control_upper_whisker = df_control['methylation'][df_control['methylation']<=control_upper_quartile+1.5*control_iqr].max()
	control_lower_whisker = df_control['methylation'][df_control['methylation']>=control_lower_quartile-1.5*control_iqr].min()
	stats_control = [control_median, control_upper_quartile, control_lower_quartile, control_upper_whisker, control_lower_whisker]

	## Extract only tumor samples into pandas dataframe
	df_tumor = data.loc[data['sample'].isin(tumors)]
	# print(df_tumor)
	
	tumor_coverage = ';'.join(str(x) for x in df_tumor['coverage'].tolist())
	tumor_methylation = ';'.join(str(x) for x in df_tumor['methylation'].tolist())
	tumor_samples = ';'.join(str(x) for x in df_tumor['sample'].tolist())

	## Calculate boxplot statistics for tumor samples
	tumor_median = np.nanmedian(df_tumor['methylation'])
	tumor_upper_quartile = np.nanpercentile(df_tumor['methylation'], 75)
	tumor_lower_quartile = np.nanpercentile(df_tumor['methylation'], 25)
	tumor_iqr = tumor_upper_quartile - tumor_lower_quartile
	tumor_upper_whisker = df_tumor['methylation'][df_tumor['methylation']<=tumor_upper_quartile+1.5*tumor_iqr].max()
	tumor_lower_whisker = df_tumor['methylation'][df_tumor['methylation']>=tumor_lower_quartile-1.5*tumor_iqr].min()
	stats_tumor = [tumor_median, tumor_upper_quartile, tumor_lower_quartile, tumor_upper_whisker, tumor_lower_whisker]
	stats = [data.loc[0]['chro'], data.loc[0]['pos'], data.loc[0]['strand']] + stats_control + stats_tumor + [control_coverage] + [control_methylation] + [control_samples] + [tumor_coverage] + [tumor_methylation] + [tumor_samples]

	return stats


def getPimerPointer(positions, maxPrimerLength, minCpGs):
	""" Function to find all possible max. length (=24nt) primer regions with min. 3 CpGs """
	primer_indices = []
	# positions = positions from boxplotdata for one chromosome and one strand
	for i in range(len(positions)-2):
		## Starting primer contains exactly three CpGs (aka minCpGs)
		pointer1 = i
		pointer2 = pointer1 + (minCpGs - 1)  ## 2 CpGs weiter 
		## check if minimal primer (3 diff. CpGs) has minimal required length (24 nucleotides)
		if(positions[pointer2]-positions[pointer1] <= maxPrimerLength):
			## extend the primer region as much as possible up to the maximal length <= 24nt
			while((pointer2 < len(positions)-1) and ((positions[pointer2] - positions[pointer1]) <= maxPrimerLength)):
				pointer2_old = pointer2
				pointer2 += 1 
			## Check if the primer exceeds the minimal required length
			if((positions[pointer2] - positions[pointer1]) <= maxPrimerLength):
				primer_indices.append([pointer1, pointer2])
			## If so, use the previous position.
			else:
				primer_indices.append([pointer1, pointer2_old])

	return primer_indices


def calcPrimerScore(boxplotData, primer, norms):
	""" 
	Calculate a score for a given primer region.
	Returns the given primer extended by the calculated score
	"""
	# print(boxplotData.head())
	# print(boxplotData.info())
	score = 0
	# print('\t'.join(str(x) for x in primer))
	# print(boxplotData.loc[(boxplotData['chro'] == primer[0]) & (boxplotData['strand'] == primer[1])].iloc[primer[2]:primer[3]+1])
	data = boxplotData.loc[(boxplotData['chro'] == primer[0]) & (boxplotData['strand'] == primer[1]) & (boxplotData['pos'] >= primer[2]) & (boxplotData['pos'] <= primer[3])]
	# print(data)
	sample = data['tumor_samples'].str.split(';').tolist()
	# print(sample)

	## Calculate average coverage over all positions in tumor samples
	coverage_tumor = data['tumor_coverage'].str.split(';').tolist()
	coverage_tumor_flat = [x for xs in coverage_tumor for x in xs]
	coverage_tumor_average = np.mean(list(map(int, coverage_tumor_flat)))

	## Calculate average coverage over all positions in control samples
	coverage_control = data['control_coverage'].str.split(';').tolist()
	coverage_control_flat = [x for xs in coverage_control for x in xs]
	coverage_control_average = np.mean(list(map(int, coverage_control_flat)))

	## Calculate average methylation over all positions in control samples
	methylation_tumor = data['tumor_methylation'].str.split(';').tolist()
	methylation_tumor_flat = [x for xs in methylation_tumor for x in xs]
	methylation_tumor_average = np.nanmean(list(map(float, methylation_tumor_flat)))

	## Calculate average methylation over all positions in control samples
	methylation_control = data['control_methylation'].str.split(';').tolist()
	methylation_control_flat = [x for xs in methylation_control for x in xs]
	methylation_control_average = np.nanmean(list(map(float, methylation_control_flat)))

	methylation = data['tumor_methylation'].str.split(';').tolist()
	# print(methylation)
	for i in range(len(coverage_tumor)): ## for every tumor sample
		for j in range(len(coverage_tumor[i])): ## for every CpG dinulceotide
			## TODO: Check if length of methylation, coverage and sample are the same? Currently would only throgh error if not the case
			# print(sample)
			# print(f"PRIMER:\t{primer[2]}")
			# print(f"sample: {sample[i][j]}")
			# print(f"coverage_tumor: {coverage_tumor[i][j]}")
			# print(f"methylation: {methylation[i][j]}")
			# print(f"norm: {norms[sample[i][j]]}")

			## Calculate score for each position and each tumor sample and sum up
			tmp = (float(coverage_tumor[i][j])/norms[sample[i][j]]) * float(methylation[i][j]) * 1000000 
			if not math.isnan(tmp):
				score += tmp
			# print(f"score:\t{score}")
	
	## Add the primer score
	primer.append(round(score,2))

	## Add the average coverage and methylation
	primer.append(round(coverage_control_average,2))
	primer.append(round(coverage_tumor_average,2))
	primer.append(round(methylation_control_average,2))
	primer.append(round(methylation_tumor_average,2))

	# print("FINAL PRIMER:")
	# print(primer)

	## Returned primer has the following format:
	## Chro, strand, start, end, #CGs, length, score, cov_control, cov_tumor, meth_control, meth_tumor
	return primer


def findPrimer(boxplotData, norms, minPrimerLength, maxPrimerLength, minCpGs):
	""" Function to find a primer region which fulfills criteria: min 18nt, max 24nt, min 3 CpGs """
	primerRegions = []

	print("Extract primer positions")
	## Go over all chromosomes
	for chro in boxplotData['chro'].unique(): 
		# print(f"Chromosome: {chro}")
		## Go over both strands
		for strand in ['+', '-']:

			## Extract the CpG Dinucleotide positions into a list
			positions = boxplotData.loc[(boxplotData['chro'] == chro) & (boxplotData['strand'] == strand)]['pos'].tolist()
			positions = list(map(int, positions))

			## extract the pointer and genomic positions of primer regions
			pointer_indices = getPimerPointer(positions, maxPrimerLength, minCpGs) # pointer in list for this chromosome & strand
			# print(pointer_indices)
			for primer in pointer_indices:
				# print(primer)
				numCGs = primer[1]-primer[0]+1  ## How many CGs are in the primer
				primerLength = positions[primer[1]]-positions[primer[0]]
				## Build an array of chro, start, stop, strand, #CGs, primer length
				if (primerLength >= minPrimerLength and primerLength <= maxPrimerLength):
					primerRegions.append([chro, strand] + [positions[primer[0]], positions[primer[1]], numCGs, primerLength])  
				# print(f"primer:\t{primerRegions}")

		# 	break ## only for + strand
		# break ## only for chr1

	## Add the score for each primer
	print("Calculating Primer score...")
	resultRegions = []
	primerRegions_df = pd.DataFrame(primerRegions, columns= ["chromosome", "strand", "start", "end", "num_CpGs", "length"])
	primerRegions_df.head()
	chros = primerRegions_df['chromosome'].unique()
	for chro in chros:
		# print(f"Calculating scores for chromosome {chro}")
		boxplotData_chro = boxplotData.loc[(boxplotData['chro'] == chro)]
		primerRegions_chro = primerRegions_df.loc[(primerRegions_df['chromosome'] == chro)]
		primerRegions_chro = primerRegions_chro.values.tolist()
		for primer in primerRegions_chro:
			primer = calcPrimerScore(boxplotData_chro, primer, norms)
			resultRegions.append(primer)

			## Returned primer has the following format:
			## Chro, strand, start, end, #CGs, length, score, cov_control, cov_tumor, meth_control, meth_tumor
			
	return resultRegions


def getPrimerStats(sequence):
	""" Calculate different statistics for a given DNA sequence: GC-content, non-CpG cytosines, CpGs, TM"""
	sequence = sequence.upper()
	composition = Counter(sequence)
	## Count CpG dinucleotides in given sequence
	cpgs = sequence.count("CG")
	##  Count cytosines which are not in a CpG context
	nonCpG_cytosines = composition["C"] - cpgs
	## Calculation GC-content
	gc = round((composition["C"] + composition["G"]) / (composition["A"] + composition["T"] + composition["C"] + composition["G"]),2)
	## Build a result string out of it
	res_array = [cpgs, gc, nonCpG_cytosines]
	return res_array


def findPrimerPairs(primerData, minAmpliconLength, maxAmpliconLength):
	""" Function to build all possible combinations of primer pairs to span a PCR product of length 400nt """
	pcr_products = []

	## Go linewise through primerData with pointer1 (= i)
	for i in range(len(primerData)):
		j = i
		## increment j until end of file or 400nt are reached on the same chromosome and on the same strand
		while((j < len(primerData)) and (primerData[j][0] == primerData[i][0]) and (primerData[j][1] == primerData[i][1]) and (int(primerData[j][3]) - int(primerData[i][2])) < maxAmpliconLength):
			chro = primerData[i][0]
			strand = primerData[i][1]
			length = int(primerData[j][3]) - int(primerData[i][2])
			score = round(float(primerData[i][6]) + float(primerData[j][6]),2)
			cov_contr = round(np.mean([float(primerData[i][7]), float(primerData[j][7])]),2)
			cov_tumor = round(np.mean([float(primerData[i][8]), float(primerData[j][8])]),2)
			meth_control = round(np.nanmean([float(primerData[i][9]), float(primerData[j][9])]),2)
			meth_tumor = round(np.nanmean([float(primerData[i][10]), float(primerData[j][10])]),2)
			# cov_mean = round(np.mean([float(primerData[i][7]), float(primerData[i][8]), float(primerData[j][7]), float(primerData[j][8])]),2)

			## Try to extract GC content from primers -- exists only if reference genome is given
			try:
				gc_fw = primerData[i][13]
				gc_rev = primerData[i][13]
				## Combine chro, fwPrimerData, revPrimerData, region length, score, cov_contr, cov_tumor, meth_mean, meth_tumor, gc_fw, gc_rev
				result = [chro, strand] + primerData[i][2:7] +  primerData[j][2:7] + [length, score, cov_contr, cov_tumor, meth_control, meth_tumor, gc_fw, gc_rev]
			except:
				## Combine chro, fwPrimerData, revPrimerData, region length, score, cov_contr, cov_tumor, meth_mean, meth_tumor
				result = [chro, strand] + primerData[i][2:7] +  primerData[j][2:7] + [length, score, cov_contr, cov_tumor, meth_control, meth_tumor]

			# ## Combine chro, fwPrimerData, revPrimerData, region length, score, cov_contr, cov_tumor, meth_mean, meth_tumor
			# result = [chro, strand] + primerData[i][2:7] +  primerData[j][2:7] + [length, score, cov_contr, cov_tumor, meth_control, meth_tumor]

			## Check if the resulting amplicon (=pcr product) has a minimal required length
			if( (int(primerData[j][3]) - int(primerData[i][2])) > minAmpliconLength ):
				pcr_products.append(result)
			## Check for further reverse primers
			j += 1

	## Convert nested list into pandas dataframe
	try: ## GC content information exists
		pcr_product_df = pd.DataFrame(pcr_products, columns= ["chromosome", "strand", "start_fw", "end_fw", "CGs_fw", "length_fw", "score_fw", "start_rev", "end_rev", "CGs_rev", "length_rev", "score_rev", "length_product", "score_product", "cov_contr", "cov_tumor", "meth_contr", "meth_tumor", "gc_fw", "gc_rev"])
	except: ## GC content information does not exist
		pcr_product_df = pd.DataFrame(pcr_products, columns= ["chromosome", "strand", "start_fw", "end_fw", "CGs_fw", "length_fw", "score_fw", "start_rev", "end_rev", "CGs_rev", "length_rev", "score_rev", "length_product", "score_product", "cov_contr", "cov_tumor", "meth_contr", "meth_tumor"])
	pcr_product_df_sorted = pcr_product_df.sort_values("score_product", ascending=False)
	pcr_product_df_sorted.to_csv(os.path.join(outFolder, 'pcrProducts.tsv'), index = False, sep='\t')

	return (pcr_products)


def annotate_pcrProducts(pcr_products, minScore, annotationFile):
	""" Extract the gene name for each gene overlapping the pcr products with a score above the given one """
	try: ## GC content information exists
		df = pd.DataFrame(pcr_products, columns= ["chromosome", "strand", "start_fw", "end_fw", "CGs_fw", "length_fw", "score_fw", "start_rev", "end_rev", "CGs_rev", "length_rev", "score_rev", "length_product", "score_product", "cov_contr", "cov_tumor", "meth_contr", "meth_tumor", "gc_fw", "gc_rev"])
	except: ## GC content information does not exist
		df = pd.DataFrame(pcr_products, columns= ["chromosome", "strand", "start_fw", "end_fw", "CGs_fw", "length_fw", "score_fw", "start_rev", "end_rev", "CGs_rev", "length_rev", "score_rev", "length_product", "score_product", "cov_contr", "cov_tumor", "meth_contr", "meth_tumor"])

	df_best = df.loc[df['score_product'] > minScore]

	## Read in the annotation file
	annotation = pd.read_csv(annotationFile, sep = '\t', index_col=False, names = ["chro", "ref", "feature", "start", "end", "score", "strand", "frame", "attribute"])
	annotated_pcrProducts = []

	## Iterate over all pcr products
	for index, row in df_best.iterrows():
		row = row.tolist()
		## Check which annotation features overlap with the pcr product
		feature = annotation.loc[(annotation['chro'] == row[0]) & (annotation['start'] <= int(row[2])) & (annotation['end'] >= int(row[3]))]
		gene = []
		for el in feature.attribute:
			# print(el)
			# print(el.split('"')[1])
			gene.append(el.split('"')[1])
		# gene = list(set(gene))
		genes = ' '.join(str(x) for x in gene)
		annotated_pcrProducts.append(row + [genes])
		# print('\t'.join(str(x) for x in (row+gene)))

		## Add the current information about the pcr product and the annotated gene name into a new list
		# annotated_pcrProducts.append(row + gene)

	try: ## GC content information exists
		pcr_product_df = pd.DataFrame(annotated_pcrProducts, columns= ["chromosome", "strand", "start_fw", "end_fw", "CGs_fw", "length_fw", "score_fw", "start_rev", "end_rev", "CGs_rev", "length_rev", "score_rev", "length_product", "score_product", "cov_contr", "cov_tumor", "meth_contr", "meth_tumor", "gc_fw", "gc_rev", "gene"])
	except: ## GC content information does not exist
		pcr_product_df = pd.DataFrame(annotated_pcrProducts, columns= ["chromosome", "strand", "start_fw", "end_fw", "CGs_fw", "length_fw", "score_fw", "start_rev", "end_rev", "CGs_rev", "length_rev", "score_rev", "length_product", "score_product", "cov_contr", "cov_tumor", "meth_contr", "meth_tumor", "gene"])
	pcr_product_df_sorted = pcr_product_df.sort_values("score_product", ascending=False)
	pcr_product_df_sorted.to_csv(os.path.join(outFolder, 'pcrProducts.tsv'), index = False, sep='\t')

	# ## Store the annotated best pcr products into a output file
	# with open(os.path.join(outFolder, 'annotated_pcrProducts.tsv'), 'w') as f:
	# 	for product in annotated_pcrProducts:
	# 		for value in product:
	# 			f.write(f"{value}\t")
	# 		f.write(f"\n")


def getSequence_pcrProducts(pcr_products, genome):
	""" Extract the DNA sequence for each MSP region and store it in a separate fasta file """
	with open(os.path.join(outFolder, 'mspRegions.fa'), 'w') as f:
		## Loop trhrough MSP regions
		for row in pcr_products:
			chro = row[0]
			strand = row[1]
			start_fw = row[2]
			end_rev = row[8]
			header = ">" + str(chro) + ":" + str(start_fw) + "-" + str(end_rev) + "(" + str(strand) + ") length: " + str(int(end_rev)-int(start_fw)) + "nt\n"
			# Extract the DNA sequence for the MSP region
			sequence = getFasta(genome, chro, strand, start_fw, end_rev)
			f.write(header)
			f.write(sequence)
			f.write("\n")
			f.flush()


### MAIN      
###############################################################################
## Parse the passed arguments
args = parser.parse_args()

bedmethylFile = args.bedmethylFile
controls = args.controls
tumors = args.tumors
outFolder = args.outfolder

minCpGs = int(args.minCpGs)   # minimum amount of differentially methylated cytosines in primer
maxMethControl = float(args.maxMethControl)   # max control methylation
minPrimerLength = int(args.minPrimerLength)	  # min primer length
maxPrimerLength = int(args.maxPrimerLength)   # max primer length
minAmpliconLength = int(args.minAmpliconLength)   # min MSP region length
maxAmpliconLength = int(args.maxAmpliconLength)   # max MSP region length

minScore = 1

## Check if optional arguments are given
reference_given = False
annotation_given = False
boxplotData_given = False
primerData_given = False
if args.reference is not None:
	genomeFile = args.reference
	reference_given = True
if args.annotation is not None:
	annotationFile = args.annotation
	annotation_given = True
if args.boxplotData is not None:
	boxplotDataFile = args.boxplotData
	boxplotData_given = True
if args.primerData is not None:
	primerDataFile = args.primerData
	primerData_given = True
print(f"reference_given: {reference_given}")
if(reference_given):
	print(f"genomeFile: {genomeFile}")
	# Read in the reference genome fasta
	print("Reading in the genome sequence...")
	genome = importFasta(genomeFile)
print(f"annotation_given: {annotation_given}")
if(annotation_given):
	(f"annotationFile: {annotationFile}")
print(f"boxplotData_given: {boxplotData_given}")
if(boxplotData_given):
	print(f"boxplotDataFile: {boxplotDataFile}")
print(f"primerData_given: {primerData_given}")
if(primerData_given):
	print(f"primerDataFile: {primerDataFile}")

# genomeFile = '/data/fass1/genomes/Eukaryots/homo_sapiens_done/ucsc/hg38.fa'
# annotationFile = pd.read_csv('/data/fass1/genomes/Eukaryots/homo_sapiens_done/ucsc/hg38_ncbiRefSeq.gtf', sep = '\t', index_col=False, names = ["chro", "ref", "feature", "start", "end", "score", "strand", "frame", "attribute"])

print(f"bedmethyl file:\t{bedmethylFile}")
print(f"maxMethControl:\t{maxMethControl}")
print(f"control samples:\t{controls}")
print(f"tumor samples:\t{tumors}")
print(f"outFolder:\t{outFolder}")

# beginMAIN_time = time.time()



## Primer Data already provided by the user? If so, resume from there.
##################################################################################################################################
if (primerData_given):
	primerRegions = []
	with open(primerDataFile, 'r') as f:
		lines = f.readlines()[1:]
		for line in lines:
			primerRegions.append(line.split())

else:
	# 1. Boxplot statistics
	##################################################################################################################################
	""" Calculate the boxplot statistics for each CpG Dinucleotide and keep information only for interesting CpGs """
	if (boxplotData_given):
		# Read in the boxplot statistics data
		print("Reading boxplotData ...")
		dtype_dict = {'chro': 'str', 'pos': 'int64', 'strand': 'str', 'median_ctrl' : 'float64', 'upper_quartile_ctrl': 'float64', 'lower_quartile_ctrl': 'float64', 'upper_whisker_ctrl': 'float64', 'lower_whisker_ctrl': 'float64', 'median_tmr': 'float64', 'upper_quartile_tmr': 'float64', 'lower_quartile_tmr': 'float64', 'upper_whisker_tmr': 'float64', 'lower_whisker_tmr': 'float64', 'control_coverage' : 'str', 'control_methylation' : 'str', 'control_samples' : 'str', 'tumor_coverage' : 'str', 'tumor_methylation' : 'str', 'tumor_samples' : 'str'}
		boxplotData = pd.read_csv(boxplotDataFile, sep = '\t', index_col=False, header = 0, dtype=dtype_dict)
		cpgInterestFile = boxplotDataFile.replace("boxplotData.", "interestingCpGs.")
	else:
		print("Calculating boxplot data...")
		boxplotData = getRegions(bedmethylFile, maxMethControl, controls, tumors)
		print("Finished getRegions, having full boxplot data")
		cpgInterestFile = os.path.join(outFolder, 'interestingCpGs.tsv')

		# ## Read in the boxplotData
		# print("Reading boxplotData ...")
		# dtype_dict = {'chro': 'str', 'pos': 'int64', 'strand': 'str', 'median_ctrl' : 'float64', 'upper_quartile_ctrl': 'float64', 'lower_quartile_ctrl': 'float64', 'upper_whisker_ctrl': 'float64', 'lower_whisker_ctrl': 'float64', 'median_tmr': 'float64', 'upper_quartile_tmr': 'float64', 'lower_quartile_tmr': 'float64', 'upper_whisker_tmr': 'float64', 'lower_whisker_tmr': 'float64', 'control_coverage' : 'str', 'control_methylation' : 'str', 'control_samples' : 'str', 'tumor_coverage' : 'str', 'tumor_methylation' : 'str', 'tumor_samples' : 'str'}
		# boxplotData = pd.read_csv(os.path.join(outFolder, 'boxplotData.tsv'), sep = '\t', index_col=False, header = 0, dtype=dtype_dict)

		# print(boxplotData)


	# 2. Get primer regions from boxplot data and interesting_CpGs
	##################################################################################################################################
	""" Get primer regions from boxplot data and interesting_CpGs """

	## Read in the list of interesting CpGs
	print("Reading interestingCpGs ...")
	interestingCpGs = pd.read_csv(cpgInterestFile, sep = '\t', index_col=False, names = ['chro', 'pos', 'strand', 'coverage', 'methylation', 'sample'])


	## Calculate coverage normalization values
	print("Calculate norms ...")
	norms = interestingCpGs.groupby(['sample'])['coverage'].sum()

	## Select regions which fulfill criteria to be PCR product region of interest
	print("Extract primer ...")
	primerRegions = findPrimer(boxplotData, norms, minPrimerLength, maxPrimerLength, minCpGs)

	## Add sequence stats if reference is given.
	if (reference_given):

		## Extract sequence information for the primers
		print("Caclulating sequence stats per Primer")
		annotated_primers = []
		## Iterate over all Primers
		for row in primerRegions:
			chro = row[0]
			strand = row[1]
			start = row[2]
			end = row[3]
			# print(row)

			# Extract the DNA sequence for the Primer
			sequence = getFasta(genome, chro, strand, start, end)
			stats = getPrimerStats(sequence)
			# print(row + [sequence] + stats)

			# Add the current information about the pcr product and the annotated gene name into a new list
			annotated_primers.append(row + [sequence] + stats)
			## Overwrite annotated_primers with the extended primerData
			primerRegions = annotated_primers

		## Write the primer data into a new file
		with open(os.path.join(outFolder, 'primerData.tsv'), 'w') as f:
			f.write(f"Chro\tstrand\tstart\tend\t#CGs\tlength\tscore\tcov_control\tcov_tumor\tmeth_control\tmeth_tumor\tsequence\tcpgs\tgc-content\tnonCpG_cytosines\n")
			for primer in annotated_primers:
				for value in primer:
					f.write(f"{value}\t")
				f.write(f"\n")
	
	else: ## no reference given
		## Write the primer data without sequence information into a new file
		with open(os.path.join(outFolder, 'primerData.tsv'), 'w') as f:
			f.write(f"Chro\tstrand\tstart\tend\t#CGs\tlength\tscore\tcov_control\tcov_tumor\tmeth_control\tmeth_tumor\n")
			for primer in primerRegions:
				for value in primer:
					f.write(f"{value}\t")
				f.write(f"\n")




# 3. Find all potential combinations of primers --> and store these as pcr products
##################################################################################################################################
print("Find primer pairs...")
pcr_products = findPrimerPairs(primerRegions, minAmpliconLength, maxAmpliconLength)




## Write the DNA sequences (+-200nt) from the PCR products into separate file
# print("Extracting DNA sequence for PCR products +-200nt")
# extractPCRfasta(pcr_products, genome)

## Extract DNA sequence for MSP region and store it in a new file mspRegions.fa
if (reference_given):
	getSequence_pcrProducts(pcr_products, genome)

## Add for each MSP region all overlapping annotated genes and store in new file annotated_pcrProducts.tsv
if (annotation_given):
	annotate_pcrProducts(pcr_products, minScore, annotationFile)


## TODO: make usage of primer data possible


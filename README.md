# diffONT
diffONT is a python-based tool for predicting methylation-specific PCR (MSP) primers, based on Nanopore sequencing data. Given multiple bedmethyl files, diffONT detects methylation-specific PCR primer regions, which can distinguish between two groups of samples (originally cancer patients vs. healthy controls).

Input for diffONT is a bedmethyl file, sorted by genomic position, with an additional column containing the sample name. This file can be generated using the script `preprocess.sh`, which extracts the sample name from the file name.
Output of diffONT is list of predicted MSP regions, containing information for the forward and reverse primer. This list might contain overlapping MSP regions, which can be collapsed with the script `groupPCRproducts.py`.

### Install:

### How to use:



# diffMONT

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.18914.svg)](https://doi.org/10.5281/zenodo.14501611)

<img align="right" width="250" alt="diffMONT_logo" src="https://github.com/user-attachments/assets/d77466f2-a333-4ffe-ba41-f9ac4476fdef" />

diffMONT is a python-based tool for predicting methylation-specific PCR (MSP) primers, based on Nanopore sequencing data. Given a merged bedmethyl file, diffMONT detects methylation-specific PCR primer regions, which can distinguish between two groups of samples (originally cancer patients vs. healthy controls).



### Install:
We recommend to use linux and miniconda for the enviroment management
1. [Download and install Conda.](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
2. Download the enviroment yml file [diffMONT.yml](https://github.com/meydaria/diffMONT_paper/blob/main/diffMONT.yml)
3. Create the Conda environment. 

    `conda env create -f diffMONT.yml`

4. Activate the Conda environment. You will need to activate the Conda environment in each terminal in which you want to use diffMONT.

    `conda activate diffMONT`

### Run:
Input for diffMONT is a bedmethyl file, sorted by genomic position, with an additional column containing the sample name. This file can be generated using the script `preprocess.sh`, which extracts the sample name from the file name.
The most basic usage of diffMONT is via:   

`python diffMONT.py --bedmethylFile mergedMethylation.bed --controls ctr_1 ctr_2 --tumors tmr_1 tmr_2 --outfolder results/`

The output of diffMONT is a list of predicted MSP regions, containing information for the forward and reverse primer. This list might contain overlapping MSP regions, which can be collapsed with the script `groupPCRproducts.py`.

#### Required parameters:
diffMONT has four required parameters.  
* `--bedmethylFile`  
The merged bedmethyl file should follow the bedmethyl [bed9+2 format](https://www.encodeproject.org/data-standards/wgbs/), with an additional last column specifying the sample names.
* `--controls`  
The sample names of the control group, separated by space.
* `--tumors`  
The sample names of the non-control group, separated by space.
* `outfolder`
Path to the where the results will be stored.

#### Optional parameters:
Optional parameters exist to adapt the requirements for the predicted MSP regions based on the user's requirements.
| command               | function                                                        |
|-----------------------|-----------------------------------------------------------------|
| --maxMethControl      | maximum average methylation allowed for the control sample      |
| --minCpGs             | minimum amount of differentially methylated cytosines in primer |
| --minPrimerLength     | minimum length required for primers                             |
| --maxPrimerLength     | maximum length required for primers                             |
| --minAmpliconLength   | minimum length required for MSP regions                         |
| --maxAmpliconLength   | maximum length required for MSP regions                         |

Additionally, with the optional parameters `--boxplotData` and `--primerData` it is possible to continue the work from intermediate results.
Further functionalities can be used by adding `--annotation` followed by an Ensembl annotation file, which will be screened for genes overlapping the predicted MSP regions (this will result in an additional output column). If the parameter `--reference` is set, followed by the reference genome, additional statistics like GC-content will be calculated and reported in additional output columns.

### Output file format:
The output file pcrProducts.tsv has the following output columns:
chromosome	strand	start_fw	end_fw	CGs_fw	length_fw	score_fw	start_rev	end_rev	CGs_rev	length_rev	score_rev	length_product	score_product	cov_contr	cov_tumor	meth_contr	meth_tumor

| **column name**       | **meaning**                                                     |
|-----------------------|-----------------------------------------------------------------|
| chromosome            | chromosome on which the MSP region was identified               |
| strand                | strand on which the MSP region was identified                   |
| start_fw              | first poistion of the first primer in pair                      |
| end_fw                | second poistion of the first primer in pair                     |
| CGs_fw                | number of CpGs of interest in the first primer in pair          |
| length_fw             | length of the first primer in pair                              |
| score_fw              | score of the first primer in pair                               |
| start_rev             | first poistion of the second primer in pair                     |
| end_rev               | second poistion of the second primer in pair                    |
| CGs_rev               | number of CpGs of interest in the second primer in pair         |
| length_rev            | length of the second primer in pair                             |
| score_rev             | score of the second primer in pair                              |
| length_product        | length of the predicted MSP region (end_rev-start_fw)           |
| score_product         | MSP region score (score_fw + score_rev)                         |
| cov_contr             | mean coverage of control samples over both primers              |
| cov_tumor             | mean coverage of control samples over both primers              |
| meth_contr            | mean methylation of control samples over both primers           |
| meth_tumor            | mean methylation of control samples over both primers           |

### Exemplary data for testing:
diffMONT can be either used on publicly available data, e.g. the Benchmark Dataset RRMS 2022.07 from Oxford Nanopore Technologies (ONT) available under AWS (s3://ont-open-data/rrms_2022.07/) and described by [ONT](https://epi2me.nanoporetech.com/rrms2022.07/).
Additionally, a minimal example of a bedmethyl file is available in the folder `/example_data/`. 

### Application and Publication
For a detailed description and benchmarking of diffMONT, see our [preprint on BioRxiv (2025)](https://www.biorxiv.org/content/10.1101/2025.02.17.638597v1.abstract).     
For the application of diffMONT on clinical data, see our [publication in Clinical Epigenetics (2025)](https://link.springer.com/article/10.1186/s13148-025-01960-7).

### Citation
When using diffMONT helps you, please cite:      
[Meyer *et al.*, "Nanopore sequencing-derived methylation biomarker prediction for methylation-specific PCR in patients with head and neck squamous cell carcinoma", Clinical Epigenetics (2025)](https://link.springer.com/article/10.1186/s13148-025-01960-7)


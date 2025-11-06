## This folder contains exemplary data
* `humanCpG_islands_sorted.bed`:   a sorted listed of human CpG islands
* `methylation_data.bed`:          exemplary (in silico) methylation data for 3 control and 3 tumor samples in bedmethyl file format

## This data can be used to test diffMONT.
`python diffMONT.py --bedmethylFile methylation_data.bed --controls ctr_1 ctr_2 ctr_3 --tumors tmr_1 tmr_2 tmr_3 --outfolder results/`

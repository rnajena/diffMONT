## Runnable exmaple data for diffMONT

This folder contains the following exemplary data
* `humanCpG_islands_sorted.bed`:   a sorted listed of human CpG islands
* `methylation_data.bed`:          exemplary (in silico) methylation data for 3 control and 3 tumor samples in bedmethyl file format

The data can be used to test diffMONT in the following way:
`python diffMONT.py --bedmethylFile example_data/methylation_data.bed --controls ctr_1 ctr_2 ctr_3 --tumors tmr_1 tmr_2 tmr_3 --outfolder example_data/`

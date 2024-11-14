## Check the PCR product results for overlapping resutls
## Currently if a region is good, multiple potential PCR primer regions are often returned for that region.
## With this script I want to exclude subotpimal PCR primer regions in close proximity

import argparse
from operator import itemgetter

# Instantiate the parser
parser = argparse.ArgumentParser(description='Extract only optimal PCR primer region per region')

parser.add_argument('--pcrProductfile',  
                    help='Text file containing the output of diffONT; a tab separated file containing informations about all PCR product regions')


def merge_overlappingIntervals(interval_list):
    merged = []
    ## Sort the list by the first element
    # sortedList = sorted(interval_list, key = itemgetter(0))
    sortedList = sorted(interval_list, key=lambda tup: tup[0])
    for higher in sortedList:
        if not merged:
            ## always add the first interval
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = [lower[0], upper_bound]  # replace by merged interval
            ## add new element if no overlap exists
            else:
                merged.append(higher)
    return merged

 
# testList = [[1,3], [5,9], [2,3], [5,9], [4,5], [9,10]]
# print(merge_overlappingIntervals(testList))
# testList = [[1,3], [2,3], [4,5], [5,9], [9,10]]

## Initialize variables
chro_dict = {}
chro_dict_merged = {}
resList = []

## Extract the information from input file
args = parser.parse_args()
pcrProductfile = args.pcrProductfile

## 1. Create a dict with all overlapping regions
with open(pcrProductfile, 'r') as infile:
    for line in infile:
        if line.startswith("chromosome"):
            continue
        else:
            # print(line)
            # print(line.split())
            chro = line.split()[0]
            start = int(line.split()[2])
            end = int(line.split()[8])
            # igv_koords = chro + ':' + str(start) + '-' + str(end)
            # line = line.strip() + '\t' + igv_koords
            if not (chro in chro_dict.keys()):
                chro_dict[chro] = []
            chro_dict[chro].append([start, end])


# print("Before Sorting:")
# for chro in chro_dict.keys():
# 	for entry in chro_dict[chro]:
# 		print(f"{chro}\t{entry}")

# print("Before Sorting:")
# for chro in chro_dict.keys():
#     print(f"{chro}\t{len(chro_dict[chro])}")
	# print(chro)
	# print(chro_dict[chro])

## Merge the overlapping regions for each chromosome
for chro in chro_dict.keys():
	## merge all overlapping intervalls
    chro_dict_merged[chro] = merge_overlappingIntervals(chro_dict[chro])
    ## add a "False" to each entry
    for entry in chro_dict_merged[chro]:
        entry.append(False)


# print("After Sorting:")
# for chro in chro_dict_merged.keys():
#     print(f"{chro}\t{len(chro_dict_merged[chro])}")
	# for entry in chro_dict_merged[chro]:
	# 	print(f"{chro}\t{entry}")

## 2. Check for each PCR region if a PCR region from the same overlap region has already been reported
with open(pcrProductfile, 'r') as infile:
    for line in infile:
        # print(line)
        # print(line.split())
        if line.startswith("chromosome"):
            short_line = line.split()[2:]
            short_line = "\t".join(short_line)
            print('chromosome\tstart_region\tend_region\tstrand\t' + short_line.strip() + '\tbest_koords\tregion_koords')
            continue
        else:
            chro = line.split()[0]
            start = int(line.split()[2])
            end = int(line.split()[8])
            igv_koords = chro + ':' + str(start) + '-' + str(end)
            line = line.strip() + '\t' + igv_koords

	        ## 2. Add only the best PCR region for each overlapping region
	        ## Check if region on this chromosome exists already
            if chro in chro_dict_merged.keys():
                ## check for all already listed regions if there is an overlap
                for entry in chro_dict_merged[chro]:
                    # print(f"{chro}\t{entry}")
                    ## check if any region is overlapped
                    if(((start <= entry[0]) and (entry[0] <= end)) or ((entry[0] <= start) and (start <= entry[1]))):
                        # print(entry)
                        # print(start)
                        ## regions are overlapping. If something from this overlap has been used before, do not add.
                        if entry[2] == False:
                            # print(entry)
                            ## Add the koordinates from the complete region
                            # line = line + '\t' + chro + ':' + str(entry[0]) + '-' + str(entry[1])
                            short_line = line.split()[2:]
                            short_line = "\t".join(short_line)
                            line = chro + '\t' + str(start) + '\t' + str(end) + '\t' + line.split()[1] + '\t' + short_line + '\t' + chro + ':' + str(entry[0]) + '-' + str(entry[1])
                            resList.append(line)
                            entry[2] = True
                        else:
                            continue

# print(resList)
for topRegion in resList:
	print(topRegion)


## TODO: write output to file
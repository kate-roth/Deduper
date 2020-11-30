#!/usr/bin/env python

import argparse
import re
import sys


# function to convert -p string to boolean
def str_to_bool(value):
    '''takes a string and converts to a boolean'''
    if value.lower() in {'false', 'f', '0', 'no', 'n'}:
        return False
    elif value.lower() in {'true', 't', '1', 'yes', 'y'}:
        return True
    raise ValueError(f'{value} is not a valid boolean value')

def get_args():
    parser = argparse.ArgumentParser(description="A program to identify and remove PCR duplicates from a SAM file. The SAM file must be sorted by refname.")
    parser.add_argument("-f", "--file", help="input sorted SAM file", required=True, type=str)
    parser.add_argument("-u", "--umi", help="file containing UMIs", required=True, type=str)
    parser.add_argument("-p", "--paired", help="-p True for paired-end reads", type=str_to_bool, nargs='?', const=True, default=False)
    return parser.parse_args()

# store args
args = get_args()
samfile = args.file
umifile = args.umi
paired = args.paired

# raise error and exit script if user requests paired-end
if paired == True:
    sys.exit("Error: pair-end capability not implemented yet\nTerminated")

# functions
def get_strandedness(bitwise_flag):
    '''takes the bitwise flag and returns a boolean for whether 
    the read is aligned on the forward strand'''
    
    if ((bitwise_flag & 16) != 16):
        return True
    else:
        return False

def adjust_pos(position, forward, cigar_string):
    '''Considers strandedness and returns an adjusted position if there
    is soft-clipping, insertions, deletions etc. denoted in the CIGAR string.
    Returns the original position if adjustment not necessary.'''
    
    # split and components from CIGAR string
    split_cigar = re.findall(r"((\d+)(\w))", cigar_string)
    # example: [('40M', '40', 'M'), ('1I', '1', 'I'), ('30M', '30', 'M'), ('3S', '3', 'S')]

    return_pos = position
    if forward:
        # check for soft-clipping at beginning of CIGAR string
        if "S" in split_cigar[0]:
            return_pos -= int(split_cigar[0][1])
    else:
        # check for insertions, deletions, or skips
        if "I" or "D" or "N" in split_cigar:
            for j in range(len(split_cigar)):
                if "I" in split_cigar[j]:       # subtract insertion
                    return_pos -= int(split_cigar[j][1])
                elif "D" in split_cigar[j]:     # add deletion
                    return_pos += int(split_cigar[j][1])
                elif "N" in split_cigar[j]:     # add skip
                    return_pos += int(split_cigar[j][1])
        # check for soft-clipping at the end of CIGAR string
        if "S" in split_cigar[-1]:              # add soft-clipping
            return_pos += int(split_cigar[-1][1])

    return return_pos


# store UMIs in a set
umis = set()
with open(umifile) as ufile:
    for line in ufile:
        line = line.strip("\n")
        umis.add(line)

# initalize global current refname and temporary dictionary
current_refname = ""
# keys are UMIs, values are set of tuples (position, fwd boolean) 
temp_dict = {}  # Example: {'CTGTTCAC': {(76814284, True), (76814287, True)}}

# create output filenames
split_path = re.split('/', samfile)
split_name = re.split("\.", split_path[-1])
unknown_filename = split_name[0] + "_deduped_unknown.sam"
output_filename = split_name[0] + "_deduped.sam"
duplicates_filename = split_name[0] + "_deduped_duplicates.sam"

# open files
out_unknown = open(unknown_filename, "w")
out_SAM = open(output_filename, "w")
out_dupl = open(duplicates_filename, "w")
input_SAM = open(samfile, "r")

# counters
ln = 0
unmapped = 0
unkn_count = 0
not_dupl_count = 0
dupl_count = 0

# main
for line in input_SAM:
    if line.startswith("@"):       # write header lines to output files
        out_unknown.write(line)
        out_SAM.write(line)
        out_dupl.write(line)
    else:
        ln += 1

        # store variables
        fields = re.split("\t", line)
        UMI = re.split(":", fields[0])[-1]  # end of column 1
        flag = int(fields[1])               # column 2
        fwd = get_strandedness(flag)        # extract from the flag
        refname = fields[2]                 # column 3
        pos = int(fields[3])                # column 4, will be adjusted based on CIGAR string
        cigarstr = fields[5]                # column 6

        # clear temp_dict if have moved onto a new chromosome/refname
        if refname != current_refname:
            temp_dict = {}
            current_refname = refname
        
        if UMI in umis:
            if ((flag & 4) != 4):   # if read is mapped
                pos = adjust_pos(pos, fwd, cigarstr)    # adjust position if necessary
                pos_strand = (pos, fwd)                 # create tuple for next steps
                if UMI in temp_dict:
                    if pos_strand in temp_dict[UMI]:    # write out duplicate reads
                        dupl_count += 1
                        out_dupl.write(line)
                    else:                               # write out unique reads
                        not_dupl_count += 1
                        out_SAM.write(line)
                        temp_dict[UMI].add(pos_strand)  # add tuple to temp_dict
                
                else:
                    not_dupl_count += 1                 # write out unique reads
                    out_SAM.write(line)
                    temp_dict[UMI] = set()              # add UMI:(pos, True/False) to temp_dict
                    temp_dict[UMI].add(pos_strand)
            else:                                       # count unmapped reads
                unmapped += 1
        else:                                           # write out reads with unknown UMIs
            unkn_count += 1
            out_unknown.write(line)


# close files
out_unknown.close()
out_SAM.close()
out_dupl.close()
input_SAM.close()


# print statements
print("total reads:", ln)
print("unmapped reads:", unmapped)
print("unknown UMIs:", unkn_count)
print("unique reads:", not_dupl_count)
print("duplicates:", dupl_count)
print()
print("percent unmapped reads:", str(format((unmapped/ln*100), '.1f')), "%")
print("percent unknown UMIs:", str(format((unkn_count/ln*100), '.1f')), "%")
print("percent unique reads:", str(format((not_dupl_count/ln*100), '.1f')), "%")
print("percent PCR duplicates:", str(format((dupl_count/ln*100), '.1f')), "%")

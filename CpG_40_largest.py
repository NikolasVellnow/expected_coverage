import time
from argparse import ArgumentParser
from Bio import SeqIO

import numpy as np
import pandas as pd

MAX_COV_VALUE = 10000



def parse_args():
    """ Parses gzipped input file from which lines are then later read"""
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True)
    parser.add_argument("-o", "--output_file", required=False, default="default_output.txt")
    args = parser.parse_args()
    return args

def get_header_from_file(file_name: str) -> list[str]:
    """ Gets the sample names from first line of gzip file"""

    with open(file_name, mode='rt') as file:
        header = file.readline()
        
    return header


def process_gzip_for(file_name: str, bin_labels, sample_list):
    """ Processes gzip file line by line using a for loop """

    bin_value_list = [np.zeros(MAX_COV_VALUE+2) for sample in sample_list]
    with gzip.open(file_name, mode='rt') as file:
        for line in file:
            splitted_line = line.strip().split('\t')
            if splitted_line[0] != "CHROM":
                cov_values = splitted_line[2::]
                cov_values = list(map(int, cov_values))
                # loops over coverage values for each sample with its corresponding np.array of so-far-accumulated bin count  


def main():
    """ Main function"""
    args = parse_args()
    data_file = args.input_file
    out_file = args.output_file

    largest_40 = []
   
    f1 = SeqIO.parse(data_file, 'fasta')
    for seq_record in f1:
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))
        largest_40.append(seq_record)
        largest_40.sort(key=lamda record: len(record))
    
 
    #with open(data_file, mode='rt') as file:
        #seq_data_list = []
        # for line in file:
            #line = line.strip()
            #largest_40.append((line, len(line)))
            #largest_40.sort(key=lambda seq_data: seq_data[1])
            #largest_40 = largest_40[0:39]
            


    print(largest_40)


if __name__ == "__main__":
    main()

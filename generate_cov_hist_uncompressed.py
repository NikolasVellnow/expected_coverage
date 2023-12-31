import time
from argparse import ArgumentParser

import numpy as np
import pandas as pd


def add_bin_count(cov_value: int, bin_labels: list[int], bin_counts: np.ndarray):
    """ Adds one count to a histogram bin when that coverage vaule was observed """
    # if coverage value is <=1000, add count to corresponding bin
    try:
        idx = bin_labels.index(cov_value)
        bin_counts[idx] += 1
    # if coverage value is >1000 (or just not in bin_counts), add count to "1001-bin"
    except ValueError:
        bin_counts[1001] += 1
    return bin_counts

def parse_args():
    """ Parses gzipped input file from which lines are then later read"""
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True)
    parser.add_argument("-o", "--output_file", required=False, default="default_output.txt")
    args = parser.parse_args()
    return args

def get_samples_from_file(file_name: str) -> list[str]:
    """ Gets the sample names from first line of file"""

    with open(file_name, mode='rt', encoding="utf-8") as file:
        line = file.readline()
        splitted_line = line.strip().split('\t')
        sample_list = splitted_line[2::]
        
    return sample_list


def process_file(file_name: str, bin_labels, sample_list):
    """ Processes file line by line using a for loop """

    bin_value_list = [np.zeros(1002) for sample in sample_list]
    with open(file_name, mode='rt', encoding="utf-8") as file:
        for line in file:
            splitted_line = line.strip().split('\t')
            if splitted_line[0] != "CHROM":
                cov_values = splitted_line[2::]
                cov_values = list(map(int, cov_values))
                for sample_value, sample_bin_values in zip(cov_values, bin_value_list):
                    add_bin_count(cov_value=sample_value,bin_labels=bin_labels, bin_counts=sample_bin_values)
            
    return bin_value_list

def main():
    """ Main function"""

    t0 = time.perf_counter()
    args = parse_args()

    data_file = args.input_file
    out_file = args.output_file

    t1 = time.perf_counter()
    print("Parsing arguments: ", t1-t0, "s")

    # create list of bins. Bin 1001 is for any value >1000
    bins = list(range(0,1002))

    t2 = time.perf_counter()
    print("Preprocessing: ", t2-t1, "s")

    
    samples = get_samples_from_file(data_file)

    bin_value_list = process_file(data_file, bins, samples)
    #bin_value_list = process_gzip_while(data_file, bins, samples)

    t3 = time.perf_counter()
    print("Data file processing: ", t3-t2, "s")

    hist_df = pd.DataFrame(bin_value_list)
    hist_df = hist_df.transpose()
    hist_df.columns = samples
    hist_df['cov_bin'] = hist_df.index

    hist_df.to_csv(out_file, sep="\t", index=False)



if __name__ == "__main__":
    main()

import gzip
import time
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np


def add_bin_count(cov_value: str, bin_labels: list[str], bin_counts: np.ndarray):
    """ Adds one count to a histogram bin when that caoverage vaule was observed """
    idx = bin_labels.index(cov_value)
    bin_counts[idx] += 1
    return bin_counts

def parse_args():
    """ Parses gzipped input file from which lines are then later read"""
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True)
    parser.add_argument("-o", "--output_file", required=False)
    args = parser.parse_args()
    return args




def main():
    """ Main function"""

    t0 = time.perf_counter()
    args = parse_args()

    

    data_file = args.input_file

    t1 = time.perf_counter()
    print("Parsing arguments: ", t1-t0, "s")

    bins = [str(x) for x in range(0,1002)]
    #bins.append('>1000')
    # create array of zeros to which
    bin_values = np.zeros(1002)
    samples=['s10', 's11', 's12', 's13', 's14', 's15', 's16', 's17', 's18', 's19', 'S1', 's20', 'S2', 's3', 's4', 's5', 's6', 's7', 's8', 's9']

    bin_value_list = [np.zeros(1002) for sample in samples]

    t2 = time.perf_counter()
    print("Preprocessing: ", t2-t1, "s")
    
    with gzip.open(data_file, mode='rt') as file:
        for line in file:
            splitted_line = line.strip().split('\t')
            cov_values = splitted_line[2::]
            if splitted_line[0] == "CHROM":
                samples = cov_values
            #print(cov_values)
            else:
                
                for sample_value, sample_bin_values in zip(cov_values, bin_value_list):
                    bin_values = np.zeros(1002)
                    add_bin_count(cov_value=sample_value,bin_labels=bins, bin_counts=sample_bin_values)


    t3 = time.perf_counter()
    print("Data file processing: ", t3-t2, "s")

    # bin_ints = [int(x) for x in bins]
    # plt.hist(bins, bins=bin_ints, weights = bin_value_list[5])
    # plt.xlim(xmin=0.0, xmax=150)
    # plt.xlabel("Coverage for individual positions")
    # plt.ylabel("Counts")
    # plt.show()




if __name__ == "__main__":
    main()

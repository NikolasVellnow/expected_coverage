import gzip
import time
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def add_bin_count(cov_value: str, bin_labels: list[str], bin_counts: np.ndarray):
    """ Adds one count to a histogram bin when that caoverage vaule was observed """
    try:
        idx = bin_labels.index(cov_value)
        bin_counts[idx] += 1
    except ValueError:
        print("Coverage value could not be found in bin labels. Maybe the coverage value was higher thatn 1000?")
        bin_counts[1001] += 1
    return bin_counts

def parse_args():
    """ Parses gzipped input file from which lines are then later read"""
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True)
    parser.add_argument("-o", "--output_file", required=False, default="default_output.txt")
    args = parser.parse_args()
    return args


def main():
    """ Main function"""

    t0 = time.perf_counter()
    args = parse_args()

    data_file = args.input_file
    out_file = args.output_file

    t1 = time.perf_counter()
    print("Parsing arguments: ", t1-t0, "s")

    bins = [str(x) for x in range(0,1002)]

    t2 = time.perf_counter()
    print("Preprocessing: ", t2-t1, "s")
    
    with gzip.open(data_file, mode='rt') as file:
        for line in file:
            splitted_line = line.strip().split('\t')
            cov_values = splitted_line[2::]
            if splitted_line[0] == "CHROM":
                samples = cov_values
                bin_value_list = [np.zeros(1002) for sample in samples]
            else:
                for sample_value, sample_bin_values in zip(cov_values, bin_value_list):
                    add_bin_count(cov_value=sample_value,bin_labels=bins, bin_counts=sample_bin_values)


    t3 = time.perf_counter()
    print("Data file processing: ", t3-t2, "s")

    hist_df = pd.DataFrame(bin_value_list)
    hist_df = hist_df.transpose()
    hist_df.columns = samples
    hist_df['cov_bin'] = hist_df.index

    hist_df.to_csv(out_file, sep="\t")

    # bin_ints = [int(x) for x in bins]
    # plt.hist(bins, bins=bin_ints, weights = bin_value_list[5])
    # plt.xlim(xmin=0.0, xmax=150)
    # plt.xlabel("Coverage for individual positions")
    # plt.ylabel("Counts")
    # plt.show()


if __name__ == "__main__":
    main()

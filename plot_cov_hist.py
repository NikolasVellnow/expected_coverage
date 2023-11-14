import time
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    """ Parses arguments """
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True)
    parser.add_argument("-o", "--output_file", required=False, default="cov_hists.txt")
    args = parser.parse_args()
    return args




def main():
    args = parse_args()
    input_file = args.input_file
    df = pd.read_csv(input_file, sep="\t")
    print(df.head())



if __name__ == "__main__":
    main()
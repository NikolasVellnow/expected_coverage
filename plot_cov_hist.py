from argparse import ArgumentParser

import matplotlib.pyplot as plt
import pandas as pd


def parse_args():
    """ Parses arguments """
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_file", required=True)
    parser.add_argument("-o", "--output_file", required=False, default="cov_hists.txt")
    args = parser.parse_args()
    return args

def plot_subplots(df, sample_list, n_rows, n_cols):

    # start multiplot figure
    figure, axis = plt.subplots(n_rows, n_cols, sharex=True, sharey=True)

    sample_idx = 0
    for row in range(0, n_rows):
        for col in range(0, n_cols):
            try:
                axis[row, col].hist(df['cov_bin'], weights = df[sample_list[sample_idx]], bins=df['cov_bin'])
                axis[row, col].set_title(sample_list[sample_idx], fontsize=8)
            except IndexError:
                pass
            sample_idx+=1
    for ax in axis.flat:
        ax.label_outer()
  
    # Set common labels
    figure.text(0.5, 0.04, 'Coverage for individual positions', ha='center', va='center')
    figure.text(0.06, 0.5, 'Counts', ha='center', va='center', rotation='vertical')

    plt.xlim(xmin=0.0, xmax=200)
    plt.tight_layout(pad=3.3, w_pad=0.1, h_pad=0.2)
    plt.show()


def main():
    """ Main function """
    args = parse_args()
    input_file = args.input_file
    df = pd.read_csv(input_file, sep="\t", dtype='int32')
    print(df.head())

    samples = sorted(list(df.columns)[:-1:])
    print(samples)

    plot_subplots(df, samples, 4, 5)


if __name__ == "__main__":
    main()
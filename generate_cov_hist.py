import gzip
from argparse import ArgumentParser

PARSER = ArgumentParser()
PARSER.add_argument("-i", "--input_file", required=True)
PARSER.add_argument("-o", "--output_file", required=False)
ARGS = PARSER.parse_args()

data_file = ARGS.input_file


with gzip.open(data_file, mode='rt') as file:
    for line in file:
        print(line)




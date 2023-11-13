import time
from argparse import ArgumentParser

import modin as md
#import pandas as pd
import modin.pandas as pd

PARSER = ArgumentParser()
PARSER.add_argument("-i", "--input_file", required=True)
PARSER.add_argument("-o", "--output_file", required=False)
ARGS = PARSER.parse_args()

data_file = ARGS.input_file
output_file = ARGS.output_file

t0 = time.time()

dtype_dict={'CHROM': 'category',
            'POS': 'int32',
            'S1': 'float32',
            'S2': 'float32',
            's3': 'float32',
            's4': 'float32',
            's5': 'float32',
            's6': 'float32',
            's7': 'float32',
            's8': 'float32',
            's9': 'float32',
            's10': 'float32',
            's11': 'float32',
            's12': 'float32',
            's13': 'float32',
            's14': 'float32',
            's15': 'float32',
            's16': 'float32',
            's17': 'float32',
            's18': 'float32',
            's19': 'float32',
            's20': 'float32'
            }

df = pd.read_csv(data_file, sep='\t', header=0, dtype=dtype_dict)


t1 = time.time()
print('Read in data: ',t1-t0, 's')

df = df.drop(columns='Unnamed: 22')

#df['CHROM'] = df['CHROM'].astype('category')
#df['POS'] = pd.to_numeric(df['POS'], downcast="unsigned")
#df['POS'] = df['POS'].astype('int32')
print(df.dtypes)
sample_list = list(df.columns.values)[2:]
#df[sample_list] = df[sample_list].astype('float32')



t2= time.time()

print("Process df: ", t2-t1, 's')

frame={'mean_cov': df[sample_list].mean(axis=0),
       'sd_cov': df[sample_list].std(axis=0)}
result_df=pd.DataFrame(frame)

t3=time.time()

print('Creating frame: ', t3-t2, 's')

result_df['cv_cov']=result_df['sd_cov']/result_df['mean_cov']

t4=time.time()

print('Calculating CV: ',t4-t3, 's')

result_df.to_csv(output_file, index_label='sample')

t5 = time.time()

print('Writing to csv file: ', t5-t4, 's')

print("Total elapsed time: ", t5-t0, "s")


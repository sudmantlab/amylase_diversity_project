import pandas as pd
import numpy as np
import argparse

def process_superpopulation(superpopulation):
    # Read the gzipped file into a pandas DataFrame
    df = pd.read_csv("genomewide.salti.lassip.hap.stats.gz", compression='gzip', delimiter='\t')
    
    # Filter the DataFrame for the specified superpopulation
    df_filtered = df[df['Superpopulation'] == superpopulation]
    
    # Rank the 'h12' and 'h2h1' columns
    df_filtered['h12_rank'] = df_filtered['h12'].rank(ascending=False)
    df_filtered['h2h1_rank'] = df_filtered['h2h1'].rank(ascending=False)
    
    # Calculate the number of windows (rows)
    nwindows = df_filtered.shape[0]
    
    # Calculate the r_h12 and r_h2h1 columns
    df_filtered['r_h12'] = np.log(df_filtered['h12_rank'] / nwindows)
    df_filtered['r_h2h1'] = np.log(df_filtered['h2h1_rank'] / nwindows)
    
    # Calculate the Fscore column
    df_filtered['Fscore'] = -2 * (df_filtered['r_h12'] + df_filtered['r_h2h1'])
    
    # Save the filtered DataFrame to a gzipped TSV file
    output_file_path = f'Fisher_{superpopulation}_windows.tsv.gz'
    df_filtered.to_csv(output_file_path, sep='\t', index=False, compression='gzip')
    
    return df_filtered.head()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process superpopulation data.')
    parser.add_argument('superpopulation', type=str, help='The superpopulation to process (e.g., WEA, AFR, CAS, EA).')
    
    args = parser.parse_args()
    superpopulation = args.superpopulation
    
    result = process_superpopulation(superpopulation)
    print(result)


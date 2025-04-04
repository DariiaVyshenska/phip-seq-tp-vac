#!/usr/bin/env python
import os
import pandas as pd
import argparse
import re

def main(output_path, meta_path, kal_out_paths):
  meta_path = os.path.abspath(meta_path)
  output_path = os.path.abspath(output_path)

  try:
    result_df = pd.read_table(meta_path, sep=',', usecols=['target_id'])
  except Exception as e:
    print(f"Error reading metadata file: {e}")
    return

  df_list = []
  for kal_out_path in kal_out_paths:
    kal_out_path = os.path.abspath(kal_out_path)
    for root, dirs, files in os.walk(kal_out_path):
      if "abundance.tsv" in files:
        file_path = os.path.join(root, "abundance.tsv")
        try:
          curr_df = pd.read_table(file_path, sep='\t', usecols=['target_id', 'est_counts'])
        except Exception as e:
          print(f"Error reading file {file_path}: {e}")
          continue
        sample_id = os.path.basename(root)
        curr_df.rename(columns={'est_counts': sample_id}, inplace=True)
        df_list.append(curr_df)

  if not df_list:
    print("No valid data files were found.")
    return

  for df in df_list:
    try:
      result_df = result_df.merge(df, on='target_id', how='left')
    except Exception as e:
      print(f"Error merging dataframe with the following column names: {df.columns}")
      continue

  os.makedirs(output_path, exist_ok=True)

  result_df.to_csv(os.path.join(output_path, 'kallisto_raw_counts_merged.csv'), index=False)
  print("All done!")
    
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Merge raw counts from kallisto output in a single csv file.',
                                   usage='kallisto_inut_parser /path/to/kal/out /path/to/meta.csv /path/to/result')
  parser.add_argument('output_path', type=str, help='Path where the merged output CSV will be saved')
  parser.add_argument('meta_path', type=str, help='Path to the metadata file (CSV format). Must contain "target_id" column that references same peptide id as were provided to kallisto')
  parser.add_argument('kal_out_paths', type=str, nargs='+', help='Paths to the kallisto output directories')
  
  
  args = parser.parse_args()
  main( args.output_path, args.meta_path, args.kal_out_paths)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 19:00:38 2023

@author: skhanal
"""

import os
import pandas as pd
import subprocess

result_dir = "/Users/skhanal/Globus_local/MOTEP-OSW/outputs/paper_draft"
# data_dir = "/Users/skhanal/Globus_local/MOTEP-OSW/inputs/ISONE8busTS_osw1_17_N5X5"
script_dir = os.path.dirname(os.path.abspath(__file__))

# A list to hold all the DataFrames
all_dfs = []

# List all subdirectories in the result_dir
subdirs = [os.path.join(result_dir, d) for d in os.listdir(result_dir) if os.path.isdir(os.path.join(result_dir, d))]

# Process each subdirectory
for subdir in subdirs:
    print(f"Processing {subdir}...")
    try:
        # Call the process_generations.py script using subprocess
        subprocess.check_call(["python", os.path.join(script_dir, "process_generations.py"), subdir])

        # After the above script runs successfully, read the generated CSV and add to the list
        csv_path = os.path.join(subdir, "Pg_reduced.csv") # Change this to the `Pg_aq.csv` if you want to compute Pg_aq across all specs
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            all_dfs.append(df)
        else:
            print(f"CSV file was not found in {subdir}")

    except subprocess.CalledProcessError as e:
        print(f"Failed to run process_generations.py on {subdir}. Error: {e}")
        continue

# Concatenate all the DataFrames if we have any
if all_dfs:
    final_df = pd.concat(all_dfs, ignore_index=True)
    # Save the final DataFrame to one CSV file
    final_df.to_csv(os.path.join(result_dir, "generations.csv"), index=False) # # Change this to the `Pg_aq_all.csv` or sth, if you used above `Pg_aq.csv` instead of `Pg_reduced.csv`
else:
    print("No data to concatenate.")

print("Done processing all directories.")

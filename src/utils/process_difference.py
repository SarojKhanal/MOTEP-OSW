#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 23:27:09 2023

@author: skhanal
"""

import os
import glob
import pandas as pd
import shutil

# Define the directories
dir_base = "/Users/skhanal/Globus_local/MOTEP-OSW/outputs/ISONE/MO"
dir_compare = "/Users/skhanal/Globus_local/MOTEP-OSW/outputs/ISONE/MO OPOI"
dir_diff = f"{os.path.basename(dir_compare)}-{os.path.basename(dir_base)}"

# Create the Differences directory if it doesn't exist
if not os.path.exists(dir_diff):
    os.makedirs(dir_diff)

# List of filenames to compare
filenames = [
    "power_output_existing_generators.csv",
    "power_output_new_ng_cc_ccs.csv",
    "power_output_new_ng_ct.csv",
    "pv_node_out.csv",
    "wind_node_out.csv",
    "pv_node_new_out.csv",
    "wind_node_new_out.csv",
    "load_curtailments.csv",
    "wind_spillage.csv",
    "ch_node.csv",
    "dis_node.csv",
    "bus_out_power.csv",
    "load_node_out.csv",
    "flexible_load.csv"
]

# Iterate over each file
for filename in filenames:
    # Construct the file paths
    path_base = os.path.join(dir_base, filename)
    path_compare = os.path.join(dir_compare, filename)
    
    # Read the files
    df_base = pd.read_csv(path_base, header=None)
    df_compare = pd.read_csv(path_compare, header=None)
    
    # Ensure the structure is identical by aligning columns
    if set(df_base.columns) != set(df_compare.columns):
        raise ValueError(f"Columns do not match in {filename}")
    
    # Calculate the difference (base - compare)
    df_diff = df_compare - df_base
    
    # Save the difference DataFrame to the Differences directory
    df_diff.to_csv(os.path.join(dir_diff, filename), index=False, header=False)

shutil.copy(glob.glob(os.path.join(dir_base, 'cases_paramsopts_*.csv'))[0], dir_diff)
shutil.copy(os.path.join(dir_base, "line_upgrade_decision.csv"), dir_diff)

print("Differences calculated and saved.")

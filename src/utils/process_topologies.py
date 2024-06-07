#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 13:49:24 2023

@author: skhanal
"""

import os
import subprocess

### !TODO: Limitations: Mannually need to run FPOI and OPOI cases.
result_dir = "/Users/skhanal/Globus_local/MOTEP-OSW/outputs/paper_draft"
# CHANGE DATA DIR FOR MOPOI SEPARATELY
data_dir = "/Users/skhanal/Globus_local/MOTEP-OSW/inputs/ISONE8busTS_osw1_17_N5X5"

# Get the directory containing the current script
script_dir = os.path.dirname(os.path.realpath(__file__))

# List all subdirectories in the result_dir
subdirs = [os.path.join(result_dir, d) for d in os.listdir(result_dir) if os.path.isdir(os.path.join(result_dir, d))]

# Process each subdirectory
for subdir in subdirs:
    print(f"Processing {subdir}...")
    try:
        # Run process_transmission.py
        subprocess.call(["python", os.path.join(script_dir, "process_transmissions.py"), subdir, data_dir])

    except subprocess.CalledProcessError as e:
        print(f"Failed to run process_transmissions.py on {subdir}. Error: {e}")
        continue

    try:
        # Run plot.py
        subprocess.call(["python", os.path.join(script_dir, "plots.py"), subdir])
    except subprocess.CalledProcessError as e:
        print(f"Failed to run plots.py on {subdir}. Error: {e}")

print("Done processing all directories.")

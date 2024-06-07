#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 19:58:39 2023

@author: skhanal
"""
import pandas as pd
import os
from pathlib import Path
# import math
import argparse

# Parameters
n_lines_existing = 12
n_buses_existing = 8
# Mannual ovveride to insconsistent zone names
zone_fix_dict = {"SEMASS" : "SEMA" , "WCMASS" : "WCMA", "NEMA/BOST": "NEMA", "CTREV":"REV", "RIREV":"REV"}

# Input and result directories
# result_dir = '/Users/skhanal/Globus_local/MOTEP-OSW/outputs/paper_draft/MO'
# data_dir = "/Users/skhanal/Globus_local/MOTEP-OSW/inputs/ISONE8busTS_osw1_17_N5X3"
parser = argparse.ArgumentParser()
parser.add_argument('result_dir', type=str, help='The result directory')
parser.add_argument('data_dir', type=str, help='The input data directory')
args = parser.parse_args()
result_dir = args.result_dir
data_dir = args.data_dir

print("result_dir: ", result_dir)
print("data_dir: ", data_dir)

epoch = len(pd.read_csv(Path(result_dir,"line_upgrade_decision.csv"), header=None).columns) # 4
t_num = int(len(pd.read_csv(Path(result_dir, "power_output_existing_generators_node.csv"), header=None).columns)/epoch) # 24 hr
n_buses = len(pd.read_csv(Path(data_dir, "nodes.csv"))) # Baseline: 8, Base SpillCost0 = 15, Optimized POI = 14
n_lines = len(pd.read_csv(Path(data_dir, "lines.csv")))
n_lines_new = n_lines - n_lines_existing 

zone_node_dict = pd.read_csv(Path(data_dir, "nodes.csv")).set_index("Zone")["node"].to_dict()
node_zone_dict = pd.read_csv(Path(data_dir, "nodes.csv")).set_index("node")["Zone"].to_dict()

df_lines = pd.read_csv(Path(data_dir, "lines.csv"))

# Existing lines
df_lines_existing = pd.read_csv(Path(data_dir, "lines.csv")).head(n_lines_existing)
df_lines_existing['Line Type'] = "Existing Overhead HVAC"
df_lines_existing['epoch'] = ""
df_lines_existing['type'] = ""


# Existing line upgrades
df_line_upgrade_decision = pd.read_csv(Path(result_dir, "line_upgrade_decision.csv"), header=None, dtype="Int64").head(n_lines_existing)
df = df_line_upgrade_decision.copy()
# find the indices of all elements that are equal to 1
indices = df.where(df == 1).stack().reset_index()
indices['type'] = "0"
indices['epoch'] = indices['level_1'] + 1
# indices['s_max'] = 1200
df_line_upgrade_decision_indices = indices.copy()

# New (offshore) lines
try:
    df_new_line_decision = pd.read_csv(Path(result_dir, "new_line_decision.csv"), header=None, dtype="Int64").tail(n_lines_new)
    c_num = int(len(df_new_line_decision.columns) / epoch) # 3 # Line type
    df = df_new_line_decision.copy()
    # Find the indices of all elements that are equal to 1
    indices = df.where(df == 1).stack().reset_index()
    # Create a new column to record the index and column name where the value is 1
    indices['type'] = (indices['level_1'] % c_num + 1).astype(str)
    indices['epoch'] = (indices['level_1']//c_num + 1).astype(str)
    # indices['s_max'] = indices['type'].map({1:400, 2:1400, 3:2200 })
    df_new_line_decision_indices = indices.copy()
except:
    df_new_line_decision_indices = pd.DataFrame(columns=df_line_upgrade_decision_indices.columns)

# Merge the new column with the original dataframe
indices = pd.concat([df_line_upgrade_decision_indices, df_new_line_decision_indices])
indices["index"] = indices["level_0"] + 1

df_lines_sliced = pd.merge(df_lines, indices[["index", "epoch", "type"]], how="left", on="index")
df_lines_sliced = df_lines_sliced[~df_lines_sliced["epoch"].isna()]
df_lines_sliced["epoch"] = df_lines_sliced["epoch"].astype(int)
df_lines_sliced["s_max"] = df_lines_sliced["type"].map({"0": 1200, "1":400, "2": 1400, "3": 2200}) 

# Proper name may be "Line Type" instead of "Line Type"
df_lines_sliced["Line Type"] = df_lines_sliced['type'].map({"0":"Overhead HVAC", "1":"Submarine HVAC", "2": "Submarine HVDC", "3":"Submarine HVDC"})

# Create new dataframe with first row repeated 5 times
df_lines_existing_placeholder = pd.concat([df_lines_existing.iloc[0]] * 5, axis=1).transpose().reset_index(drop=True)
df_lines_existing_placeholder["Line Type"] = pd.Series(["Submarine HVDC", "Submarine HVDC", "Submarine HVAC","Overhead HVAC", "Existing Overhead HVAC"])
df_lines_existing_placeholder["From Zone"] = df_lines_existing_placeholder["To Zone"]
df_lines_existing_placeholder["s_max"]  = pd.Series([2200, 1400, 400 , 1200, 1200])

# Add existing lines for existing network topology visualization
df_lines_sliced = pd.concat([df_lines_existing_placeholder, df_lines_existing, df_lines_sliced])


df_lines_sliced['Line Type'] = df_lines_sliced['Line Type'] + ' ' + df_lines_sliced['s_max'].astype(str) + " MW"

# Fix zones to have a single node for RIREV and CTREV
df_lines_sliced['From Zone'] = df_lines_sliced['From Zone'].astype(str).replace(zone_fix_dict)
df_lines_sliced['To Zone'] = df_lines_sliced['To Zone'].astype(str).replace(zone_fix_dict)


df_lines_sliced.to_csv(Path(result_dir, "line_existing_upgrades_new.csv"), index = False, date_format="%m/%d/%Y")

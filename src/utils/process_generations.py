#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 00:34:38 2023

@author: skhanal
"""

import pandas as pd
import numpy as np
import os
# from pathlib import Path
# import math
import glob
import argparse


## JUST FOR AIR QUALITY DIFF

# # Diff


### CHANGE THESE BEFORE RUNS ###
GSw_PJM = False
# result_dir = "/Users/skhanal/Globus_local/MOTEP-OSW/outputs/paper_draft/runs/SO"

# base_case_dirname = "SO" # Base on which differentials should be computed

parser = argparse.ArgumentParser()
parser.add_argument('result_dir', type=str, help='The result directory')
args = parser.parse_args()
result_dir = args.result_dir

print("result_dir: ", result_dir)

#
params_dict = pd.read_csv(glob.glob(f"{result_dir}/cases_paramsopts_*.csv")[0]).set_index('Keys')['Values'].to_dict()
data_dir = params_dict['datadir']

# Hour => Epoch : CEILING(([Hour]+1) / 24)
n_lines_existing = int(params_dict["n_lines_existing"])
n_buses_existing = int(params_dict["n_buses_existing"])

epoch = len(pd.read_csv(f"{result_dir}/line_upgrade_decision.csv", header=None).columns) # 4
epoch_length = int(int(params_dict["Y"])/epoch)
t_num = int(params_dict["t_num"])
n_buses = len(pd.read_csv(f"{data_dir}/nodes.csv")) # Baseline: 8, Fixed POI = 15, Optimized POI = 14
n_lines = len(pd.read_csv(f"{data_dir}/lines.csv"))
e_num = int(params_dict["e_num"]) + (int(params_dict["e_numx"]) if bool(int(params_dict["GSw_ExtremeScenarios"])) == True else 0)

# Check if e_num is valid

# if int(len(pd.read_csv(f"{result_dir}/load_curtailments.csv", header=None))/n_buses) != e_num:
#     raise ValueError(f"Issue with e_num input and output files in result_dir = {result_dir}")

n_lines_new = n_lines - n_lines_existing

existing_generator_data = f"{data_dir}/generators_loc_emission_damages.csv"

zone_node_dict = pd.read_csv(f"{data_dir}/nodes.csv").set_index("Zone")["node"].to_dict()
node_zone_dict = pd.read_csv(f"{data_dir}/nodes.csv").set_index("node")["Zone"].to_dict()

# Mannual ovveride to insconsistent zone names
zone_fix_dict = {"SEMASS" : "SEMA" , "WCMASS" : "WCMA", "NEMA/BOST": "NEMA", "CTREV":"REV", "RIREV":"REV"}



# Elicit weights of representative days count for annual generation contribution
alpha_E = float(params_dict["alpha_E"])  # Replace with your actual alpha_E value
tau = 365 # float(params_dict["tau"])  # Replace with your actual tau value, if it's variable

# Read the CSV file into a pandas DataFrame
weights_of_scenarios_df = pd.read_csv(os.path.join(data_dir, "weights_of_scenarios.csv"), header=None)

# Convert the DataFrame to a numpy matrix
omega_e = weights_of_scenarios_df.values

# Split the matrix into two halves
half_length = len(omega_e) // 2
omega_N = omega_e[:half_length]
omega_E = omega_e[half_length:]

# Conditional logic as per the Julia code
if bool(int(params_dict["GSw_ExtremeScenarios"])):
    if alpha_E == -1.0:
        omega_e = omega_e
    else:
        omega_N = omega_N + omega_E
        omega_E = omega_E / omega_E.sum()
        omega_e = np.concatenate(((1 - alpha_E) * omega_N, alpha_E * omega_E))
else:
    omega_e = omega_N + omega_E

# Apply tau to omega_e
tau_e = tau * omega_e
tau_e = tau_e.flatten()


# Existing thermal power plant information
df_thermal_plants_existing = pd.read_csv(existing_generator_data)
n_generators = len(df_thermal_plants_existing)

# Dataset specific, might need change for other input data set
type_tech_dict = {'Coal' : 'Coal',
                  'Gas_CCGT' : 'NG-CC',
                  'Gas_GT' : 'NG-CT',
                  'Gas_IC' : 'NG-CT',
                  'Gas_Other' : 'NG-CC',
                  'Gas_Steam' : 'NG-CC',
                  'Nuclear' : 'Nuclear',
                  'Oil' : 'Oil'}

type_fuel_dict = {'Coal' : 'Coal',
                  'Gas_CCGT' : 'Gas',
                  'Gas_GT' : 'Gas',
                  'Gas_IC' : 'Gas',
                  'Gas_Other' : 'Gas',
                  'Gas_Steam' : 'Gas',
                  'Nuclear' : 'Nuclear',
                  'Oil' : 'Oil'}

if GSw_PJM:
    type_tech_dict = {"Combustion Turbine": "NG-CT",
                    "O/G Steam": "NG-CC",
                    "Nuclear": "Nuclear",
                    "Coal Steam": "Coal",
                    "Combined Cycle": "NG-CC"}

    type_fuel_dict = {"Combustion Turbine": "Gas",
                    "O/G Steam": "Gas",
                    "Nuclear": "Nuclear",
                    "Coal Steam": "Coal",
                    "Combined Cycle": "Gas"}

df_thermal_plants_existing["status"] = "Existing"
df_thermal_plants_existing["zone"] = df_thermal_plants_existing["Zone"].astype(str).replace(zone_fix_dict)
df_thermal_plants_existing["node"] = df_thermal_plants_existing["zone"].map(zone_node_dict)
df_thermal_plants_existing["tech"] = df_thermal_plants_existing["type"].map(type_tech_dict)
df_thermal_plants_existing["fuel"] = df_thermal_plants_existing["type"].map(type_fuel_dict) 
df_thermal_plants_existing["latitude"] = df_thermal_plants_existing["Latitude"]
df_thermal_plants_existing["longitude"] = df_thermal_plants_existing["Longitude"]
df_thermal_plants_existing["index_existing"] = df_thermal_plants_existing["index"]

zone_latitude_dict = {"CTREV":41.149972,"RIREV":41.149972,"REV":41.149972,"VINE":41.03325,"PKCTY":40.859972,"COMW":40.77634,"MFLR1":40.797153,"MFLR2":40.806968,"NH":43.685,"NEMA":42.498045,"SEMA":41.850783,"CT":41.6184324,"RI":41.671667,"WCMA":42.395356,"VT":43.926667,"ME":45.253333}
zone_longitude_dict = {"CTREV":-71.069972, "RIREV":-71.069972, "REV":-71.069972, "VINE":-70.61667, "PKCTY":-70.709972, "COMW":-70.695545, "MFLR1":-70.616667, "MFLR2":-70.547001, "NH":-71.577222, "NEMA":-71.361273, "SEMA":-70.780095, "CT":-72.7137063, "RI":-71.576667, "WCMA":-72.506761, "VT":-72.671667, "ME":-69.233333} 

df_thermal_plants_existing_new = pd.concat([df_thermal_plants_existing]*e_num, ignore_index=True)
df_thermal_plants_existing_gen = pd.read_csv(f"{result_dir}/power_output_existing_generators.csv", header=None)

p_gen_existing = pd.concat([df_thermal_plants_existing_new,df_thermal_plants_existing_gen], axis=1)
p_gen_existing["e_num"] = pd.Series([j+1 for j in range(e_num) for i in range(0,n_generators)])

p_gen_existing['unit'] = p_gen_existing["status"] + " " + p_gen_existing["tech"] + " " + p_gen_existing["index"].astype(str)

# Mapping of desired variable names to their respective file names
file_mapping = {
    "p_gen_ng_cc_ccs": "power_output_new_ng_cc_ccs.csv",
    "p_gen_ng_ct": "power_output_new_ng_ct.csv",
    "p_gen_pv": "pv_node_out.csv",
    "p_gen_wind": "wind_node_out.csv",
    "p_gen_pv_new": "pv_node_new_out.csv",
    "p_gen_wind_new": "wind_node_new_out.csv",
    "ls": "load_curtailments.csv",
    "ws": "wind_spillage.csv",
    "ch_node": "ch_node.csv",
    "dis_node": "dis_node.csv",
    "bus_out_power": "bus_out_power.csv",
    "load_node_out": "load_node_out.csv",
    "flexible_load": "flexible_load.csv"
}

# Template DataFrame for structure
template_file = "power_output_new_ng_cc_ccs.csv"  # Assuming this file exists for structure reference

# Check if the template file exists in the directory
template_path = f"{result_dir}/{template_file}"
if not glob.glob(template_path):
    raise FileNotFoundError(f"Template file {template_file} not found in {result_dir}.")

# Read the columns from a known file to use as a template
template_df = pd.read_csv(template_path, header=None)

# Read files and assign to dataframes
dfs = {}
for var_name, file in file_mapping.items():
    file_path = f"{result_dir}/{file}"
    try:
        # Attempt to read the CSV file
        dfs[var_name] = pd.read_csv(file_path, header=None)
    except FileNotFoundError:
        # If the file is not found, log and create an empty DataFrame using the template
        print(f"File not found: {file}. Creating empty DataFrame for {var_name}.")
        dfs[var_name] = pd.DataFrame(columns=template_df.columns)
    except Exception as e:
        # Log any other exception that occurs
        print(f"Error reading {file}: {e}")
        dfs[var_name] = pd.DataFrame(columns=template_df.columns)

# Dynamically create variables from the keys in dfs
# for var_name, df in dfs.items():
#     # Use exec to create dynamic variables - Note: Generally not recommended
#     exec(f"{var_name} = df")

# Manually unpack dataframes from the dictionary into variables
p_gen_ng_cc_ccs = dfs.get("p_gen_ng_cc_ccs", pd.DataFrame())
p_gen_ng_ct = dfs.get("p_gen_ng_ct", pd.DataFrame())
p_gen_pv = dfs.get("p_gen_pv", pd.DataFrame())
p_gen_wind = dfs.get("p_gen_wind", pd.DataFrame())
p_gen_pv_new = dfs.get("p_gen_pv_new", pd.DataFrame())
p_gen_wind_new = dfs.get("p_gen_wind_new", pd.DataFrame())
ls = dfs.get("ls", pd.DataFrame())
ws = dfs.get("ws", pd.DataFrame())
ch_node = dfs.get("ch_node", pd.DataFrame())
dis_node = dfs.get("dis_node", pd.DataFrame())
bus_out_power = dfs.get("bus_out_power", pd.DataFrame())
load_node_out = dfs.get("load_node_out", pd.DataFrame())
flexible_load = dfs.get("flexible_load", pd.DataFrame())

# CG: Put wind spillage, demand response, charging, lost load as negative generation
ws = -ws
ls = -ls
flexible_load = -abs(flexible_load) # Accouting both direction's contribution of demand response as generation
ch_node = -ch_node

# Define a helper function for column assignment
def assign_cols(df, tech, status, fuel, e_num, node):
    if not df.empty:
        df["tech"] = tech
        df["status"] = status
        df["fuel"] = fuel
        df["e_num"] = pd.Series([j+1 for j in range(e_num) for i in range(n_buses)])
        df["node"] = pd.Series([i+1 for j in range(e_num) for i in range(n_buses)])
        df["zone"] = df["node"].map(node_zone_dict).astype(str).replace(zone_fix_dict)
        df["latitude"] = df["zone"].map(zone_latitude_dict)
        df["longitude"] = df["zone"].map(zone_longitude_dict)
        df["unit"] = df["status"] + " " + df["tech"] + " " + df["node"].astype(str)
        df["ave_marg_damages_isrm_LePeule"] = 0
        # df["index_existing"] = 0

# Assign columns to dataframes
assign_cols(p_gen_ng_cc_ccs, "NG-CC", "New", "Gas", e_num, n_buses)
assign_cols(p_gen_ng_ct, "NG-CT", "New", "Gas", e_num, n_buses)
assign_cols(p_gen_pv, "Solar PV", "Existing", "Solar PV", e_num, n_buses)
assign_cols(p_gen_wind, "Wind", "Existing", "Wind", e_num, n_buses)
assign_cols(p_gen_pv_new, "Solar PV", "New", "Solar PV", e_num, n_buses)
assign_cols(p_gen_wind_new, "Wind", "New", "Wind", e_num, n_buses)
assign_cols(ls, "Lost Load", "", "", e_num, n_buses)
assign_cols(ws, "Renewables Curtailment", "", "", e_num, n_buses)
assign_cols(ch_node, "Battery Charging", "New", "", e_num, n_buses)
assign_cols(dis_node, "Battery Discharging", "New", "Battery", e_num, n_buses)
assign_cols(bus_out_power, "Flow Out", "", "", e_num, n_buses)
assign_cols(load_node_out, "Load", "", "", e_num, n_buses)
assign_cols(flexible_load, "Demand Response", "New", "Demand Response", e_num, n_buses)


p_gen_ng_cc_ccs["ave_marg_damages_isrm_LePeule"] = 22.415 * 1.08
p_gen_ng_ct["ave_marg_damages_isrm_LePeule"] = 26.22 * 1.08


final_cols = list(bus_out_power.columns)
final_cols.append('index_existing')
p_gen_existing = p_gen_existing[final_cols]

p_gen = pd.concat([p_gen_existing, p_gen_ng_cc_ccs, p_gen_ng_ct, p_gen_pv, p_gen_wind, p_gen_pv_new, p_gen_wind_new, ls, ws, ch_node, dis_node, bus_out_power, load_node_out, flexible_load], axis=0)

# Replace this with the list of column names you want to pivot, e.g., ['1', '2', '3', '4']
columns_to_pivot = list(pd.Series([j for j in range(t_num * epoch)]))

# Create a list of column names to keep as identifiers
id_vars = ["index_existing", "unit", "tech","status","fuel","e_num","node","zone","latitude","longitude", "ave_marg_damages_isrm_LePeule"]

# Pivot the DataFrame
p_gen = p_gen.melt(id_vars=id_vars, value_vars=columns_to_pivot, var_name='Hour', value_name='Generation (MWh) Unadjusted')
p_gen['Epoch'] = np.ceil((p_gen['Hour'] + 1) / t_num)

p_gen['Generation (MWh)'] = p_gen.apply(lambda row: row['Generation (MWh) Unadjusted'] * tau_e[int(row['e_num']) - 1] * epoch_length, axis=1)

# Output generation data for power balance equation
p_gen.to_csv(f"{result_dir}/Pg.csv", index = False, date_format="%m/%d/%Y")


# List of columns you want to remove
columns_to_remove = ['index_existing', 'unit', 'latitude', 'longitude']

# Drop the unwanted columns
p_gen_reduced = p_gen.drop(columns=columns_to_remove)

columns_to_groupby = list(set(p_gen.columns) - set(columns_to_remove) - set(['Generation (MWh) Unadjusted','Generation (MWh)']))

# Now group by 'e_num' and sum all the remaining numeric columns
p_gen_reduced = p_gen_reduced.groupby(columns_to_groupby, as_index=False).sum()
p_gen_reduced["Specification"] = os.path.basename(result_dir)

p_gen_reduced['Generation (TWh)'] = p_gen_reduced['Generation (MWh)']/1e6

p_gen_reduced.to_csv(f"{result_dir}/Pg_reduced.csv", index = False, date_format="%m/%d/%Y")


p_gen_aq = pd.concat([p_gen_existing, p_gen_ng_cc_ccs, p_gen_ng_ct], axis=0)

# Pivot the DataFrame
p_gen_aq = p_gen_aq.melt(id_vars=id_vars, value_vars=columns_to_pivot, var_name='Hour', value_name='Generation (MWh) Unadjusted')
p_gen_aq['Epoch'] = np.ceil((p_gen_aq['Hour'] + 1) / t_num)

p_gen_aq['Generation (MWh)'] = p_gen_aq.apply(lambda row: row['Generation (MWh) Unadjusted'] * tau_e[int(row['e_num']) - 1] * epoch_length, axis=1)

epoch2year_start = {1: 1, 2: 6, 3: 11, 4: 16}
epoch2year_end = {1: 5, 2: 10, 3: 15, 4: 20}

p_gen_aq['aq_costs'] = p_gen_aq.apply(lambda row: 1 / ((1 + float(params_dict["r"])) ** (epoch2year_end[row['Epoch']] - 1)) * row['Generation (MWh)'] * row["ave_marg_damages_isrm_LePeule"], axis=1)

columns_to_remove = ['index_existing','Epoch', 'Hour', 'e_num']
p_gen_aq = p_gen_aq.drop(columns=columns_to_remove)

columns_to_groupby = list(set(p_gen_aq.columns) - set(columns_to_remove) - set(['Generation (MWh) Unadjusted','Generation (MWh)', 'aq_costs']))
p_gen_aq = p_gen_aq.groupby(columns_to_groupby, as_index=False).sum()

p_gen_aq = p_gen_aq[p_gen_aq['aq_costs'] != 0]

p_gen_aq["Specification"] = os.path.basename(result_dir)


p_gen_aq.to_csv(f"{result_dir}/Pg_aq.csv", index = False, date_format="%m/%d/%Y")




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 17:17:04 2023

@author: skhanal
"""

import pandas as pd
# from pathlib import Path
import os

# Parameters
n_lines_existing = 12 # 13 for PJM
n_buses_existing = 8 # 9 for PJM
epoch = 4
t_num = 24
epoch_length = 5 # 5 years represent an epoch

### CHANGE THESE BEFORE RUNS ###
result_dir = "/Users/skhanal/Globus_local/MOTEP-OSW/outputs/paper_draft/"
data_dir = "/Users/skhanal/Globus_local/MOTEP-OSW/inputs/ISONE8busTS_osw1_17_N5X5/"
base_case_dirname = "SO" # Base on which differentials should be computed

n_buses = len(pd.read_csv(f"{data_dir}/nodes.csv")) # SOline: 8, SO = 15, Optimized POI = 14
# n_lines = len(pd.read_csv(f"{data_dir}/lines.csv"))
# e_num = int(len(pd.read_csv(f"{data_dir}/wind_power_nodes.csv", header=None))/t_num)
# n_lines_new = n_lines - n_lines_existing

zone_node_dict = pd.read_csv(f"{data_dir}nodes.csv").set_index("Zone")["node"].to_dict()
node_zone_dict = pd.read_csv(f"{data_dir}nodes.csv").set_index("node")["Zone"].to_dict()

# Mannual ovveride to insconsistent zone names
zone_fix_dict = {"SEMASS" : "SEMA" , "WCMASS" : "WCMA", "NEMA/BOST": "NEMA", "CTREV":"REV", "RIREV":"REV"}

# Annual costs
cost_components = {
    'BESS Investment': 'annual_cost_bess_investment.csv',
    'Line Investment': 'annual_cost_line_investment.csv',
    'Operation': 'annual_cost_operation.csv',
    'Air Quality': 'annual_cost_externalities_airquality.csv',
    'Emission': 'annual_cost_externalities_emission.csv',
    'Generation Investment': 'annual_cost_gen_investment.csv'
}

df_collect_cost = []
for file in [f for f in os.listdir(result_dir) if os.path.isdir(os.path.join(result_dir, f))]: 
    for component, filename in cost_components.items():
        try:
            cost_df = pd.read_csv(os.path.join(result_dir, file, filename), header=None, names=["Cost ($)"])
            cost_df["Cost ($)"] = cost_df["Cost ($)"] * epoch_length
            cost_df["Scenario"] = file
            cost_df["Epoch"] = pd.Series([j+1 for j in range(epoch)])
            cost_df["Cost Components"] = component
            df_collect_cost.append(cost_df)
        except FileNotFoundError:
            print(f"File not found: {os.path.join(result_dir, file, filename)}")

df_costs = pd.concat(df_collect_cost)
df_costs = df_costs[['Scenario', 'Cost Components', 'Epoch', "Cost ($)"]]

def calculate_percentage_difference(df, cost_field):

    first_epoch_cost = df.loc[df['Epoch'] == df['Epoch'].min(), cost_field].iloc[0]
    if round(first_epoch_cost) == 0:
        first_epoch_cost = 0

    # Calculate the percentage difference with the locally modified cost field, handling division by zero
    def compute_difference(cost, first_epoch_cost):
        if abs(cost) == 0:
            cost = 0
        if first_epoch_cost == 0:
            return 0
        else:
            return ((cost - first_epoch_cost) / first_epoch_cost) * 100

    df[f'{cost_field} Percentage Difference (%)'] = df[cost_field].apply(lambda cost: compute_difference(cost, first_epoch_cost)).fillna(0)

    return df
df_costs = df_costs.groupby(['Scenario', 'Cost Components']).apply(lambda group: calculate_percentage_difference(group, 'Cost ($)')).reset_index(drop=True)
df_costs.to_csv(f"/../{result_dir}/costs.csv", index = False, date_format="%m/%d/%Y")

def process_dataframe_collapseddiff(df, base_scenario, value_column, group_columns):
    # Group by group_columns, then sum the value_column
    df_sum = df.groupby(group_columns + ['Epoch'])[value_column].sum().reset_index()

    # Collapse 'Epoch' column
    df_sum = df_sum.groupby(group_columns)[value_column].sum().reset_index()

    # Calculate the base scenario values
    df_base = df_sum[df_sum[group_columns[0]] == base_scenario].set_index(group_columns[1])[value_column]

    # Calculate percentage difference for each scenario's value compared to the base scenario
    percentage_column = f"Percentage Difference (%) ({value_column})"
    df_sum[percentage_column] = df_sum.apply(lambda row: ((row[value_column] - df_base.loc[row[group_columns[1]]]) / df_base.loc[row[group_columns[1]]]) * 100, axis=1).fillna(0)

    # Reset index and return the final dataframe
    df_sum = df_sum.reset_index(drop=True)
    return df_sum

process_dataframe_collapseddiff(df_costs, base_case_dirname, "Cost ($)", ['Scenario', 'Cost Components']).to_csv(f"/../{result_dir}/costs_diff.csv", index = False, date_format="%m/%d/%Y")

# Operations cost breakdown
df_collect_cost_operation = []
for file in [f for f in os.listdir(result_dir) if os.path.isdir(os.path.join(result_dir, f))]: 
    try:
        annual_cost_operation_demand_flexibility = pd.read_csv(f"{result_dir}{file}/annual_cost_operation_demand_flexibility.csv", header=None, names=["Cost ($)"])
        annual_cost_operation_demand_flexibility["Cost ($)"] = annual_cost_operation_demand_flexibility["Cost ($)"] * epoch_length
        annual_cost_operation_demand_flexibility["Scenario"] = file
        annual_cost_operation_demand_flexibility["Epoch"] = pd.Series([j+1 for j in range(epoch)])
        annual_cost_operation_demand_flexibility["Cost Components"] = "Demand Flexibility"
        df_collect_cost_operation.append(annual_cost_operation_demand_flexibility)
    except Exception as e:
        print(f"An error occurred: {e}")
        annual_cost_operation_demand_flexibility["Cost ($)"] = 0
        annual_cost_operation_demand_flexibility["Scenario"] = file
        annual_cost_operation_demand_flexibility["Epoch"] = pd.Series([j+1 for j in range(epoch)])
        annual_cost_operation_demand_flexibility["Cost Components"] = "Demand Flexibility"
        df_collect_cost_operation.append(annual_cost_operation_demand_flexibility)

for file in [f for f in os.listdir(result_dir) if os.path.isdir(os.path.join(result_dir, f))]: 
    try:
        annual_cost_operation_over_generation = pd.read_csv(f"{result_dir}{file}/annual_cost_operation_over_generation.csv", header=None, names=["Cost ($)"])
        annual_cost_operation_over_generation["Cost ($)"] = annual_cost_operation_over_generation["Cost ($)"] * epoch_length
        annual_cost_operation_over_generation["Scenario"] = file
        annual_cost_operation_over_generation["Epoch"] = pd.Series([j+1 for j in range(epoch)])
        annual_cost_operation_over_generation["Cost Components"] = "Over Generation Penalty"
        df_collect_cost_operation.append(annual_cost_operation_over_generation)
    except Exception as e:
        print(f"An error occurred: {e}")
        annual_cost_operation_over_generation["Cost ($)"] = 0
        annual_cost_operation_over_generation["Scenario"] = file
        annual_cost_operation_over_generation["Epoch"] = pd.Series([j+1 for j in range(epoch)])
        annual_cost_operation_over_generation["Cost Components"] = "Over Generation Penalty"
        df_collect_cost_operation.append(annual_cost_operation_over_generation)

for file in [f for f in os.listdir(result_dir) if os.path.isdir(os.path.join(result_dir, f))]: 
    annual_cost_operation_under_generation = pd.read_csv(f"{result_dir}{file}/annual_cost_operation_under_generation.csv", header=None, names=["Cost ($)"])
    annual_cost_operation_under_generation["Cost ($)"] = annual_cost_operation_under_generation["Cost ($)"] * epoch_length
    annual_cost_operation_under_generation["Scenario"] = file
    annual_cost_operation_under_generation["Epoch"] = pd.Series([j+1 for j in range(epoch)])
    annual_cost_operation_under_generation["Cost Components"] = "Under Generation Penalty"
    df_collect_cost_operation.append(annual_cost_operation_under_generation)

for file in [f for f in os.listdir(result_dir) if os.path.isdir(os.path.join(result_dir, f))]: 
    annual_cost_operation_rps_noncompliance = pd.read_csv(f"{result_dir}{file}/annual_cost_operation_rps_noncompliance.csv", header=None, names=["Cost ($)"])
    annual_cost_operation_rps_noncompliance["Cost ($)"] = annual_cost_operation_rps_noncompliance["Cost ($)"] * epoch_length
    annual_cost_operation_rps_noncompliance["Scenario"] = file
    annual_cost_operation_rps_noncompliance["Epoch"] = pd.Series([j+1 for j in range(epoch)])
    annual_cost_operation_rps_noncompliance["Cost Components"] = "RPS Noncompliance"
    df_collect_cost_operation.append(annual_cost_operation_rps_noncompliance)

for file in [f for f in os.listdir(result_dir) if os.path.isdir(os.path.join(result_dir, f))]: 
    annual_cost_operation_demand_flexibility = pd.read_csv(f"{result_dir}{file}/annual_cost_operation_demand_flexibility.csv", header=None, names=["Cost ($)"])
    annual_cost_operation_over_generation = pd.read_csv(f"{result_dir}{file}/annual_cost_operation_over_generation.csv", header=None, names=["Cost ($)"])
    annual_cost_operation_under_generation = pd.read_csv(f"{result_dir}{file}/annual_cost_operation_under_generation.csv", header=None, names=["Cost ($)"])
    annual_cost_operation_rps_noncompliance = pd.read_csv(f"{result_dir}{file}/annual_cost_operation_rps_noncompliance.csv", header=None, names=["Cost ($)"])
    annual_cost_operation = pd.read_csv(f"{result_dir}{file}/annual_cost_operation.csv", header=None, names=["Cost ($)"])
    annual_cost_operation["Cost ($)"] = annual_cost_operation["Cost ($)"] - annual_cost_operation_demand_flexibility["Cost ($)"] - annual_cost_operation_over_generation["Cost ($)"] - annual_cost_operation_under_generation["Cost ($)"] - annual_cost_operation_rps_noncompliance["Cost ($)"]
    annual_cost_operation["Cost ($)"] = annual_cost_operation["Cost ($)"] * epoch_length 
    annual_cost_operation["Scenario"] = file
    annual_cost_operation["Epoch"] = pd.Series([j+1 for j in range(epoch)])
    annual_cost_operation["Cost Components"] = "Generator O&M"
    df_collect_cost_operation.append(annual_cost_operation)

df_costs_operation = pd.concat(df_collect_cost_operation)
df_costs_operation = df_costs_operation[['Scenario', 'Cost Components', 'Epoch', "Cost ($)"]]
df_costs_operation = df_costs_operation.groupby(['Scenario', 'Cost Components']).apply(lambda group: calculate_percentage_difference(group, 'Cost ($)')).reset_index(drop=True)
df_costs_operation = df_costs_operation
df_costs_operation.to_csv(f"/../{result_dir}costs_operation.csv", index = False, date_format="%m/%d/%Y")
process_dataframe_collapseddiff(df_costs_operation, base_case_dirname, "Cost ($)", ['Scenario', 'Cost Components']).to_csv(f"/../{result_dir}/costs_operation_diff.csv", index = False, date_format="%m/%d/%Y")

# Annual capacity additions
df_collect_capacity = []
# Create a list of technology names and corresponding filenames
tech_names = [("Battery", "new_build_gen_cap_bess_power.csv"),
              ("NG-CC-CCS", "new_build_gen_cap_ng_cc_ccs.csv"),
              ("NG-CT", "new_build_gen_cap_ng_ct.csv"),
              ("Wind", "new_build_gen_cap_wind.csv"),
              ("Solar PV", "new_build_gen_cap_pv.csv")]

# Loop through the technology names and filenames to extract data
for tech_name, filename in tech_names:
    for file in [f for f in os.listdir(result_dir) if os.path.isdir(os.path.join(result_dir, f))]:
        new_build_gen_cap = pd.read_csv(f"{result_dir}{file}/{filename}", header=None, names=[j+1 for j in range(epoch)])
        new_build_gen_cap = new_build_gen_cap#.head(n_buses_existing)
        new_build_gen_cap["Scenario"] = file
        new_build_gen_cap["Zone Number"] = pd.Series([j+1 for j in range(n_buses)])#range(n_buses_existing)
        new_build_gen_cap["Zone"] = new_build_gen_cap["Zone Number"].map(node_zone_dict).replace(zone_fix_dict)
        new_build_gen_cap["Technology"] = tech_name
        df_collect_capacity.append(new_build_gen_cap)

df_capacity = pd.concat(df_collect_capacity)
df_capacity = df_capacity[["Scenario", "Zone Number", "Zone"] + [j+1 for j in range(epoch)] + ["Technology"]]

# Replace this with the list of column names you want to pivot, e.g., ['1', '2', '3', '4']
columns_to_pivot = list(pd.Series([j+1 for j in range(epoch)]))

# Create a list of column names to keep as identifiers
id_vars = ['Scenario', 'Zone Number', 'Zone', 'Technology']

# Pivot the DataFrame
df_capacity = df_capacity.melt(id_vars=id_vars, value_vars=columns_to_pivot, var_name='Epoch', value_name='Capacity (MW)')
df_capacity = df_capacity.sort_values(by=['Scenario', 'Zone Number', 'Zone', 'Technology', 'Epoch']).reset_index(drop=True)
df_capacity.loc[(df_capacity['Zone Number'] > n_buses_existing) & (df_capacity['Technology'] == 'Battery'), 'Technology'] = 'Offshore Storage'

df_capacity.to_csv(f"/../{result_dir}capacity_additions.csv", index = False, date_format="%m/%d/%Y")
process_dataframe_collapseddiff(df_capacity, base_case_dirname, "Capacity (MW)", ['Scenario', 'Technology']).to_csv(f"/../{result_dir}/capacity_additions_diff.csv", index = False, date_format="%m/%d/%Y")

# Annual emissions in nominal units
emission_components = {
    'SO2' : 'Emissions_SO2.csv',
    'CO2' : 'Emissions_CO2.csv',
    'NOx' : 'Emissions_NOx.csv',
    'PM25' : 'Emissions_PM25.csv'
}

df_collect_emission = []
for file in [f for f in os.listdir(result_dir) if os.path.isdir(os.path.join(result_dir, f))]: 
    for component, filename in emission_components.items():
        try:
            emission_df = pd.read_csv(os.path.join(result_dir, file, filename), header=None, names=["Emission (mt)"])
            emission_df["Scenario"] = file
            emission_df["Epoch"] = pd.Series([j+1 for j in range(epoch)])
            emission_df["Emission Components"] = component
            df_collect_emission.append(emission_df)
        except FileNotFoundError:
            print(f"File not found: {os.path.join(result_dir, file, filename)}")

df_emissions = pd.concat(df_collect_emission)
df_emissions = df_emissions[['Scenario', 'Emission Components', 'Epoch', "Emission (mt)"]]
df_emissions = df_emissions.groupby(['Scenario', 'Emission Components']).apply(lambda group: calculate_percentage_difference(group, "Emission (mt)")).reset_index(drop=True)
df_emissions.to_csv(f"/../{result_dir}/emissions.csv", index = False, date_format="%m/%d/%Y")
process_dataframe_collapseddiff(df_emissions, base_case_dirname, "Emission (mt)", ['Scenario', 'Emission Components']).to_csv(f"/../{result_dir}/emissions_diff.csv", index = False, date_format="%m/%d/%Y")

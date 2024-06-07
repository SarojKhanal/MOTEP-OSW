import shutil
from io import BytesIO, TextIOWrapper
from zipfile import ZipFile
import urllib.request
import ssl
import certifi
import csv
from shapely.geometry import Point
from shapely.ops import nearest_points
import geopandas as gpd
import pandas as pd
import numpy as np
from pathlib import Path
import time

# This allows us to use the 'run_sr' function in the 'sr_util.py' file in this same directory.
# from src.utils.inmap.sr_util import run_sr, _tmpdir
from sr_util import run_sr, _tmpdir
import matplotlib.pyplot as plt
import requests


#############################################################
# PARAMETERS                                                #
#############################################################
# https://github.com/USEPA/cam-api-examples/blob/main/Python/account_data_demo.py
API_KEY = pd.read_table(Path("/mnt/c/users/c/keys/data_gov.txt")).columns[0]
use_existing_camd_data = True  # file ./data/epa_camd.csv must exist

#############################################################
# FUNCTIONS                                                 #
#############################################################

# api parameters for the streaming account holdings endpoint
parameters = {
    "api_key": API_KEY,
}


def get_power_plant_annual_camd_emissions(facility_id, year, parameters):
    # aggregate (over units) gross load [MWh], so2Mass [short ton], co2Mass [short ton], noxMass [short ton], heatInput [short ton]
    # output: 'grossLoad', 'so2Mass', 'co2Mass', 'noxMass', 'heatInput'
    page_number = 1  # unit
    error = False
    results_df = pd.DataFrame()
    while error is False:
        url = (
            "https://api.epa.gov/easey/emissions-mgmt/emissions/apportioned/annual?facilityId="
            + str(facility_id)
            + "&year="
            + str(year)
            + "&page="
            + str(page_number)
            + "&perPage=1"
        )
        results = requests.get(url, params=parameters)
        temp = pd.DataFrame(results.json())
        time.sleep(6.0)
        if len(temp) > 0:
            results_df = pd.concat([results_df, temp])
            page_number += 1
        else:
            error = True

    return results_df


def get_power_plant_emissions(url):
    VOC, NOx, NH3, SOx, PM2_5 = [], [], [], [], []
    height, diam, temp, velocity = [], [], [], []
    coords, fips, facility, facility_id, oris_facility_id, oris_boiler_id = (
        [],
        [],
        [],
        [],
        [],
        [],
    )

    def add_record(row):
        """Process one row of the emissions file"""
        pol = row[
            12
        ]  # The pollutant is in the 13th column of the CSV file (In Python, the first column is called column 0.)
        # We are only extracting annual total emissions here. If monthly emissions are reported, we'll miss them. Emissions are short tons/year.
        emis = row[13]
        if emis == "":
            return
        if pol in [
            "VOC",
            "VOC_INV",
            "XYL",
            "TOL",
            "TERP",
            "PAR",
            "OLE",
            "NVOL",
            "MEOH",
            "ISOP",
            "IOLE",
            "FORM",
            "ETOH",
            "ETHA",
            "ETH",
            "ALD2",
            "ALDX",
            "CB05_ALD2",
            "CB05_ALDX",
            "CB05_BENZENE",
            "CB05_ETH",
            "CB05_ETHA",
            "CB05_ETOH",
            "CB05_FORM",
            "CB05_IOLE",
            "CB05_ISOP",
            "CB05_MEOH",
            "CB05_OLE",
            "CB05_PAR",
            "CB05_TERP",
            "CB05_TOL",
            "CB05_XYL",
            "ETHANOL",
            "NHTOG",
            "NMOG",
            "VOC_INV",
        ]:
            VOC.append(float(emis))
            NOx.append(0)
            NH3.append(0)
            SOx.append(0)
            PM2_5.append(0)
        elif pol in [
            "PM25-PRI",
            "PM2_5",
            "DIESEL-PM25",
            "PAL",
            "PCA",
            "PCL",
            "PEC",
            "PFE",
            "PK",
            "PMG",
            "PMN",
            "PMOTHR",
            "PNH4",
            "PNO3",
            "POC",
            "PSI",
            "PSO4",
            "PTI",
        ]:
            VOC.append(0)
            NOx.append(0)
            NH3.append(0)
            SOx.append(0)
            PM2_5.append(float(emis))
        elif pol in ["NOX", "HONO", "NO", "NO2"]:
            VOC.append(0)
            NOx.append(float(emis))
            NH3.append(0)
            SOx.append(0)
            PM2_5.append(0)
        elif pol == "NH3":
            VOC.append(0)
            NOx.append(0)
            NH3.append(float(emis))
            SOx.append(0)
            PM2_5.append(0)
        elif pol == "SO2":
            VOC.append(0)
            NOx.append(0)
            NH3.append(0)
            SOx.append(float(emis))
            PM2_5.append(0)
        else:
            return

        h = row[17]
        height.append(float(h) * 0.3048) if h != "" else height.append(0)

        d = row[18]
        diam.append(float(d) * 0.3048) if d != "" else diam.append(0)

        t = row[19]
        temp.append((float(t) - 32) * 5.0 / 9.0 + 273.15) if t != "" else temp.append(0)

        v = row[21]
        velocity.append(float(v) * 0.3048) if v != "" else velocity.append(0)

        coords.append(Point(float(row[23]), float(row[24])))

        # facility_name, id, unit id
        facility.append(str(row[15]))
        facility_id.append(int(row[3]))
        oris_facility_id.append(str(row[41]))
        oris_boiler_id.append(str(row[42]))

        # region_cd
        fips.append(int(row[1]))

    with ZipFile(BytesIO(url.read())) as zf:
        for contained_file in zf.namelist():
            # Only process files with electricity generating unit (EGU) emissions.
            if "egu" in contained_file:
                for row in csv.reader(
                    TextIOWrapper(zf.open(contained_file, "r"), newline="")
                ):
                    if (len(row) == 0) or (len(row[0]) == 0) or (row[0][0] == "#"):
                        continue
                    add_record(row)

    emis = gpd.GeoDataFrame(
        {
            "FACILITY_ID": facility_id,
            "FACILITY": facility,
            "ORIS_FACILITY_ID": oris_facility_id,
            "ORIS_BOILER_ID": oris_boiler_id,
            "FIPS": fips,
            "VOC": VOC,
            "NOx": NOx,
            "NH3": NH3,
            "SOx": SOx,
            "PM2_5": PM2_5,
            "height": height,
            "diam": diam,
            "temp": temp,
            "velocity": velocity,
        },
        geometry=coords,
        crs="EPSG:4269",
    )

    return emis


def get_fips():
    # FIPS (Federal Information Processing Standards) codes. There are state codes and county codes: the 2016 state and county FIPS codes can be found at the US Census Website.
    # https://plotly.com/python/county-choropleth/
    dt = pd.read_csv(
        "https://raw.githubusercontent.com/plotly/datasets/master/minoritymajority.csv"
    )
    # dt = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/laucnty16.csv')
    # dt['State FIPS Code'] = dt['State FIPS Code'].apply(lambda x: str(x).zfill(2))
    # dt['County FIPS Code'] = dt['County FIPS Code'].apply(lambda x: str(x).zfill(3))
    # dt['FIPS'] = dt['State FIPS Code'] + dt['County FIPS Code']
    return dt


def get_pollution_deaths_per_year(emis, model="isrm"):
    # Source-receptor (SR) matrices
    # model="isrm", "apsca_q0"
    output_variables = {
        "TotalPM25": "PrimaryPM25 + pNH4 + pSO4 + pNO3 + SOA",
        "deathsK": "(exp(log(1.06)/10 * TotalPM25) - 1) * TotalPop * 1.0465819687408728 * MortalityRate / 100000 * 1.025229357798165",
        "deathsL": "(exp(log(1.14)/10 * TotalPM25) - 1) * TotalPop * 1.0465819687408728 * MortalityRate / 100000 * 1.025229357798165",
    }

    results = run_sr(
        emis, model, emis_units="tons/year", output_variables=output_variables
    )

    return results


def summarize_results(results, model="isrm", vsl=9.0e6):
    deaths = pd.DataFrame.from_dict(
        {
            "Model": [model],
            "Krewski Deaths": [results.deathsK.sum()],
            "LePeule Deaths": [results.deathsL.sum()],
        }
    )

    deaths_dollar = pd.DataFrame.from_dict(
        {
            "Model": [model],
            "Krewski Damages": deaths["Krewski Deaths"] * vsl,
            "LePeule Damages": deaths["LePeule Deaths"] * vsl,
        }
    )

    return [deaths, deaths_dollar]


def plot_us_map(emis):
    pols = ["SOx", "NOx", "PM2_5", "VOC", "NH3"]
    pol_names = ["SO$_2$", "NO$_x$", "PM$_{2.5}$", "VOC", "NH$_3$"]

    fig, axes = plt.subplots(figsize=(7, 3), nrows=2, ncols=3, sharex=True, sharey=True)
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0.1, hspace=0.1)

    i = 0
    for x in axes:
        for ax in x:
            if i < len(pols):
                emis.plot(ax=ax, markersize=emis[pols[i]] ** 0.5 / 5)
                ax.set_title(pol_names[i])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.axis("off")
            i = i + 1
    plt.show()


def match_epa_eia_plants(eia, epa_eia_crosswalk):
    dt = eia.copy()

    # merge eia data with epa data using crosswalk
    epa_eia_crosswalk = epa_eia_crosswalk.loc[~epa_eia_crosswalk["EIA_PLANT_ID"].isna()]

    epa_eia_crosswalk = epa_eia_crosswalk.astype(
        {"EIA_PLANT_ID": int, "CAMD_PLANT_ID": int}
    )

    epa_eia_crosswalk = epa_eia_crosswalk.loc[
        epa_eia_crosswalk["EIA_PLANT_ID"].isin(eia["EIA_Plant_Code"].unique())
    ]
    epa_eia_crosswalk = epa_eia_crosswalk.loc[
        epa_eia_crosswalk["EIA_PLANT_ID"].isin(eia["EIA_Plant_Code"].unique())
    ]
    epa_eia_crosswalk = epa_eia_crosswalk[
        ["EIA_PLANT_ID", "CAMD_PLANT_ID", "CAMD_FUEL_TYPE"]
    ].drop_duplicates()

    dt["EIA_FUEL_TYPE"] = np.nan
    for j in dt["EIA_Type"].unique():
        dt.loc[dt["EIA_Type"] == j, "EIA_FUEL_TYPE"] = j.split("_")[0]

    epa_eia_crosswalk["CAMD_FUEL"] = np.nan
    for j in dt["EIA_FUEL_TYPE"].unique():
        epa_eia_crosswalk.loc[
            epa_eia_crosswalk["CAMD_FUEL_TYPE"].str.find(j) > -1, "CAMD_FUEL"
        ] = j

    epa_eia_crosswalk = epa_eia_crosswalk[
        ["EIA_PLANT_ID", "CAMD_PLANT_ID", "CAMD_FUEL"]
    ].drop_duplicates()

    epa_eia_crosswalk = epa_eia_crosswalk.dropna()

    dt = dt.merge(
        epa_eia_crosswalk[["EIA_PLANT_ID", "CAMD_PLANT_ID", "CAMD_FUEL"]],
        how="inner",
        left_on=["EIA_Plant_Code", "EIA_FUEL_TYPE"],
        right_on=["EIA_PLANT_ID", "CAMD_FUEL"],
    )
    dt = dt.drop("EIA_PLANT_ID", axis=1)

    return dt


def match_epa_emissions_to_eia_plants(emis, eia):
    # inner match between eia plant level data and emis based on EIA Plant Code and ORIS Facility ID
    # test = gpd.sjoin_nearest(left_df=emis_isone, right_df=eia, how='left')
    # test = test.merge(eia[['EIA_Plant_Code', 'EIA_Type', 'geometry']].rename({'geometry': 'EIA_geometry'}, axis=1), how='left', on=['EIA_Plant_Code', 'EIA_Type'])
    # https://github.com/USEPA/camd-eia-crosswalk
    # CAMD_PLANT_ID: 	The unique ID (also known as ORIS code/ORISPL code) of the facility in CAMD's data. (CAMD key)
    # https://www.cmascenter.org/smoke/documentation/4.5/html/ch08s02s08.html#sect_input_ptinv_ff10
    dt = eia.copy()

    emis = emis.loc[emis.ORIS_FACILITY_ID != ""]
    emis = emis.astype({"ORIS_FACILITY_ID": int})

    dt = dt.merge(
        emis, how="inner", left_on="EIA_Plant_Code", right_on="ORIS_FACILITY_ID"
    )
    dt = dt.rename(
        {"geometry_x": "EIA_geometry", "geometry_y": "ORIS_geometry"}, axis=1
    )

    return dt


def plot_marg_loc_emissions(
    dt,
    ave_marg_damages_col,
    fname="ave_marg_damages_map.png",
    x_lim=(-80, -66),
    y_lim=(38, 48),
    markerscaler=20.0,
):
    # us_counties_path = "https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_county_500k.zip"
    # us_counties = gpd.read_file(us_counties_path)

    us_states_path = (
        "https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_state_500k.zip"
    )
    us_states = gpd.read_file(us_states_path)

    # dt = pd.read_csv(Path(Path.cwd(), 'ave_marg_damages.csv'))
    gdt = gpd.GeoDataFrame(dt[~dt["Latitude"].isna()])
    gdt["geometry"] = gpd.GeoSeries.from_xy(
        x=gdt["Longitude"], y=gdt["Latitude"], crs=us_states.crs
    )
    # file:///C:/Users/c/Downloads/plant_power_eia_v6.pdf
    # dt = gpd.GeoDataFrame(dt, geometry="geometry", crs="EPSG:3857")
    # dt = gdt.to_crs('EPSG:4269') # us_states.crs
    gdt["markersize"] = (
        markerscaler * gdt[ave_marg_damages_col] / max(gdt[ave_marg_damages_col])
    )

    ax = us_states.boundary.plot(edgecolor="black", linewidth=0.2, figsize=(10, 5))
    # add ave_marg_damages for each power plant
    gdt.plot(ax=ax, marker="o", color="red", markersize=gdt.markersize)

    plt.tight_layout()
    plt.xlim(x_lim)
    plt.ylim(y_lim)
    plt.axis("off")
    plt.savefig(
        Path(Path.cwd(), "figures", fname),
        dpi=800,
        format="png",
        metadata=None,
        bbox_inches=None,
        pad_inches=0.1,
        facecolor="auto",
        edgecolor="auto",
        backend=None,
    )


def compute_marg_loc_emissions(emis_facility_boiler, model="isrm", vsl=9.0e6):
    # compute marginal local emissions for each ORIS_FACILITY_ID/ORIS_BOILER_ID
    # emis_facility_boiler shoud be a boiler specific emissions data from epa (emis file)
    resultsISRM = get_pollution_deaths_per_year(emis_facility_boiler, model="isrm")
    time.sleep(5)
    deaths, deaths_dollar = summarize_results(resultsISRM, model=model, vsl=vsl)
    annual_gen = emis_facility_boiler["grossLoad"].mean()
    annual_em = emis_facility_boiler[["VOC", "NOx", "NH3", "SOx", "PM2_5"]].sum()

    # increase emission by assume 1 MWh more grossLoad of facility/boiler
    temp = emis_facility_boiler.copy()

    for jj in annual_em.index:
        non_zero_elements = len(
            emis_facility_boiler.loc[(emis_facility_boiler[jj] > 0.0), jj]
        )

        if annual_em[jj] > 0:
            temp.loc[(temp[jj] > 0.0), jj] += (
                annual_em[jj] / annual_gen
            ) / non_zero_elements

    # re-run increasing emission by 1MWh of electric output more
    resultsISRM_marg_incr = get_pollution_deaths_per_year(temp, model=model)
    time.sleep(5)
    deaths_marg_incr, deaths_dollar_marg_incr = summarize_results(
        resultsISRM_marg_incr, model=model, vsl=vsl
    )

    out = (
        emis_facility_boiler[
            ["FACILITY_ID", "FACILITY", "ORIS_FACILITY_ID", "ORIS_BOILER_ID", "FIPS"]
        ]
        .drop_duplicates()
        .copy()
    )
    out["ave_marg_damages_" + model + "_" + "Krewski"] = 0.0
    out["ave_marg_damages_" + model + "_" + "LePeule"] = 0.0

    for damage_type in ["Krewski", "LePeule"]:
        out["ave_marg_damages_" + model + "_" + damage_type] = float(
            deaths_dollar_marg_incr.loc[
                deaths_dollar_marg_incr.Model == model, damage_type + " Damages"
            ]
            - deaths_dollar.loc[deaths_dollar.Model == model, damage_type + " Damages"]
        )
    return out


def main():
    # see https://www.inmap.run/blog/2019/04/20/sr/

    # Download file from EPA website.
    # url = urllib.request.urlopen("ftp://newftp.epa.gov/air/emismod/2016/alpha/2016fd/emissions/2016fd_inputs_point.zip")
    # or alternatively
    # url = urllib.request.urlopen("https://gaftp.epa.gov/Air/emismod/2016/alpha/2016fd/emissions/2016fd_inputs_point.zip")
    url = urllib.request.urlopen(
        "file://" + str(Path(Path.cwd(), "inputs", "2016fd_inputs_point.zip"))
    )

    emis = get_power_plant_emissions(url)
    emis = emis.loc[emis.ORIS_FACILITY_ID != ""]
    emis["ORIS_FACILITY_ID"] = emis["ORIS_FACILITY_ID"].astype(int)
    emis["ORIS_BOILER_ID_CLEANED"] = np.nan
    emis.loc[emis.ORIS_BOILER_ID.str.isnumeric(), "ORIS_BOILER_ID_CLEANED"] = emis.loc[
        emis.ORIS_BOILER_ID.str.isnumeric(), "ORIS_BOILER_ID"
    ]
    emis.loc[emis.ORIS_BOILER_ID.str.isalpha(), "ORIS_BOILER_ID_CLEANED"] = emis.loc[
        emis.ORIS_BOILER_ID.str.isalpha(), "ORIS_BOILER_ID"
    ]

    # load epa_resource data
    for syst_name, topo_name in zip(
        ["isone", "pjm"],
        [
            "ISONE3busTS_N5",
            "PJM9busTS_N5",
        ],
    ):
        dt = pd.read_csv(Path(Path.cwd(), "inputs", topo_name, "generators.csv"))
        emis_syst = emis.loc[
            emis.ORIS_FACILITY_ID.isin(np.unique(dt["ORIS Plant Code"]))
        ]

        if use_existing_camd_data:
            epa_dt = pd.read_csv(
                Path(Path.cwd(), "inputs", "epa_camd_" + syst_name + ".csv")
            )
        else:
            epa_dt = pd.DataFrame()
            for j in np.unique(dt["ORIS Plant Code"]):
                print(j)
                temp = get_power_plant_annual_camd_emissions(
                    facility_id=j, year=2016, parameters=parameters
                )
                epa_dt = pd.concat([epa_dt, temp])

            epa_dt = epa_dt.reset_index(drop=True)
            epa_dt.grossLoad = epa_dt.grossLoad.astype("float")
            epa_dt.to_csv(
                Path(Path.cwd(), "inputs", "epa_camd_" + syst_name + ".csv"),
                index=False,
            )

        # add grossLoad to emissions data
        epa_grossload = epa_dt[["facilityId", "unitId", "grossLoad"]]
        emis_syst = emis_syst.merge(
            epa_grossload.rename(
                {"facilityId": "ORIS_FACILITY_ID", "unitId": "ORIS_BOILER_ID"}, axis=1
            ),
            on=["ORIS_FACILITY_ID", "ORIS_BOILER_ID"],
            how="inner",
        )
        emis_syst = emis_syst.loc[
            (~emis_syst.grossLoad.isna()) & (emis_syst.grossLoad > 0.0)
        ]

        # compute local marginal emissions for each ORIS_FACILITY_ID / ORIS_BOILER_ID
        try:
            epa_margin_emis = pd.read_csv(
                Path(
                    Path.cwd(), "inputs", "epa_camd_marg_emission_" + syst_name + ".csv"
                )
            )
        except Exception as e:
            epa_margin_emis = pd.DataFrame(columns=["ORIS_FACILITY_ID"])

        for j1, j2 in (
            emis_syst.loc[
                ~emis_syst.ORIS_FACILITY_ID.isin(epa_margin_emis.ORIS_FACILITY_ID),
                ["ORIS_FACILITY_ID", "ORIS_BOILER_ID"],
            ]
            .drop_duplicates()
            .values
        ):
            print(j1, j2)
            temp = compute_marg_loc_emissions(
                emis_syst.loc[
                    (emis_syst["ORIS_FACILITY_ID"] == j1)
                    & (emis_syst["ORIS_BOILER_ID"] == j2)
                ],
                model="isrm",
                vsl=9.0e6,
            )
            epa_margin_emis = pd.concat([epa_margin_emis, temp])

            epa_margin_emis = epa_margin_emis.drop_duplicates()
            epa_margin_emis = epa_margin_emis.reset_index(drop=True)
            epa_margin_emis.to_csv(
                Path(
                    Path.cwd(), "inputs", "epa_camd_marg_emission_" + syst_name + ".csv"
                ),
                index=False,
            )

        # merge with generator database
        epa_margin_emis = epa_margin_emis.loc[
            epa_margin_emis["ave_marg_damages_isrm_Krewski"] != np.inf
        ]
        epa_eia_crosswalk = pd.read_csv(
            "https://github.com/USEPA/camd-eia-crosswalk/releases/download/v0.3/epa_eia_crosswalk.csv"
        )

        epa_margin_emis = epa_margin_emis.merge(
            epa_eia_crosswalk[
                ["CAMD_PLANT_ID", "CAMD_UNIT_ID", "CAMD_GENERATOR_ID"]
            ].drop_duplicates(),
            left_on=["ORIS_FACILITY_ID", "ORIS_BOILER_ID"],
            right_on=["CAMD_PLANT_ID", "CAMD_UNIT_ID"],
            how="left",
        )
        # add lon / lat
        emis_syst["lon"] = emis_syst["geometry"].x
        emis_syst["lat"] = emis_syst["geometry"].y
        epa_margin_emis = epa_margin_emis.merge(
            emis_syst[
                ["ORIS_FACILITY_ID", "ORIS_BOILER_ID", "lon", "lat"]
            ].drop_duplicates(),
            on=["ORIS_FACILITY_ID", "ORIS_BOILER_ID"],
            how="left",
        )

        dt[["ave_marg_damages_isrm_Krewski", "ave_marg_damages_isrm_LePeule"]] = np.nan
        for j in dt.index:
            (
                temp_plant,
                temp_unit,
            ) = dt.loc[j, ["ORIS Plant Code", "Unit ID"]].values
            if (
                len(
                    epa_margin_emis.loc[
                        (epa_margin_emis.CAMD_PLANT_ID == temp_plant)
                        & (epa_margin_emis.CAMD_GENERATOR_ID == temp_unit)
                    ]
                )
                > 1
            ):
                dt.loc[
                    j,
                    ["ave_marg_damages_isrm_Krewski", "ave_marg_damages_isrm_LePeule"],
                ] = epa_margin_emis.loc[
                    (epa_margin_emis.CAMD_PLANT_ID == temp_plant)
                    & (epa_margin_emis.CAMD_GENERATOR_ID == temp_unit),
                    ["ave_marg_damages_isrm_Krewski", "ave_marg_damages_isrm_LePeule"],
                ].values.mean(
                    0
                )
                dt.loc[j, ["Latitude", "Longitude"]] = epa_margin_emis.loc[
                    (epa_margin_emis.CAMD_PLANT_ID == temp_plant)
                    & (epa_margin_emis.CAMD_GENERATOR_ID == temp_unit),
                    ["lat", "lon"],
                ].values.mean(0)
                print(
                    j,
                    temp_plant,
                    temp_unit,
                    dt.loc[
                        j,
                        [
                            "ave_marg_damages_isrm_Krewski",
                            "ave_marg_damages_isrm_LePeule",
                            "Latitude",
                            "Longitude",
                        ],
                    ].values,
                )

            elif (
                len(
                    epa_margin_emis.loc[
                        (epa_margin_emis.CAMD_PLANT_ID == temp_plant)
                        & (epa_margin_emis.CAMD_GENERATOR_ID == temp_unit)
                    ]
                )
                == 1
            ):
                dt.loc[
                    j,
                    ["ave_marg_damages_isrm_Krewski", "ave_marg_damages_isrm_LePeule"],
                ] = epa_margin_emis.loc[
                    (epa_margin_emis.CAMD_PLANT_ID == temp_plant)
                    & (epa_margin_emis.CAMD_GENERATOR_ID == temp_unit),
                    ["ave_marg_damages_isrm_Krewski", "ave_marg_damages_isrm_LePeule"],
                ].values[
                    0
                ]
                dt.loc[j, ["Latitude", "Longitude"]] = epa_margin_emis.loc[
                    (epa_margin_emis.CAMD_PLANT_ID == temp_plant)
                    & (epa_margin_emis.CAMD_GENERATOR_ID == temp_unit),
                    ["lat", "lon"],
                ].values[0]

        # plot marginal loc emissions
        if syst_name == "isone":
            x_lim = (-80, -66)
            y_lim = (38, 48)
        elif syst_name == "pjm":
            x_lim = (-93, -73)
            y_lim = (36, 43)
        else:
            x_lim = (-125, -66)
            y_lim = (25, 50)

        plot_marg_loc_emissions(
            dt=dt,
            ave_marg_damages_col="ave_marg_damages_isrm_LePeule",
            fname="ave_marg_damages_" + syst_name + "_map.png",
            x_lim=x_lim,
            y_lim=y_lim,
            markerscaler=20.0,
        )

        # impute missing values with Zone/Type/Fuel Averages
        dt.loc[(dt.Type == "Nuclear"), "ave_marg_damages_isrm_Krewski"] = 0.0
        dt.loc[(dt.Type == "Nuclear"), "ave_marg_damages_isrm_LePeule"] = 0.0
        dt["ave_marg_damages_imputed"] = False
        dt.loc[
            dt["ave_marg_damages_isrm_Krewski"].isna(), "ave_marg_damages_imputed"
        ] = True
        for t, z, f in (
            dt.loc[
                dt["ave_marg_damages_isrm_Krewski"].isna(),
                ["Type", "Zone", "Modeled Fuels"],
            ]
            .drop_duplicates()
            .values
        ):
            dt.loc[
                (dt["ave_marg_damages_isrm_Krewski"].isna())
                & (dt["Type"] == t)
                & (dt["Zone"] == z)
                & (dt["Modeled Fuels"] == f),
                ["ave_marg_damages_isrm_Krewski", "ave_marg_damages_isrm_LePeule"],
            ] = (
                dt.loc[
                    (~dt["ave_marg_damages_isrm_Krewski"].isna())
                    & (dt["Type"] == t)
                    & (dt["Zone"] == z)
                    & (dt["Modeled Fuels"] == f),
                    ["ave_marg_damages_isrm_Krewski", "ave_marg_damages_isrm_LePeule"],
                ]
                .mean()
                .values
            )

        # impute missing values with Type/Fuel Averages
        for t, f in (
            dt.loc[
                dt["ave_marg_damages_isrm_Krewski"].isna(), ["Type", "Modeled Fuels"]
            ]
            .drop_duplicates()
            .values
        ):
            dt.loc[
                (dt["ave_marg_damages_isrm_Krewski"].isna())
                & (dt["Type"] == t)
                & (dt["Modeled Fuels"] == f),
                ["ave_marg_damages_isrm_Krewski", "ave_marg_damages_isrm_LePeule"],
            ] = (
                dt.loc[
                    (~dt["ave_marg_damages_isrm_Krewski"].isna())
                    & (dt["Type"] == t)
                    & (dt["Modeled Fuels"] == f),
                    ["ave_marg_damages_isrm_Krewski", "ave_marg_damages_isrm_LePeule"],
                ]
                .mean()
                .values
            )

        # impute missing values with Type Averages
        for t in (
            dt.loc[dt["ave_marg_damages_isrm_Krewski"].isna(), "Type"]
            .drop_duplicates()
            .values
        ):
            dt.loc[
                (dt["ave_marg_damages_isrm_Krewski"].isna()) & (dt["Type"] == t),
                ["ave_marg_damages_isrm_Krewski", "ave_marg_damages_isrm_LePeule"],
            ] = (
                dt.loc[
                    (~dt["ave_marg_damages_isrm_Krewski"].isna()) & (dt["Type"] == t),
                    ["ave_marg_damages_isrm_Krewski", "ave_marg_damages_isrm_LePeule"],
                ]
                .mean()
                .values
            )

        # save file
        dt.to_csv(
            Path(
                Path.cwd(), "inputs", topo_name, "generators_loc_emission_damages.csv"
            ),
            index=False,
        )

        # delete temp file
        try:
            shutil.rmtree(_tmpdir)
        except Exception as e:
            print(e)


if __name__ == "__main__":
    # calling main function
    main()

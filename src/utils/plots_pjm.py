#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  13 11:46:47 2023

@author: Christoph Graf
"""

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiPolygon, shape
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import matplotlib.patches as patches
from pathlib import Path
import os
import argparse

# CHANGE THIS ACCORDING TO THE DATABASE
GSw_PJM = True

# result_dir = '/Users/skhanal/Globus_local/MOTEP-OSW/outputs/paper_draft/MO'
# COMMENT OUT LINES BELOW AND UNCOMMENT LINE ABOVE IF MANNUALLY RUNNING THIS SCRIPT
parser = argparse.ArgumentParser()
parser.add_argument('result_dir', type=str, help='The result directory')
args = parser.parse_args()
result_dir = args.result_dir

def plot_marg_loc_emissions(dt, ave_marg_damages_col):
    # us_counties_path = "https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_county_500k.zip"
    # us_counties = gpd.read_file(us_counties_path)

    us_states_path = (
        "https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_state_500k.zip"
    )
    us_states = gpd.read_file(us_states_path)

    # dt = pd.read_csv(Path(Path.cwd(), 'ave_marg_damages.csv'))
    dt = dt[~dt[ave_marg_damages_col].isna()]
    # dt["geometry"] = gpd.GeoSeries.from_wkt(dt["geometry"])
    # file:///C:/Users/c/Downloads/plant_power_eia_v6.pdf
    # dt = gpd.GeoDataFrame(dt, geometry="geometry", crs="EPSG:3857")
    # dt = dt.to_crs('EPSG:4269') # us_states.crs
    dt["markersize"] = 10 * dt[ave_marg_damages_col] / max(dt[ave_marg_damages_col])

    ax = us_states.boundary.plot(edgecolor="black", linewidth=0.2, figsize=(10, 5))
    # add ave_marg_damages for each power plant
    dt.plot(ax=ax, marker="o", color="red", markersize=dt.markersize)

    plt.tight_layout()
    plt.xlim(-80, -66)
    plt.ylim(38, 48)
    plt.axis("off")
    plt.savefig(
        Path(Path.cwd(), "figures", "ave_marg_damages_map.png"),
        dpi=800,
        format="png",
        metadata=None,
        bbox_inches=None,
        pad_inches=0.1,
        facecolor="auto",
        edgecolor="auto",
        backend=None,
    )


def make_spatial_transmission_plot(
    dt_zones,
    dt_lines,
    bb_adder,
    bb_adder_focus,
    trans_line_scaler,
    fig_name,
):
    """
    inputs:
        dt_zones is a pd.Dataframe with the following cols:
            name (str); geometry (gpd)

        dt_lines is a pd.Dataframe with the following cols:
            name (str); from (str); to (str); mw (float), exisiting (bool), line_type (str)

    output:
        spatial figure with transmission lines

    """
    dt_zones = dt_zones.copy()
    dt_lines = dt_lines.copy()

    # sort data
    dt_zones = dt_zones.sort_values(["focus", "name"]).reset_index(drop=True)
    dt_lines = dt_lines.sort_values("line_type").reset_index(drop=True)

    # add centroid to zonal data
    dt_zones["centroid"] = dt_zones.geometry.centroid
    dt_zones["centroid_long"] = [
        dt_zones["centroid"].values[j].coords[0][0] for j in range(len(dt_zones))
    ]
    dt_zones["centroid_lat"] = [
        dt_zones["centroid"].values[j].coords[0][1] for j in range(len(dt_zones))
    ]

    # zone color
    # TODO: Breaks if there are more than 8 zones
    dt_zones["zone_color"] = [
        mpl.cm.Pastel2(j) if j < 8 else mpl.cm.Pastel1(j - 8)
        for j in dt_zones["name"].factorize()[0]
    ]

    # dt_zones["zone_color"] = np.nan

    # different color scheme for focus area (offshore zones)
    dt_zones.loc[dt_zones.focus == True, "zone_color"] = pd.Series(
        [
            mpl.cm.Dark2(j)
            for j in dt_zones.loc[dt_zones.focus == True, "name"].factorize()[0]
        ],
        index=dt_zones[dt_zones.focus == True].index,
    )

    # min/max long/lat
    x_lim = (
        dt_zones.bounds["minx"].min() - bb_adder,
        dt_zones.bounds["maxx"].max() + bb_adder,
    )
    y_lim = (
        dt_zones.bounds["miny"].min() - bb_adder,
        dt_zones.bounds["maxy"].max() + bb_adder,
    )

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax1 = dt_zones.boundary.plot(edgecolor="black", linewidth=0.2, ax=ax)
    # ax = dt_zones.boundary.plot(edgecolor="black", linewidth=0.2, figsize=(10, 5))

    """
    for k, v in enumerate(dt_zones.loc[dt_zones.focus == False].index):
        gpd.GeoSeries(dt_zones.loc[v, "geometry"]).plot(
                ax=ax,
                color=mpl.cm.Pastel2(k),
                alpha=0.7,
                label=dt_zones.loc[v, "name"],
                edgecolor="k",
            )
        

    ax.legend(loc="lower right", fontsize=12, frameon=False)

    """
    ax2 = dt_zones.loc[dt_zones.focus == False].plot(
        column=dt_zones.loc[dt_zones.focus == False, "name"],
        # color=dt_zones.loc[dt_zones.focus == False, "zone_color"],
        cmap=colors.ListedColormap(
            list(dt_zones.loc[dt_zones.focus == False, "zone_color"].values)
        ),
        # cmap="Pastel2",
        alpha=0.7,
        edgecolor="k",
        legend=True,
        legend_kwds={
            "loc": "center left",
            "bbox_to_anchor": (1, 0.5),
            "ncol": 2,
            #    "labels": list(dt_zones.loc[dt_zones.focus == False, "name"].values),
        },
        ax=ax,
    )

    """
    legend1 = ax.legend(
        scatterpoints=1,
        loc="upper right",
        bbox_to_anchor=(1.0, 0.5),
    )

    fig.gca().add_artist(legend1)
    """

    leg = ax.get_legend()
    ax.add_artist(leg)

    # add zone colors
    for j in range(len(leg.get_lines())):
        temp_row = dt_zones.loc[dt_zones.focus == False].index[j]
        dt_zones._set_value(
            temp_row, "zone_color", leg.get_lines()[j].get_markerfacecolor()
        )

    # # handles, labels = ax.get_legend_handles_labels()
    # print(leg, handles, labels)

    # leg = plt.legend(loc=(1.03, 0), title="Year")
    # ax.add_artist(leg)

    # ax.legend()
    # handles, labels = ax.get_legend_handles_labels()
    # handles, labels = plt.gca().get_legend_handles_labels()
    # print(handles, labels)
    # handles, labels = plt.gca().get_legend_handles_labels()

    # add lines
    dt_lines = dt_lines.merge(
        dt_zones[["name", "centroid_long", "centroid_lat"]].rename(
            {
                "name": "from",
                "centroid_long": "from_centroid_long",
                "centroid_lat": "from_centroid_lat",
            },
            axis=1,
        ),
        how="left",
        on="from",
    )
    dt_lines = dt_lines.merge(
        dt_zones[["name", "centroid_long", "centroid_lat"]].rename(
            {
                "name": "to",
                "centroid_long": "to_centroid_long",
                "centroid_lat": "to_centroid_lat",
            },
            axis=1,
        ),
        how="left",
        on="to",
    )
    dt_lines = dt_lines.reset_index(drop=True)
    dt_lines["linewidth"] = trans_line_scaler * (dt_lines["mw"] / max(dt_lines["mw"]))
    custom_colors = ["yellow", "green", "red", "blue", "orange", "purple"]  # Define your colors here
    dt_lines["line_type_color"] = [
        custom_colors[j % len(custom_colors)] for j in dt_lines["line_type"].factorize()[0]
    ]
    # dt_lines["line_type_color"] = [
    #     mpl.cm.Set1(1 + j) for j in dt_lines["line_type"].factorize()[0]
    # ]

    for j in dt_lines.index:
        from_temp = (
            dt_lines.loc[j, "from_centroid_long"],
            dt_lines.loc[j, "to_centroid_long"],
        )
        to_temp = (
            dt_lines.loc[j, "from_centroid_lat"],
            dt_lines.loc[j, "to_centroid_lat"],
        )

        ax.plot(
            from_temp,
            to_temp,
            # color="black",
            color=dt_lines.loc[j, "line_type_color"],
            linewidth=dt_lines.loc[j, "linewidth"],
            alpha=0.7,
            # label=dt_lines.loc[j, "line_type"],
        )

    # adapt legend
    cust_legend_elements = [
        Line2D([0], [0], lw=trans_line_scaler / 2, color=col, label=typ)
        for typ, col in dt_lines[["line_type", "line_type_color"]]
        .drop_duplicates()
        .values
    ]
    ax.legend(
        handles=cust_legend_elements,
        # loc="upper left",
        bbox_to_anchor=(1.05, 1),
        loc=2,
        borderaxespad=0.0,
    )

    # plot centroid of each zone
    dt_zones["centroid"].plot(ax=ax, marker="o", color="red", markersize=2, zorder=2)

    # cust_legend = plt.legend(cust_legend_elements, loc=4)
    # ax.add_artist(cust_legend_elements)
    # plt.legend(handles=h, labels=range(5, 13), loc=(1.03, 0.5), title="Quality")
    # ax.legend(handles=cust_legend_elements, loc="center")

    # add bounding box over focus area
    if np.any(dt_zones.focus):
        x_lim_focus = (
            dt_zones.loc[dt_zones.focus == True].bounds["minx"].min()
            - bb_adder_focus / 2.0,
            dt_zones.loc[dt_zones.focus == True].bounds["maxx"].max()
            + bb_adder_focus / 2.0,
        )

        y_lim_focus = (
            dt_zones.loc[dt_zones.focus == True].bounds["miny"].min()
            - bb_adder_focus / 2.0,
            dt_zones.loc[dt_zones.focus == True].bounds["maxy"].max()
            + bb_adder_focus / 2.0,
        )

        rect = patches.Rectangle(
            (x_lim_focus[0], y_lim_focus[0]),
            x_lim_focus[1] - x_lim_focus[0],
            y_lim_focus[1] - y_lim_focus[0],
            linewidth=0.2,
            linestyle="dashed",
            edgecolor="0.2",
            facecolor="none",
        )

        # Add the patch to the Axes
        ax.add_patch(rect)

    # save figure
    plt.tight_layout()
    plt.xlim(x_lim)
    plt.ylim(y_lim)
    plt.axis("off")
    plt.savefig(
        Path(Path.cwd(), "figures", fig_name + ".png"),
        dpi=800,
        format="png",
        metadata=None,
        bbox_inches=None,
        pad_inches=0.1,
        facecolor="auto",
        edgecolor="auto",
        backend=None,
    )

    plt.close()

    if np.any(dt_zones.focus):
        dt_zones_focus = dt_zones.copy()
        dt_lines_focus = dt_lines.copy()
        focus_box_poly = Polygon(
            [
                (x_lim_focus[0], y_lim_focus[0]),
                (x_lim_focus[1], y_lim_focus[0]),
                (x_lim_focus[1], y_lim_focus[1]),
                (x_lim_focus[0], y_lim_focus[1]),
            ]
        )
        for k in dt_zones_focus.index:
            if focus_box_poly.disjoint(dt_zones_focus.loc[k, "geometry"]):
                dt_zones_focus = dt_zones_focus.drop(k, axis=0)

        dt_zones_focus = dt_zones_focus.sort_values("name").reset_index(drop=True)

        # dt_lines_focus = dt_lines_focus.loc[dt_lines_focus["from"].isin(dt_zones_focus["name"]) | dt_lines_focus["to"].isin(dt_zones_focus["name"])]
        # dt_lines_focus = dt_lines_focus.sort_values("line_type").reset_index(drop=True)

        # plot
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        ax1 = dt_zones_focus.boundary.plot(edgecolor="black", linewidth=0.2, ax=ax)
        ax2 = dt_zones_focus.plot(
            column=dt_zones_focus["name"].values,
            cmap=colors.ListedColormap(list(dt_zones_focus["zone_color"].values)),
            # cmap="Pastel2",
            alpha=0.7,
            edgecolor="k",
            legend=True,
            legend_kwds={"loc": "center left", "bbox_to_anchor": (1, 0.5)},
            ax=ax,
        )

        # plot centroid of each zone
        dt_zones_focus.loc[dt_zones_focus.focus == False, "centroid"].plot(
            ax=ax,
            marker="o",
            color="red",
            markersize=2,
            zorder=3,
        )

        # redraw offshore zones
        point_geo = np.array([j.type == "Point" for j in dt_zones_focus.geometry])
        dt_zones_focus.loc[point_geo].plot(
            column=dt_zones_focus.loc[point_geo, "name"].values,
            cmap=colors.ListedColormap(
                list(dt_zones_focus.loc[point_geo, "zone_color"].values)
            ),
            alpha=1.0,
            edgecolor="k",
            legend=False,
            ax=ax,
            zorder=4,
        )

        leg = ax.get_legend()
        for k, v in enumerate(leg.legendHandles):
            if point_geo[k]:
                leg.legendHandles[k] = v.set_alpha(1.0)
        ax.add_artist(leg)

        for j in dt_lines_focus.index:
            from_temp = (
                dt_lines.loc[j, "from_centroid_long"],
                dt_lines.loc[j, "to_centroid_long"],
            )
            to_temp = (
                dt_lines.loc[j, "from_centroid_lat"],
                dt_lines.loc[j, "to_centroid_lat"],
            )

            ax.plot(
                from_temp,
                to_temp,
                # color="black",
                color=dt_lines.loc[j, "line_type_color"],
                linewidth=dt_lines.loc[j, "linewidth"],
                alpha=0.7,
                # label=dt_lines.loc[j, "line_type"],
            )

        # adapt legend
        cust_legend_elements = [
            Line2D([0], [0], lw=trans_line_scaler / 2, color=col, label=typ)
            for typ, col in dt_lines[["line_type", "line_type_color"]]
            .drop_duplicates()
            .values
        ]

        ax.legend(handles=cust_legend_elements, loc="lower left")

        plt.tight_layout()
        plt.xlim(x_lim_focus)
        plt.ylim(y_lim_focus)
        plt.axis("off")

        plt.savefig(
            Path(Path.cwd(), "figures", fig_name + "_focus.png"),
            dpi=800,
            format="png",
            metadata=None,
            bbox_inches=None,
            pad_inches=0.1,
            facecolor="auto",
            edgecolor="auto",
            backend=None,
        )
        plt.close()

    """
    if np.any(dt_zones.focus):
        # save figure zooming into focus area
        plt.tight_layout()
        plt.xlim(x_lim_focus)
        plt.ylim(y_lim_focus)
        plt.axis("off")
        # add wind power labels
        for j in dt_zones[dt_zones.focus == True].index:
            temp_name = dt_zones.loc[j, "name"]
            temp_x = dt_zones.loc[j, "geometry"].x
            temp_y = dt_zones.loc[j, "geometry"].y
            plt.text(temp_x, temp_y, temp_name, bbox={"pad": 0.2})

        plt.savefig(
            Path(Path.cwd(), "figures", fig_name + "_focus.png"),
            dpi=800,
            format="png",
            metadata=None,
            bbox_inches=None,
            pad_inches=0.1,
            facecolor="auto",
            edgecolor="auto",
            backend=None,
        )
    """


def main():
    us_states_path = (
        "https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_state_500k.zip"
    )
    us_states = gpd.read_file(us_states_path)

    us_counties_path = (
        "https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_county_500k.zip"
    )
    us_counties = gpd.read_file(us_counties_path)
    states2drop = us_states.loc[
        us_states["NAME"].isin(
            [
                "Alaska",
                "Hawaii",
                "Puerto Rico",
                "Commonwealth of the Northern Mariana Islands",
                "United States Virgin Islands",
                "American Samoa",
                "Guam",
            ]
        ),
        "STATEFP",
    ].values
    us_counties = us_counties.dropna()
    us_counties = us_counties[~us_counties["STATEFP"].isin(states2drop)]
    us_counties = us_counties.merge(
        us_states[["NAME", "STUSPS", "STATEFP", "geometry"]].rename(
            {"NAME": "STATE", "geometry": "geometry_zone"}, axis=1
        ),
        on="STATEFP",
        how="left",
    )

    isone_states = [
        "Connecticut",
        "Maine",
        "Massachusetts",
        "New Hampshire",
        "Rhode Island",
        "Vermont",
    ]

    dt_zones = us_counties.loc[us_counties.STATE.isin(isone_states)]
    dt_zones["zone"] = pd.NA

    # FIXME: IS THIS MAPPING CORRECT?
    # dt_zones.loc[dt_zones.zone==z, 'geometry_zone'].unary_union
    dt_zones.loc[
        (
            dt_zones.NAME.isin(
                ["Berkshire", "Franklin", "Hampshire", "Hampden", "Worcester"]
            )
        )
        & (dt_zones.STATE == "Massachusetts"),
        "zone",
    ] = "WCMA"
    dt_zones.loc[
        (
            dt_zones.NAME.isin(
                ["Plymouth", "Bristol", "Barnstable", "Dukes", "Nantucket"]
            )
        )
        & (dt_zones.STATE == "Massachusetts"),
        "zone",
    ] = "SEMA"
    dt_zones.loc[
        (dt_zones.zone.isna()) & (dt_zones.STATE == "Massachusetts"), "zone"
    ] = "NEMA"

    dt_zones.loc[dt_zones.zone.isna(), "zone"] = dt_zones.loc[
        dt_zones.zone.isna(), "STUSPS"
    ]

    dt_zones.loc[
        dt_zones["STUSPS"] != dt_zones["zone"], "geometry_zone"
    ] = dt_zones.loc[dt_zones["STUSPS"] != dt_zones["zone"], "geometry"]
    dt_zones = dt_zones[["zone", "geometry_zone"]].drop_duplicates()
    dt_zones = dt_zones[["zone", "geometry_zone"]].rename(
        {"zone": "name", "geometry_zone": "geometry"}, axis=1
    )
    dt_zones = dt_zones.set_geometry("geometry", crs=us_counties.crs)
    # merge counties to zones
    for z in dt_zones.name.unique():
        if sum(dt_zones.name == z) > 1:
            temp_index = dt_zones.loc[dt_zones.name == z].index
            temp_geo_union = dt_zones.loc[temp_index, "geometry"].unary_union

            # drop existing items
            dt_zones = dt_zones.drop(temp_index, axis=0)
            # replace by merged zone
            dt_zones = pd.concat(
                [
                    dt_zones,
                    gpd.GeoDataFrame(
                        {"name": [z], "geometry": [shape(temp_geo_union)]},
                        crs=dt_zones.crs,
                    ),
                ]
            )

    dt_zones = dt_zones.reset_index(drop=True)
    dt_zones["focus"] = False

    # add offshore zones
    offshore = pd.DataFrame(
        {
            "name": [
                "REV",
                "VINE",
                "PKCTY",
                "COMW",
                "MFLR1",
                "MFLR2",
            ],  # "Revolution Wind", "Vineyard Wind 1", "Park City Wind", "Commonwealth Wind", "Mayflower Wind 1", "Mayflower Wind 2"
            "geometry": [
                Point(-71.069972, 41.149972),
                Point(-70.61667, 41.03325),
                Point(-70.709972, 40.859972),
                Point(-70.695545, 40.77634),
                Point(-70.616667, 40.797153),
                Point(-70.547001, 40.806968),
            ],
        }
    )

    offshore["focus"] = True

    dt_zones = pd.concat([dt_zones, offshore], axis=0).reset_index(drop=True)

    dt_lines = pd.DataFrame(
        {
            "name": [
                "CT_REV",
                "RI_REV",
                "SEMA_VINE",
                "SEMA_PKCTY",
                "SEMA_COMW",
                "SEMA_MFLR1",
                "SEMA_MFLR2",
                "ME_NH",
                "NH_VT",
                "NH_WCMA",
                "NH_NEMA",
                "WCMA_NEMA",
                "WCMA_CT",
                "WCMA_RI",
                "RI_SEMA",
                "SEMA_NEMA",
                "CT_RI",
            ],
            "from": [
                "CT",
                "RI",
                "SEMA",
                "SEMA",
                "SEMA",
                "SEMA",
                "SEMA",
                "ME",
                "NH",
                "NH",
                "NH",
                "WCMA",
                "WCMA",
                "WCMA",
                "RI",
                "SEMA",
                "CT",
            ],
            "to": [
                "REV",
                "REV",
                "VINE",
                "PKCTY",
                "COMW",
                "MFLR1",
                "MFLR2",
                "NH",
                "VT",
                "WCMA",
                "NEMA",
                "NEMA",
                "CT",
                "RI",
                "SEMA",
                "NEMA",
                "RI",
            ],
            "mw": [300.0, 400.0, 800.0, 800.0, 1200.0, 800.0, 400] + [1000.0] * 10,
            "existing": [False, False, False, False, False, False, False] + [True] * 10,
            "line_type": ["AC", "DC", "AC", "AC", "DC", "AC", "AC"] + ["AC"] * 10,
        }
    )

    # read in external lines file
    dt_lines = pd.read_csv(Path(result_dir, "line_existing_upgrades_new.csv"))
    dt_lines = dt_lines[["From Zone", "To Zone", "s_max", "Line Type"]].rename(
        {
            "From Zone": "from",
            "To Zone": "to",
            "s_max": "mw",
            "Line Type": "line_type",
        },
        axis=1,
    )
    dt_lines.loc[dt_lines["from"] == "NEMA/BOST", "from"] = "NEMA"
    dt_lines.loc[dt_lines["to"] == "NEMA/BOST", "to"] = "NEMA"
    dt_lines["name"] = dt_lines["from"] + "_" + dt_lines["to"]

    make_spatial_transmission_plot(
        dt_zones,
        dt_lines,
        bb_adder=2.0,
        bb_adder_focus=1.5,
        trans_line_scaler=5.0,
        fig_name=os.path.basename(result_dir),
    )

    if GSw_PJM:
        # make plots for PJM
        epa_ipm_shape_path = Path(Path.cwd(), "inputs", "ipm_v6_regions.zip")
        epa = gpd.read_file(epa_ipm_shape_path, crs="EPSG:4326")
        epa = epa.to_crs("EPSG:4269")
        dt_zones_pjm = epa.loc[epa.IPM_Region.str.startswith("PJM_")].rename(
            {"IPM_Region": "name"}, axis=1
        )
        dt_zones_pjm["focus"] = False
        # add offshore wind locations
        dt_zones_pjm = pd.concat(
            [
                dt_zones_pjm,
                pd.DataFrame(
                    {
                        "name": [
                            "ATLSHR",
                            "OCN1",
                            "OCN2",
                            "SKPJK1",
                            "SKPJK2",
                            "MRWN",
                            "MMTM",
                            "CVOW",
                        ],
                        "geometry": [
                            Point(-74.200000, 39.189972),
                            Point(-74.299972, 39.100000),
                            Point(-74.424, 39.071),
                            Point(-74.700000, 38.649972),
                            Point(-74.700000, 38.670000),
                            Point(-74.778, 38.252),
                            Point(-74.760000, 38.340000),
                            Point(-75.217, 36.947),
                        ],
                        "focus": [True] * 8,
                    }
                ),
            ]
        ).reset_index(drop=False)

        dt_lines_pjm = dt_lines = pd.read_csv(Path(result_dir, "line_existing_upgrades_new.csv"))
        dt_lines_pjm["name"] = (
            dt_lines_pjm["From Zone"].str.replace("PJM_", "")
            + "_"
            + dt_lines_pjm["To Zone"].str.replace("PJM_", "")
        )
        dt_lines_pjm["existing"] = True
        dt_lines_pjm = dt_lines_pjm.rename(
            {"From Zone": "from", "To Zone": "to", "s_max": "mw", "Line Type": "line_type"},
            axis=1,
        )
        dt_lines_pjm = dt_lines_pjm[["name", "from", "to", "mw", "existing", "line_type"]]

        make_spatial_transmission_plot(
            dt_zones_pjm,
            dt_lines_pjm,
            bb_adder=2.5,
            bb_adder_focus=3.5,
            trans_line_scaler=15.0,
            fig_name=os.path.basename(result_dir),
        )

if __name__ == "__main__":
    # calling main function
    main()

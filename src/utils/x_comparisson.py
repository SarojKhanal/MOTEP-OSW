#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  13 11:46:47 2023

@author: Christoph Graf
"""

import numpy as np
import pandas as pd
import itertools
from pathlib import Path


def x_compare_annual_results(results_path_list, attr_list, aggr=False):
    # results_paths is a list of paths pointing to folder where results are stored and that will be compared to each other
    # attr is a list of attributes that will be compared
    # aggr (binary); true then annual values are summed up)
    dt = pd.concat(
        [
            pd.read_csv(Path(j1, j2), header=None)
            .assign()
            .assign(Case=j1.name)
            .assign(Attribute=j2.replace(".csv", ""))
            for j1, j2 in itertools.product(results_path_list, attr_list)
        ]
    )
    dt = dt.rename({0: "Value"}, axis=1)
    dt = dt.reset_index().rename({"index": "Year"}, axis=1)
    if aggr:
        dt = dt.groupby(["Case", "Attribute"], as_index=False)["Value"].sum()

    return dt


def main():
    x_compare_annual_results(
        results_path_list=[
            Path(
                Path.cwd().parent,
                "20230622_Biweekly_meeting",
                "paper_draft",
                "Fixed POI",
            ),
            Path(
                Path.cwd().parent,
                "20230622_Biweekly_meeting",
                "paper_draft",
                "Fixed POI UC",
            ),
        ],
        attr_list=["annual_cost_operation.csv"],
        aggr=True,
    )

    x_compare_annual_results(
        results_path_list=[
            Path(
                Path.cwd().parent,
                "20230622_Biweekly_meeting",
                "paper_draft",
                "Fixed POI CostOnlyObj",
            ),
            Path(
                Path.cwd().parent,
                "20230622_Biweekly_meeting",
                "paper_draft",
                "Fixed POI",
            ),
        ],
        attr_list=["annual_cost_operation.csv"],
        aggr=True,
    )

    x_compare_annual_results(
        results_path_list=[
            Path(
                Path.cwd().parent,
                "20230622_Biweekly_meeting",
                "paper_draft",
                "Fixed POI CostOnlyObj UC",
            ),
            Path(
                Path.cwd().parent,
                "20230622_Biweekly_meeting",
                "paper_draft",
                "Fixed POI UC",
            ),
        ],
        attr_list=[
            "annual_cost_operation.csv",
            "annual_cost_gen_investment.csv",
            "annual_cost_bess_investment.csv",
            "annual_cost_line_investment.csv",
            "annual_cost_externalities.csv",
            "annual_cost_operation_under_generation.csv",
        ],
        aggr=True,
    )

    x_compare_annual_results(
        results_path_list=[
            Path(
                Path.cwd().parent,
                "20230622_Biweekly_meeting",
                "paper_draft",
                "Fixed POI CostOnlyObj",
            ),
            Path(
                Path.cwd().parent,
                "20230622_Biweekly_meeting",
                "paper_draft",
                "Fixed POI",
            ),
        ],
        attr_list=[
            "annual_cost_operation.csv",
            "annual_cost_gen_investment.csv",
            "annual_cost_bess_investment.csv",
            "annual_cost_line_investment.csv",
            "annual_cost_externalities.csv",
            "annual_cost_operation_under_generation.csv",
        ],
        aggr=True,
    )

    x_compare_annual_results(
        results_path_list=[
            Path(
                Path.cwd().parent,
                "20230622_Biweekly_meeting",
                "paper_draft",
                "Fixed POI",
            ),
            Path(
                Path.cwd().parent,
                "20230622_Biweekly_meeting",
                "paper_draft",
                "Fixed POI UC",
            ),
        ],
        attr_list=[
            "new_build_gen_cap_bess_power.csv",
            "new_build_gen_cap_ng_cc_ccs.csv",
            "new_build_gen_cap_ng_ct.csv",
            "new_build_gen_cap_pv.csv",
            "new_build_gen_cap_wind.csv",
            "new_line_decision.csv",
        ],
        aggr=True,
    )

    x_compare_annual_results(
        results_path_list=[
            # Path(Path.cwd().parent, "20230417_Results", "_RESULTS", "Base SpillCost0"),
            Path(
                Path.cwd().parent,
                "20230622_Biweekly_meeting",
                "paper_draft",
                "Fixed POI",
            ),
            Path(
                Path.cwd().parent,
                "20230622_Biweekly_meeting",
                "paper_draft",
                "Fixed POI UC",
            ),
            # Path(Path.cwd().parent, "20230417_Results", "_RESULTS", "Fixed POI FullUC SpillCost0"),
        ],
        attr_list=[
            "annual_cost_operation.csv",
            "annual_cost_bess_investment.csv",
            "annual_cost_gen_investment.csv",
            "annual_cost_line_investment.csv",
            "annual_cost_externalities.csv",
            "annual_cost_externalities_airquality.csv",
            "annual_cost_externalities_emission.csv",
            "annual_cost_operation_demand_flexibility.csv",
            "annual_cost_total.csv",
            # "new_build_gen_cap_bess_energy.csv",
            # "annual_cost_externalities_airquality.csv",
            # "annual_cost_externalities_emission.csv"
        ],
        aggr=True,
    )

    x_compare_annual_results(
        results_path_list=[
            # Path(Path.cwd().parent, "20230417_Results", "_RESULTS", "Base SpillCost0"),
            # Path(Path.cwd().parent, "20230417_Results", "_RESULTS", "Base FullUC SpillCost0"),
            Path(
                Path.cwd().parent,
                "20230417_Results",
                "_RESULTS",
                "Fixed POI SpillCost0",
            ),
            Path(
                Path.cwd().parent,
                "20230417_Results",
                "_RESULTS X",
                "Fixed POI SpillCost0 X3",
            ),
        ],
        attr_list=[
            "annual_cost_operation.csv",
            "annual_cost_bess_investment.csv",
            "annual_cost_gen_investment.csv",
            "annual_cost_line_investment.csv",
            "annual_cost_externalities.csv",
            "annual_cost_externalities_airquality.csv",
            "annual_cost_externalities_emission.csv",
            "annual_cost_operation_demand_flexibility.csv",
            "annual_cost_total.csv",
            # "new_build_gen_cap_bess_energy.csv",
            # "annual_cost_externalities_airquality.csv",
            # "annual_cost_externalities_emission.csv"
        ],
        aggr=True,
    )

    x_compare_annual_results(
        results_path_list=[
            # Path(Path.cwd().parent, "20230417_Results", "_RESULTS", "Base SpillCost0"),
            Path(
                Path.cwd().parent,
                "20230417_Results",
                "_RESULTS",
                "Base FullUC SpillCost0",
            ),
            # Path(Path.cwd().parent, "20230417_Results", "_RESULTS", "Optimized POI SpillCost0"),
            Path(
                Path.cwd().parent,
                "20230417_Results",
                "_RESULTS",
                "Optimized POI FullUC SpillCost_On0_Off0",
            ),
        ],
        attr_list=[
            "annual_cost_externalities_airquality.csv",
            "annual_cost_externalities_emission.csv",
        ],
        aggr=True,
    )

    x_compare_annual_results(
        results_path_list=[
            Path(
                Path.cwd().parent,
                "20230417_Results",
                "_RESULTS",
                "Base FullUC SpillCost0",
            ),
            Path(
                Path.cwd().parent,
                "20230417_Results",
                "_RESULTS",
                "Fixed POI FullUC SpillCost0",
            ),
            Path(
                Path.cwd().parent,
                "20230417_Results",
                "_RESULTS",
                "Optimized POI FullUC SpillCost_On0_Off0",
            ),
        ],
        attr_list=["annual_cost_operation.csv"],
        aggr=True,
    )

    x_compare_annual_results(
        results_path_list=[
            Path(
                Path.cwd().parent,
                "20230417_Results",
                "_RESULTS",
                "Base FullUC SpillCost0",
            ),
            Path(
                Path.cwd().parent,
                "20230417_Results",
                "_RESULTS",
                "Fixed POI FullUC SpillCost0",
            ),
            Path(
                Path.cwd().parent,
                "20230417_Results",
                "_RESULTS",
                "Fixed POI FullUC SpillCost0 FlexRPS",
            ),
            Path(
                Path.cwd().parent,
                "20230417_Results",
                "_RESULTS",
                "Fixed POI FullUC SpillCost0 DemFlex50",
            ),
        ],
        attr_list=["annual_cost_operation.csv"],
        aggr=True,
    )


if __name__ == "__main__":
    # calling main function
    main()

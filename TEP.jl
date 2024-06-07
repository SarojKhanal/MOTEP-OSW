# Set up environment from .toml files
# ]add DataFrames, CSV, JuMP, Gurobi, Plots, IJulia, Conda
import Pkg
Pkg.activate(".")
Pkg.instantiate()

# Load necessary packages
using LinearAlgebra, DataFrames, CSV, JuMP, Gurobi, Plots

# Load functions and model
include(joinpath(@__DIR__, "src", "input.jl")) # Type definitions and read-in functions
include(joinpath(@__DIR__, "src", "model.jl")) # Model definition # NOTE: change "model.jl" to "model_pjm.jl" if you're running for PJM test system
include(joinpath(@__DIR__, "src", "output.jl")) # Postprocessing of solved model # NOTE: change "output.jl" to "output_pjm.jl" if you're running for PJM test system

function diagnose_infeasibility(model)
    if termination_status(model) == MOI.INFEASIBLE
        println("Model is infeasible")

        # Access the Gurobi model
        grb_model = model.moi_backend.optimizer.model

        # Call Gurobi's computeIIS function
        Gurobi.GRBcomputeIIS(grb_model)

        # Write the IIS to a file
        Gurobi.GRBwrite(grb_model, "report_infeasibility.ilp")
        println("IIS written to report_infeasibility.ilp")

        # Export the model to MPS format
        # Gurobi.GRBwrite(grb_model, "model.mps")
        # println("Model exported to model.mps")
    end
end

cases = CSV.read(joinpath(@__DIR__, "cases.csv"), DataFrame, header=true, types=String)
for case in names(cases)[5:end]
    @info "################$(repeat("#", length(case)))################"
    @info "##### RUNNING CASE: \"$case\" ... #####"
    @info "################$(repeat("#", length(case)))################"
    # Scenario name for output folder
    global output_dir = isempty(case) ? joinpath(@__DIR__, "outputs","base") : joinpath(@__DIR__, "outputs",case)
    if !ispath(output_dir)
        mkpath(output_dir; mode = 0o777)
    end
    # Formulation switches and settings
    cases[!,case] .= coalesce.(cases[!,case], cases[!,"Default Value"])
    local cases_paramsopts = Dict(row[Symbol("Index")] => row[(case)] for row in eachrow(cases))
    # Write original "cases.csv" file the dictionary to a CSV file
    CSV.write(joinpath(output_dir, "case.csv"), cases) #!TODO: Only output corresponding case.csv, remove other columns
    CSV.write(joinpath(output_dir, "cases_paramsopts_" * case * ".csv"), pairs(cases_paramsopts), header=["Keys", "Values"])
    global datadir = abspath(joinpath(String(cases_paramsopts["datadir"])))
    global MIPGap = parse(Float64, String(cases_paramsopts["MIPGap"]))
    global ObjScale = parse(Float64, String(cases_paramsopts["ObjScale"]))
    global GSw_AnnualResolution = parse(Bool, String(cases_paramsopts["GSw_AnnualResolution"]))
    global GSw_RPS = parse(Bool, String(cases_paramsopts["GSw_RPS"]))
    global GSw_MultiObj = parse(Bool, String(cases_paramsopts["GSw_MultiObj"]))
    global GSw_Emission = parse(Bool, String(cases_paramsopts["GSw_Emission"]))
    global P_CO2 = parse(Float64, String(cases_paramsopts["P_CO2"]))
    global GSw_AirQuality = parse(Bool, String(cases_paramsopts["GSw_AirQuality"]))
    global k_scale_AQ = parse(Float64, String(cases_paramsopts["k_scale_AQ"]))
    global GSw_DemandFlexibility = parse(Bool, String(cases_paramsopts["GSw_DemandFlexibility"]))
    global αᶠˡᵉˣ = parse(Float64, String(cases_paramsopts["alpha_flex"]))
    global GSw_Battery = parse(Bool, String(cases_paramsopts["GSw_Battery"]))
    global GSw_BatteryOSW = parse(Bool, String(cases_paramsopts["GSw_BatteryOSW"]))
    global GSw_OSW = parse(Bool, String(cases_paramsopts["GSw_OSW"]))
    global n_lines_existing = parse(Int, String(cases_paramsopts["n_lines_existing"])) # !TODO: is this necessary? Automate this.
    global ηₛᶜʰ = parse(Float64, String(cases_paramsopts["eff_charging"]))
    global ηₛᵈⁱˢ = parse(Float64, String(cases_paramsopts["eff_discharging"]))
    global k_cal = parse(Float64, String(cases_paramsopts["k_cal"]))
    global k_cycle = parse(Float64, String(cases_paramsopts["k_cycle"]))
    global battery_size_hr = parse(Float64, String(cases_paramsopts["battery_size_hr"]))
    global GSw_WindPTC = parse(Bool, String(cases_paramsopts["GSw_WindPTC"]))
    global GSw_PVITC = parse(Bool, String(cases_paramsopts["GSw_PVITC"]))
    global GSw_ExogenousRetirements = parse(Bool, String(cases_paramsopts["GSw_ExogenousRetirements"]))
    global k_scale_OnshoreLineCost = parse(Float64, String(cases_paramsopts["k_scale_OnshoreLineCost"]))
    global GSw_RPSNC = parse(Bool, String(cases_paramsopts["GSw_RPSNC"]))
    global RPSNC = parse(Float64, String(cases_paramsopts["RPSNC"]))
    global GSw_FullUC = parse(Bool, String(cases_paramsopts["GSw_FullUC"]))
    global VoLL = parse(Float64, String(cases_paramsopts["VoLL"]))
    global VoWS = parse(Float64, String(cases_paramsopts["VoWS"]))
    global VoOSWS = parse(Float64, String(cases_paramsopts["VoOSWS"]))
    global r = parse(Float64, String(cases_paramsopts["r"]))
    global δt = parse(Float64, String(cases_paramsopts["delta_t"]))
    global LT = parse(Int, String(cases_paramsopts["LT"]))
    global GT = parse(Int, String(cases_paramsopts["GT"]))
    global ST = parse(Int, String(cases_paramsopts["ST"]))
    global n_buses_existing = parse(Int, String(cases_paramsopts["n_buses_existing"]))
    global Y = parse(Int, String(cases_paramsopts["Y"]))
    global Y_start = parse(Int, String(cases_paramsopts["Y_start"]))
    global e_num = parse(Int, String(cases_paramsopts["e_num"]))
    global GSw_ExtremeScenarios = parse(Bool, String(cases_paramsopts["GSw_ExtremeScenarios"]))
    global e_numx = parse(Int, String(cases_paramsopts["e_numx"]))
    global VoLLX = parse(Float64, String(cases_paramsopts["VoLLX"]))
    global αᴱ = parse(Float64, String(cases_paramsopts["alpha_E"]))
    global α⁺ = parse(Float64, String(cases_paramsopts["alpha_plus"]))

    # Load data
    local network = load_network(datadir)
    local wind_power = Matrix(CSV.read(joinpath(datadir, "wind.csv"), DataFrame, header=false))'
    local pv_power = Matrix(CSV.read(joinpath(datadir, "solar.csv"), DataFrame, header=false))'
    local load = Matrix(CSV.read(joinpath(datadir, "load.csv"), DataFrame, header=true, transpose=true))[:, 2:end]
    local load_num, t_num  = size(load)
    local wind_num, t_num = size(wind_power)
    local pv_num, t_num = size(pv_power)
    t_num = parse(Int, String(cases_paramsopts["t_num"]))

    if GSw_ExtremeScenarios
        e_num = e_num + e_numx
        println("e_num: $e_num")
        println("e_numx: $e_numx")
    end

    local m = build_model(network, wind_power, pv_power, load, wind_num, pv_num, t_num)
    local solvetime = @elapsed optimize!(m)
    @show termination_status(m)
    diagnose_infeasibility(m)
    @show solution_summary(m; verbose = false)
    @show raw_status(m)
    @show objective_value(m)    
    save_results(m, network, t_num)
end

function build_model(network, wind_power, pv_power, load, wind_num, pv_num, t_num)
    buses = network.buses
    lines = network.lines
    G = network.generators
    n_buses = network.n_buses
    n_lines = network.n_lines
    n_generators = network.n_generators
    slack_bus = network.slack_bus
    global S_base = 1 # NOTE: line characteristics data in `lines.csv` file in the data folder have data that are in per unit is based on S_base = 1 MVA.
    global Δt = 1 # hr # Model resolution
    if GSw_AnnualResolution
        # For run with annual resolution
        Y_epochs = 19
        epoch2year_start = Dict(1 => 1, 2 => 2, 3 => 3, 4 => 4, 5 => 5, 6 => 6, 7 => 7, 8 => 8, 9 => 9, 10 => 10, 11 => 11, 12 => 12, 13 => 13, 14 => 14, 15 => 15, 16 => 16, 17 => 17, 18 => 18, 19 => 19)
        epoch2year_end = Dict(1 => 2, 2 => 3, 3 => 4, 4 => 5, 5 => 6, 6 => 7, 7 => 8, 8 => 9, 9 => 10, 10 => 11, 11 => 12, 12 => 13, 13 => 14, 14 => 15, 15 => 16, 16 => 17, 17 => 18, 18 => 19, 19 => 20)
    else
        # For run with 5-year epoch resolution
        Y_epochs = 4 # [1]: year 1 to 5, [2]: year 6 to 10, [3]: year 11 to 15, [4]: year 16 to 20
        epoch2year_start = Dict(1 => 1, 2 => 6, 3 => 11, 4 => 16)
        epoch2year_end = Dict(1 => 5, 2 => 10, 3 => 15, 4 => 20)
    end

    global y_list = collect(1:Y_epochs) # indices of study years

    # Offshore wind (OSW) units integration specific switch and settings
    n_lines_new_osw = n_lines - n_lines_existing
    n_lines_new = GSw_OSW ? n_lines_new_osw : 0

    if !(n_lines_new_osw in [0, 17, 51]) # NOTE: This part needs to be upated when using model with new sets of data. e.g. in `model_pjm.jl` When running with EPA-derived data this is updated. 1) ISO-NE 3 zone model use [0, 17, 33] 2). 
        error("n_lines_new_osw must be one of [0, 17, 51] to be consistent with input files")
        stop()
    end

    if GSw_OSW & (n_lines_new_osw == 17) # NOTE: This part needs to be upated when using model with new sets of data. e.g. `model_pjm.jl`
        if GSw_AnnualResolution
            wind_online_year =  [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4, 2, 1, 1] # NOTE: This part needs to be upated when using model with new sets of data. e.g. `model_pjm.jl`
        else
            wind_online_year = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] # NOTE: This part needs to be upated when using model with new sets of data. e.g. `model_pjm.jl`
        end
    elseif GSw_OSW & (n_lines_new_osw == 51) # NOTE: This part needs to be upated when using model with new sets of data. e.g. `model_pjm.jl`
        if GSw_AnnualResolution
            wind_online_year = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4, 2, 1] # NOTE: This part needs to be upated when using model with new sets of data. e.g. `model_pjm.jl`
        else
            wind_online_year = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] # NOTE: This part needs to be upated when using model with new sets of data. e.g. `model_pjm.jl`
        end
    end

    _wind_node = zeros(n_buses, e_num, t_num)
    wind_node = zeros(Y, n_buses, e_num, t_num) # wind node power indexed at y, n, e, t
    @info "Y = $Y, n_bus = $n_buses, e_num = $e_num, t_num = $t_num"
    for y = 1:Y,n = 1:n_buses, e = 1:e_num, t = 1:t_num
            if GSw_OSW
                _wind_node[n, e, t] = (y >= wind_online_year[n]) ? (wind_power[n, t+(e-1)*t_num]) : 0
            else
                _wind_node[n, e, t] = wind_power[n, t+(e-1)*t_num]
            end
            wind_node[y, :, :, :] = _wind_node[:, :, :] + (y - 1) * 0.00 * _wind_node[:, :, :] # Existing wind generation capacity stays constant
    end

    global wind_node_out = wind_node[[epoch2year_end[y] for y in y_list], :, :, :]

    _pv_node = zeros(n_buses, e_num, t_num)
    pv_node = zeros(Y, n_buses, e_num, t_num) # pv node power indexed at y, n, e, t
    for y = 1:Y,n = 1:n_buses, e = 1:e_num, t = 1:t_num
            _pv_node[n, e, t] = pv_power[n, t+(e-1)*t_num]
            pv_node[y, :, :, :] = _pv_node[:, :, :] + (y - 1) * 0.00 * _pv_node[:, :, :]
    end

    global pv_node_out = pv_node[[epoch2year_end[y] for y in y_list], :, :, :]

    # calculate net wind and pv injection
    wind_node = wind_node + pv_node

    _load_node = zeros(n_buses, e_num, t_num)
    load_node = zeros(Y, n_buses, e_num, t_num) # demand node power indexed at y, n, e, t
    for y = 1:Y,n = 1:n_buses, e = 1:e_num, t = 1:t_num
            _load_node[n, e, t] = load[n, t+(e-1)*t_num]
            load_node[y, :, :, :] = _load_node[:, :, :] .* ((1 + 0.023)^(y-1)) # "ISO-NE: Overall electricity use is expected to increase, by an average of 2.3% annually"
    end

    global load_node_out = load_node[[epoch2year_end[y] for y in y_list], :, :, :]

    # normalized profiles for new wind and solar
    wind_power_new = Matrix(CSV.read(joinpath(datadir, "wind_normalized.csv"), DataFrame, header=false))'
    pv_power_new = Matrix(CSV.read(joinpath(datadir, "solar_normalized.csv"), DataFrame, header=false))'

    _wind_node_new = zeros(n_buses, e_num, t_num)
    wind_node_new = zeros(Y, n_buses, e_num, t_num) # wind node power indexed at y,n, e, t
    for y = 1:Y,n = 1:n_buses, e = 1:e_num, t = 1:t_num
            _wind_node_new[n, e, t] = wind_power_new[n, t+(e-1)*t_num]
            wind_node_new[y, :, :, :] = _wind_node_new[:, :, :] + (y - 1) * 0.00 * _wind_node_new[:, :, :]
    end

    global wind_node_new_out = wind_node_new[[epoch2year_end[y] for y in y_list], :, :, :]

    _pv_node_new = zeros(n_buses, e_num, t_num)
    pv_node_new = zeros(Y, n_buses, e_num, t_num) # pv node power indexed at y, n, e, t
    for y = 1:Y,n = 1:n_buses, e = 1:e_num, t = 1:t_num
            _pv_node_new[n, e, t] = pv_power_new[n, t+(e-1)*t_num]
            pv_node_new[y, :, :, :] = _pv_node_new[:, :, :] + (y - 1) * 0.00 * _pv_node_new[:, :, :]
    end
    global pv_node_new_out = pv_node_new[[epoch2year_end[y] for y in y_list], :, :, :]

    # ============== #
    # Local switches #
    # ============== #
    # Isolating constraints for for debugging/testing
    GSw_Ramping = true
    GSw_Reserve = true 

    # ============ #
    #  Parameters  #
    # ============ #
    τ = 365 # days
    _VoWS = vcat(VoWS * ones(n_buses_existing), VoOSWS * ones(n_buses - n_buses_existing))
    t_list = collect(1:t_num) # Indices to the time interval, from 1 to T
    s_list = collect(1:n_buses) # Indices to the nodes, from 1 to S
    s_list_existing = s_list[1:n_buses_existing]
    s_list_new = s_list[n_buses_existing+1:end]
    l_list = collect(1:n_lines) # Indices to the transmission lines
    i_list = collect(1:n_generators) # Indices to the generators, from 1 to I
    e_list = collect(1:e_num) # Indices to the representative days
    l_list_candidate = collect(1:n_lines_existing+n_lines_new)
    l_list_existing = collect(1:n_lines_existing)
    l_list_candidate_upgrade = l_list_existing
    l_list_candidate_new = l_list_candidate_upgrade[end] .+ collect(1:n_lines_new)
    c_list = collect(1:3) # Indices to cable type, in increasing discrete sizes

    # =================== # 
    #  Model Definitions  #
    # =================== #
    m = Model(Gurobi.Optimizer)
    # General settings:
    set_optimizer_attribute(m, "OutputFlag", 1) # Display solver logs. Default: 1
    set_optimizer_attribute(m, "Seed", 1234) # Random seed. Default: 0; Set a fixed random seed. The random seed affects the pseudo-random numbers used by the solver.
    # Method and Numeric settings:
    set_optimizer_attribute(m, "Method", 1) # Method. Default: -1=automatic; 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent, and 5=deterministic concurrent simplex.
    # set_optimizer_attribute(m, "BarHomogeneous", 1) # To be used with Barrier algorithm, to deal with the numerical challenges?
    # ObjScale == 1 && set_optimizer_attribute(m, "NumericFocus", 3) # Improve focus on numeric stability. Default: 0. Set NumericFocus parameter to avoid numerical issues
    # set_optimizer_attribute(m, "NumericFocus", 3) # Improve focus on numeric stability. Default: 0. Set NumericFocus parameter to avoid numerical issues
    # set_optimizer_attribute(m, "ObjScale", 1e-3) # Type: double , Default value: 0.0, Minimum value: -1, Maximum value: Infinity. Scale objective function
    # MIP settings:
    set_optimizer_attribute(m, "MIPGap", MIPGap)
    # set_optimizer_attribute(m, "MIPFocus", 0) # MIP Strategy. Default: 0, Quick feasible solutions: 1, Quality solutions to prove optimality: 2, best bound moving very slow or not: 3)
    # set_optimizer_attribute(m, "Heuristics", 0.10) # [0,1]: the desired fraction of total MIP runtime devoted to heuristics (so by default, we aim to spend 5% of runtime on heuristics).
    # set_optimizer_attribute(m, "Cuts", 0) # Disable all cuts (better for determinism, but severely impact performance. Not recommneded)
    # Advanced and tuning settings:
    # set_optimizer_attribute(m, "Threads", 1) # Threads. Limit Gurobi to use only one thread. # Default: 0, generally use as many threads as there are virtual processors.
    # set_optimizer_attribute(m, "Presolve", 0) # Disable presolve (better for determinism, but severely impact performance. Not recommneded). # A value of -1 corresponds to an automatic setting. Other options are off (0), conservative (1), or aggressive (2). More aggressive application of presolve takes more time, but can sometimes lead to a significantly tighter model.
    # set_optimizer_attribute(m, "ScaleFlag", 2) # Scale entire model. 0: Turn off automatic scaling (not recommended), 1 (default): Use Gurobi's standard model scaling, and 2: Use more aggressive model scaling.
    # set_optimizer_attribute(m, "DualReductions", 0) # Default: 1.
    set_optimizer_attribute(m, "FeasibilityTol", 1e-9) # Primal feasibility tolerance.
    set_optimizer_attribute(m, "OptimalityTol", 1e-9) # Dual feasibility tolerance. Default: 1e-6

    # =================== #
    #  Decision variables
    # =================== #
    # Define Binary variables
    @variable(m, Iᴸ[l_list, y_list], Bin) # Line decision, 1 - start construction, 0 - no construction.
    @variable(m, Zᴸ[l_list, y_list], Bin) # New Line status, 1 - in service, 0 - not in service.
    # @variable(m, Sᴸ[l_list, y_list] >=0 ) # fraction of original size
    # # NOTE: This part needs to be upated when using model with new sets of data. e.g. `model_pjm.jl` for EPA-derived data: Below two lines used when transmission line capacity on land are of varying capacity and upgraded as a fraction of existing capacity. When running with EPA-derived data this needs to commented out and these should reflected on new capacity on upgraded lines.
    # @variable(m, Sᴸₓ[1:10], Bin) 
    # @constraint(m, [l=l_list, y=y_list], Sᴸ[l, y] == sum(i * 0.10 * lines[l].s_max * Sᴸₓ[i] for i in 1:10)) # fraction of original size

    # Define Variables
    @variable(m, p[y_list, i_list, e_list, t_list] >= 0) # Total generation of each thermal generator.
    @variable(m, wp[y_list, s_list, e_list, t_list] >= 0) # Wind power generation at year y,  at node s,  in representative day e,  at time stamp t.
    @variable(m, rg[y_list, i_list, e_list, t_list] >= 0) # Reserve provided by a generator
    @variable(m, re[y_list, s_list, e_list, t_list] >= 0) # Reserve provided by a BESS
    @variable(m, rg_new_ng_cc_ccs[y_list, s_list, e_list, t_list] >= 0) # Reserve provided by a new NG-CC-CCS
    @variable(m, rg_new_ng_ct[y_list, s_list, e_list, t_list] >= 0) # Reserve provided by a new NG-CT
    @variable(m, ls[y_list, s_list, e_list, t_list] >= 0) # Load slack
    @variable(m, lsx[y_list, s_list, e_list, t_list] >= 0) # Load slack for extreme days
    @variable(m, ws[y_list, s_list, e_list, t_list] >= 0) # Overgeneration or renewables curtailment slack
    @variable(m, θ[y_list, s_list, e_list, t_list])
    @variable(m, bus_out_power[y_list, s_list, e_list, t_list]) # Power out from bus
    @variable(m, DC_flow[y_list, l_list, e_list, t_list]) # DC line power flows
    @variable(m, new_DC_flow[y_list, l_list, e_list, t_list]) # DC line power flows on new transmission lines
    @variable(m, new_DC_flow_osw[y_list, l_list, e_list, t_list]) # DC line power flows on new transmission lines
    @variable(m, p_new_cap_wind[s_list, y_list] >= 0) # new wind
    @variable(m, p_new_cap_pv[s_list, y_list] >= 0) # new pv
    @variable(m, p_new_cap_ng_cc_ccs[s_list, y_list] >= 0) # New NG-CC-CCS generation capacity
    @variable(m, p_new_ng_cc_ccs[y_list, s_list, e_list, t_list] >= 0) # Generation from new NG-CC-CCS generation capacity
    @variable(m, p_new_cap_ng_ct[s_list, y_list] >= 0) # New NG-CC-CCS generation capacity
    @variable(m, p_new_ng_ct[y_list, s_list, e_list, t_list] >= 0) # Generation from new NG-CT generation capacity
    @variable(m, ch[y_list, s_list, e_list, t_list] >= 0)
    @variable(m, dis[y_list, s_list, e_list, t_list] >= 0)
    @variable(m, SoC[y_list, s_list, e_list, t_list] >= 0)
    @variable(m, Eᵐᵃˣ[s_list, y_list] >= 0) # MWh size of BESS
    @variable(m, Pᵐᵃˣ[s_list, y_list] >= 0) # MW size of BESS
    @variable(m, APe[s_list, y_list] >= 0) # Available power capacity by the endof a year accounting degradation
    @variable(m, AE[s_list, y_list] >= 0) # Available energy capacity by the endof a year accounting degradation
    @variable(m, Eᵐⁱⁿ[s_list, y_list] >= 0) # min MWh size of BESS
    @variable(m, rps_nc[y_list, s_list] >= 0) # Noncompliance of renewable target
    @variable(m, Iᴸ_osw[l_list, c_list, y_list], Bin)
    @variable(m, Zᴸ_osw[l_list, c_list, y_list], Bin)

    # Remove offshore builds to see how configuration changes
    if ~GSw_BatteryOSW
        @constraint(m, [s = s_list_new, y = y_list], Pᵐᵃˣ[s, y] == 0)
        @constraint(m, [s = s_list_new, y = y_list], Eᵐᵃˣ[s, y] == 0)
    end
    # Don't build any new generation resources endogenously on new (offshore) nodes
    @constraint(m, [s = s_list_new, y = y_list], p_new_cap_wind[s, y] == 0)
    @constraint(m, [s = s_list_new, y = y_list], p_new_cap_pv[s, y] == 0)
    @constraint(m, [s = s_list_new, y = y_list], p_new_cap_ng_cc_ccs[s, y] == 0)
    @constraint(m, [s = s_list_new, y = y_list], p_new_cap_ng_ct[s, y] == 0)

    # ================== #
    # Objective function #
    # ================== #
    costs_atb2022 = CSV.read(joinpath(datadir, "costs_atb2022.csv"), DataFrame)
    vom_ng_cc_ccs = costs_atb2022[!, :vom_ng_cc_ccs]

    # Annual operation cost
    vom_ng_cc_ccs = costs_atb2022[!, :vom_ng_cc_ccs] # $/MWh, 2020 USD
    fom_ng_cc_ccs = costs_atb2022[!, :fom_ng_cc_ccs] * 1000 # * 1000 $/kW-yr =>  $/MW-yr
    vom_ng_ct = costs_atb2022[!, :vom_ng_ct] # $/MWh
    fom_ng_ct = costs_atb2022[!, :fom_ng_ct] * 1000 # * 1000 $/kW-yr =>  $/MW-yr
    fom_wind = costs_atb2022[!, :fom_wind] * 1000 # * 1000 $/kW-yr =>  $/MW-yr
    ptc_wind = costs_atb2022[!, :ptc_wind]
    ptc_wind = GSw_WindPTC ? ptc_wind : ptc_wind * 0
    fom_pv = costs_atb2022[!, :fom_pv] * 1000 # * 1000 $/kW-yr =>  $/MW-yr
    itc_pv = costs_atb2022[!, :itc_pv]
    itc_pv = GSw_PVITC ? itc_pv : itc_pv * 0
    fom_battery = costs_atb2022[!, :fom_battery] * 1000 # * 1000 $/kW-yr =>  $/MW-yr

    ωₑ = Matrix(CSV.read(joinpath(datadir, "weights_of_scenarios.csv"),DataFrame, header=false))
    ωₑᴺ = ωₑ[1:Int(length(ωₑ)/2)]
    ωₑᴱ = ωₑ[Int(length(ωₑ)/2)+1:end]

    if GSw_ExtremeScenarios
        if αᴱ == -1.0
            ωₑ = ωₑ
        else
            ωₑᴺ = ωₑᴺ .+ ωₑᴱ
            ωₑᴱ = ωₑᴱ ./= sum(ωₑᴱ)
            ωₑ = vcat((1 - αᴱ) .* ωₑᴺ, αᴱ .* ωₑᴱ)
        end
    else
        ωₑ = ωₑᴺ .+ ωₑᴱ
    end
    τₑ = τ * ωₑ
    if GSw_DemandFlexibility
        # Define parameters related to Demand flexibility model and flexible demand as a decision variable
        c₁ᶠˡᵉˣ = 400 # "Dispatch Cost Coefficient a ($/MWh)" of the most expensive generators
        c₂ᶠˡᵉˣ = 0.0035 # ""Dispatch Cost Coefficient b ($/MW^2h)" of the most expensive generators
        Pₘₐₓ = 300 # Pₘₐₓ of the most expensive generators
        @variable(m, Δd[y_list, s_list, e_list, t_list]) # Flexibility contribution, 0 := no contribution, ℝ₊₊ := flexibility up, ℝ₋₋ := flexibility down
        @variable(m, Δd₁[y_list, s_list, e_list, t_list]) # The least expensive block of flexible demand
        @variable(m, Δd₂[y_list, s_list, e_list, t_list]) # The third most expensive block of flexible demand
        @variable(m, Δd₃[y_list, s_list, e_list, t_list]) # The second most expensive block of flexible demand
        @variable(m, Δd₄[y_list, s_list, e_list, t_list]) # The most expensive block of flexible demand
        # Calulate absolute value of flexible demand blocks
        @variable(m, eΔd₁[y_list, s_list, e_list, t_list] >= 0) # Absolute value of Δd₁
        @variable(m, eΔd₂[y_list, s_list, e_list, t_list] >= 0) # Absolute value of Δd₂
        @variable(m, eΔd₃[y_list, s_list, e_list, t_list] >= 0) # Absolute value of Δd₃
        @variable(m, eΔd₄[y_list, s_list, e_list, t_list] >= 0) # Absolute value of Δd₄
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], eΔd₁[y, s, e, t] >= -Δd₁[y, s, e, t])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], eΔd₁[y, s, e, t] >= Δd₁[y, s, e, t])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], eΔd₂[y, s, e, t] >= -Δd₂[y, s, e, t])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], eΔd₂[y, s, e, t] >= Δd₂[y, s, e, t])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], eΔd₃[y, s, e, t] >= -Δd₃[y, s, e, t])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], eΔd₃[y, s, e, t] >= Δd₃[y, s, e, t])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], eΔd₄[y, s, e, t] >= -Δd₄[y, s, e, t])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], eΔd₄[y, s, e, t] >= Δd₄[y, s, e, t])
        # Divide each flexible demand into four blocks
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], -αᶠˡᵉˣ * load_node[epoch2year_end[y], s, e, t] / S_base / 4 <= Δd₁[y, s, e, t] <= αᶠˡᵉˣ * load_node[epoch2year_end[y], s, e, t] / S_base / 4)
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], -αᶠˡᵉˣ * load_node[epoch2year_end[y], s, e, t] / S_base / 4 <= Δd₂[y, s, e, t] <= αᶠˡᵉˣ * load_node[epoch2year_end[y], s, e, t] / S_base / 4)
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], -αᶠˡᵉˣ * load_node[epoch2year_end[y], s, e, t] / S_base / 4 <= Δd₃[y, s, e, t] <= αᶠˡᵉˣ * load_node[epoch2year_end[y], s, e, t] / S_base / 4)
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], -αᶠˡᵉˣ * load_node[epoch2year_end[y], s, e, t] / S_base / 4 <= Δd₄[y, s, e, t] <= αᶠˡᵉˣ * load_node[epoch2year_end[y], s, e, t] / S_base / 4)
        # Collect all flexible demand blocks
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], Δd[y, s, e, t] == Δd₁[y, s, e, t] + Δd₂[y, s, e, t] + Δd₃[y, s, e, t] + Δd₄[y, s, e, t])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], - αᶠˡᵉˣ * load_node[epoch2year_end[y], s, e, t] / S_base <= Δd[y, s, e, t] <= αᶠˡᵉˣ * load_node[epoch2year_end[y], s, e, t] / S_base)
        # Operation cost with demand flexibility
        @expression(m, eC[y = y_list], (1 / ((1 + r)^(epoch2year_end[y] - 1))) * (sum(τₑ[e] *
                                                                                            (sum(((G[i].c1 * p[y, i, e, t] * S_base)) for i in i_list,t in t_list)
                                                                                            + sum((_VoWS[s] * ws[y, s, e, t] * S_base) for s in s_list,t in t_list)
                                                                                            + sum((VoLL * ls[y, s, e, t] * S_base) for s in s_list,t in t_list)
                                                                                            + (GSw_ExtremeScenarios ? sum((VoLLX * lsx[y, s, e, t] * S_base) for s in s_list,t in t_list; init = 0) : 0)
                                                                                            + sum(vom_ng_cc_ccs[epoch2year_end[r]] * p_new_ng_cc_ccs[r, s, e, t] * S_base for r in 1:y, s in s_list, t in t_list)
                                                                                            + sum(vom_ng_ct[epoch2year_end[r]] * p_new_ng_ct[r, s, e, t] * S_base for r in 1:y, s in s_list, t in t_list)
                                                                                            + sum(p_new_cap_wind[s, y] * S_base * wind_node_new[epoch2year_end[y], s, e, t] * 0 for s in s_list, t in t_list)
                                                                                            - (GSw_WindPTC ? (sum(p_new_cap_wind[s, y] * wind_node_new[epoch2year_end[y], s, e, t] * ptc_wind[y] for s in s_list, t in t_list)) : 0)
                                                                                            + sum(p_new_cap_pv[s, y] * S_base * pv_node_new[epoch2year_end[y], s, e, t] * 0 for s in s_list, t in t_list)
                                                                                            + sum(VoLL * eΔd₄[y, s, e, t] * S_base + (VoLL - (c₁ᶠˡᵉˣ + 2 * c₂ᶠˡᵉˣ * Pₘₐₓ))/4 * eΔd₃[y, s, e, t] * S_base + (VoLL - (c₁ᶠˡᵉˣ + 2 * c₂ᶠˡᵉˣ * Pₘₐₓ))/8 * eΔd₂[y, s, e, t] * S_base + (VoLL - (c₁ᶠˡᵉˣ + 2 * c₂ᶠˡᵉˣ * Pₘₐₓ))/12 * eΔd₁[y, s, e, t] * S_base for s in s_list,t in t_list)) for e in e_list)
                                                                                + sum(fom_ng_cc_ccs[epoch2year_end[y]] * p_new_cap_ng_cc_ccs[s, y] * S_base for s in s_list)
                                                                                + sum(fom_ng_ct[epoch2year_end[y]] * p_new_cap_ng_ct[s, y] * S_base for s in s_list)
                                                                                + sum(fom_wind[epoch2year_end[y]] * p_new_cap_wind[s, y] * S_base for s in s_list)
                                                                                + sum(fom_pv[epoch2year_end[y]] * p_new_cap_pv[s, y] * S_base for s in s_list)
                                                                                + (GSw_Battery ? sum(fom_battery[epoch2year_end[y]] * Pᵐᵃˣ[s, y] * S_base for s in s_list; init = 0) : 0)
                                                                                + sum((RPSNC * rps_nc[y, s] * S_base) for s in s_list)
                                                                                )
        )
        @expression(m, eCOverGenPenalty[y = y_list], (1 / ((1 + r)^(epoch2year_end[y] - 1))) * sum(τₑ[e] * (sum((_VoWS[s] * ws[y, s, e, t] * S_base) for s in s_list,t in t_list)) for e in e_list))
        @expression(m, eCUnderGenPenalty[y = y_list], (1 / ((1 + r)^(epoch2year_end[y] - 1))) * sum(τₑ[e] * (sum((VoLL * ls[y, s, e, t] * S_base + VoLLX * lsx[y, s, e, t] * S_base) for s in s_list,t in t_list)) for e in e_list))
        @expression(m, eCDemandFlexibility[y = y_list], (1 / ((1 + r)^(epoch2year_end[y] - 1))) * sum(τₑ[e] * (sum(VoLL * eΔd₄[y, s, e, t] * S_base + (VoLL - (c₁ᶠˡᵉˣ + 2 * c₂ᶠˡᵉˣ * Pₘₐₓ))/4 * eΔd₃[y, s, e, t] * S_base + (VoLL - (c₁ᶠˡᵉˣ + 2 * c₂ᶠˡᵉˣ * Pₘₐₓ))/8 * eΔd₂[y, s, e, t] * S_base + (VoLL - (c₁ᶠˡᵉˣ + 2 * c₂ᶠˡᵉˣ * Pₘₐₓ))/12 * eΔd₁[y, s, e, t] * S_base for s in s_list,t in t_list)) for e in e_list))
        @expression(m, eCRPSNCPenalty[y = y_list], (1 / ((1 + r)^(epoch2year_end[y] - 1))) * sum((RPSNC * rps_nc[y, s] * S_base) for s in s_list))
    else
        # Operation cost
        @expression(m, eC[y = y_list], (1 / ((1 + r)^(epoch2year_end[y] - 1))) * (sum(τₑ[e] *
                                                                                            (sum(((G[i].c1 * p[y, i, e, t] * S_base)) for i in i_list,t in t_list)
                                                                                            + sum((_VoWS[s] * ws[y, s, e, t] * S_base) for s in s_list,t in t_list)
                                                                                            + sum((VoLL * ls[y, s, e, t] * S_base) for s in s_list,t in t_list)
                                                                                            + (GSw_ExtremeScenarios ? sum((VoLLX * lsx[y, s, e, t] * S_base) for s in s_list,t in t_list; init = 0) : 0)
                                                                                            + sum(vom_ng_cc_ccs[epoch2year_end[r]] * p_new_ng_cc_ccs[r, s, e, t] * S_base for r in 1:y, s in s_list, t in t_list)
                                                                                            + sum(vom_ng_ct[epoch2year_end[r]] * p_new_ng_ct[r, s, e, t] * S_base for r in 1:y, s in s_list, t in t_list)
                                                                                            + sum(p_new_cap_wind[s, y] * S_base * wind_node_new[epoch2year_end[y], s, e, t] * 0 for s in s_list, t in t_list)
                                                                                            - (GSw_WindPTC ? (sum(p_new_cap_wind[s, y] * S_base * wind_node_new[epoch2year_end[y], s, e, t] * ptc_wind[y] for s in s_list, t in t_list)) : 0)
                                                                                            + sum(p_new_cap_pv[s, y] * S_base * pv_node_new[epoch2year_end[y], s, e, t] * 0 for s in s_list, t in t_list)) for e in e_list)
                                                                                + sum(fom_ng_cc_ccs[epoch2year_end[y]] * p_new_cap_ng_cc_ccs[s, y] * S_base for s in s_list)
                                                                                + sum(fom_ng_ct[epoch2year_end[y]] * p_new_cap_ng_ct[s, y] * S_base for s in s_list)
                                                                                + sum(fom_wind[epoch2year_end[y]] * p_new_cap_wind[s, y] * S_base for s in s_list)
                                                                                + sum(fom_pv[epoch2year_end[y]] * p_new_cap_pv[s, y] * S_base for s in s_list)
                                                                                + (GSw_Battery ? sum(fom_battery[epoch2year_end[y]] * Pᵐᵃˣ[s, y] * S_base for s in s_list; init = 0) : 0) 
                                                                                + sum((RPSNC * rps_nc[y, s] * S_base) for s in s_list)
                                                                                )
        )
        @expression(m, eCOverGenPenalty[y = y_list], (1 / ((1 + r)^(epoch2year_end[y] - 1))) * sum(τₑ[e] * (sum((_VoWS[s] * ws[y, s, e, t] * S_base) for s in s_list,t in t_list)) for e in e_list))
        @expression(m, eCUnderGenPenalty[y = y_list], (1 / ((1 + r)^(epoch2year_end[y] - 1))) * sum(τₑ[e] * (sum((VoLL * ls[y, s, e, t] * S_base) for s in s_list,t in t_list)) for e in e_list))
        @expression(m, eCRPSNCPenalty[y = y_list], (1 / ((1 + r)^(epoch2year_end[y] - 1))) * sum((RPSNC * rps_nc[y, s] * S_base) for s in s_list))
    end

    # Annual BESS investment cost #  30-year life (values in 2020$)
    Cᴱ = costs_atb2022[!, :CE] * 1e3 # $/MWh # Battery Energy Capital Cost ($/kWh)
    Cᴾ = costs_atb2022[!, :CP] * 1e3 # $/MW # Battery Power Capital Cost ($/kW)
    kᵉ = Cᴱ * r * (1 + r)^ST / ((1 + r)^ST - 1) # Annualized cost per MWh of BESS
    kᴾ = Cᴾ * r * (1 + r)^ST / ((1 + r)^ST - 1) # Annualized cost per MW of BESS
    γ(n) = (n <= ST) ? 1 / (1 + r)^(n - 1) : 0
    @expression(m, eCBESS_INV[y = y_list], sum(γ(epoch2year_start[y-n+1]) * (sum(kᵉ[epoch2year_end[n]] * Eᵐᵃˣ[s,n] * S_base + kᴾ[epoch2year_end[n]] * Pᵐᵃˣ[s,n] * S_base for s in s_list)) for n in 1:y))

    # Annual line investment cost
    # Paramters above calculated using formula below
    # Cₓ = 1.24301e-7 # Capacitance C (F/mile) of 431.2 MVA cables # xₓ = 1 ./ (2 .* π .* f .* Cₓ .* Lˡ) # f = 60 # Hz
    Lˡ = CSV.read(joinpath(datadir, "lines.csv"), DataFrame)[:, "Distance (miles)"]
    Sˡ = CSV.read(joinpath(datadir, "lines.csv"), DataFrame)[:, "s_max"] # NOTE: This part needs to be upated when using model with new sets of data. e.g. `model_pjm.jl` for EPA-derived data: Use Sᴸ[l, y] instead of Sˡ and comment this line out.
    if GSw_OSW
        # Source: Cost analysis and comparison of HVAC, LFAC and HVDC for offshore wind power connection (https://ieeexplore.ieee.org/document/7787099)
        # 1£ = 1.5285 USD @ 2015, 1 USD at 2015 =  1.09 at 2020, 1 mile = 1.609 km
        # Cˡ_600_HVAC = (0.0119 .* (Lˡ .* 1.609).^2 + 0.8103 .* (Lˡ .* 1.609) .+ 59.41) .* 1e6 .* (1.5285 * 1.09) ./(0.6 * 1000)
        # Cˡ_1400_HVAC = (0.0306 .* (Lˡ .* 1.609).^2 + 0.9835 .* (Lˡ .* 1.609) .+ 165.8) .* 1e6 .* (1.5285 * 1.09) #./(1.4 * 1000)
        # Coefficient-wise interpolation: Cˡ_400_HVAC = (0.007225 .* (Lˡ .* 1.609).^2 + 0.767 .* (Lˡ .* 1.609) .+ 32.8125) .* 1e6 .* (1.5285 * 1.09) ./(0.4 * 1000)
        # Cˡ_600_HVDC = (0.92 .* (Lˡ .* 1.609) .+ 171.48) .* 1e6 .* (1.5285 * 1.09) #./ (0.6 * 1000)
        # Coefficient-wise interpolation: Cˡ_1000_HVDC = (1.14 .* (Lˡ .* 1.609) .+ 269.13) .* 1e6 .* (1.5285 * 1.09) #./ (1 * 1000)
        # Cˡ_1400_HVDC = (1.36 .* (Lˡ .* 1.609) .+ 366.78) .* 1e6 .* (1.5285 * 1.09) #./ (1.4 * 1000)
        # Coefficient-wise interpolation: Cˡ_2200_HVDC = (1.8 .* (Lˡ .* 1.609) .+ 562.08) .* 1e6 .* (1.5285 * 1.09) #./ (2.2 * 1000)
        k_scale_OffshoreLineCost =  0.7340821889389222 # (320000 * 1000 * 1.2772 * 1.01)/Cˡ_1000_HVDC ## https://guidetoanoffshorewindfarm.com/wind-farm-costs $ 1 2020 USD = $1.01 USD, Average exchange rate in 2019: 1.2772 USD, £/MW: 320000, 37.2823 miles (60 km) # (£/MW 320000) * 1000 MW * 1.2772 $/pound * 1.01
        Cˡ_400 = k_scale_OffshoreLineCost * (0.007225 .* (Lˡ .* 1.609).^2 + 0.767 .* (Lˡ .* 1.609) .+ 32.8125) .* 1e6 .* (1.5285 * 1.09) #./(0.4 * 1000)
        Cˡ_1400 = k_scale_OffshoreLineCost * (1.36 .* (Lˡ .* 1.609) .+ 366.78) .* 1e6 .* (1.5285 * 1.09) #./ (1.4 * 1000)
        Cˡ_2200 = k_scale_OffshoreLineCost * (1.8 .* (Lˡ .* 1.609) .+ 562.08) .* 1e6 .* (1.5285 * 1.09) #./ (2.2 * 1000)

        kˡ_400 = Cˡ_400 * r * (1 + r)^LT / ((1 + r)^LT - 1) # Annualized cost per MW of a transmission line
        kˡ_1400 = Cˡ_1400 * r * (1 + r)^LT / ((1 + r)^LT - 1) # Annualized cost per MW of a transmission line
        kˡ_2200 = Cˡ_2200 * r * (1 + r)^LT / ((1 + r)^LT - 1) # Annualized cost per MW of a transmission line
    else
        Cˡ = k_scale_OnshoreLineCost * 3500 * 1.111 * Sˡ[1:n_lines_existing] .* Lˡ[1:n_lines_existing] # https://www.nrel.gov/docs/fy21osti/78195.pdf pp. $1 2013 = $1.11 2020 USD 11.1%
        Lˡ = Lˡ[1:n_lines_existing]
        kˡ = Cˡ * r * (1 + r)^LT / ((1 + r)^LT - 1) # Annualized cost per MW of a transmission line
    end

    # Override cost for first n_lines_existing
    Cˡ = k_scale_OnshoreLineCost * 3500 * 1.111 * Sˡ[1:n_lines_existing] .* Lˡ[1:n_lines_existing] # https://www.nrel.gov/docs/fy21osti/78195.pdf pp. $1 2013 = $1.11 2020 USD 11.1%
    kˡ = Cˡ * r * (1 + r)^LT / ((1 + r)^LT - 1) # Annualized cost per MW of a transmission line
    β(n) = (n <= LT) ? 1 / (1 + r)^(n - 1) : 0
    # NOTE: This part needs to be upated when using model with new sets of data. e.g. `model_pjm.jl` for EPA-derived data:: Use below kˡ[l] * Iᴸ[l, n] * Sᴸ[l, y], and use Lˡ[1:n_lines_existing] instead of Sˡ[1:n_lines_existing] .* Lˡ[1:n_lines_existing] above
    @expression(m, eCL_INV[y = y_list], sum(β(epoch2year_start[y-n+1]) * sum(kˡ[l] * Iᴸ[l, n] for l in l_list_candidate_upgrade) for n in 1:y)
                                        +
                                        (GSw_OSW ? (sum(β(epoch2year_start[y-n+1]) * sum(kˡ_400[l] * Iᴸ_osw[l, 1, n] for l in l_list_candidate_new) for n in 1:y) + sum(β(epoch2year_start[y-n+1]) * sum(kˡ_1400[l] * Iᴸ_osw[l,2,n] for l in l_list_candidate_new) for n in 1:y) + sum(β(epoch2year_start[y-n+1]) * sum(kˡ_2200[l] * Iᴸ_osw[l,3,n] for l in l_list_candidate_new) for n in 1:y)) : 0))

    # Annual generator investment cost
    Cᵖ_ng_cc_ccs = costs_atb2022[!, :Cp_ng_cc_ccs] * 1000 # * 1000 because of CAPEX ($/kW) =>  $/MW
    Cᵖ_ng_ct = costs_atb2022[!, :Cp_ng_ct] * 1000 # * 1000 because of CAPEX ($/kW) =>  $/MW
    Cᵖ_wind = costs_atb2022[!, :Cp_wind] * 1000 #* 1000 because of CAPEX ($/kW) =>  $/MW
    Cᵖ_pv = costs_atb2022[!, :Cp_pv] * 1000 #* 1000 because of CAPEX ($/kW) =>  $/MW
    Cᵖ_pv = -Cᵖ_pv .* (itc_pv .- 1) # ITC for PV
    kᵖ_ng_cc_ccs = Cᵖ_ng_cc_ccs * r * (1 + r)^GT / ((1 + r)^GT - 1)
    kᵖ_ng_ct = Cᵖ_ng_ct * r * (1 + r)^GT / ((1 + r)^GT - 1)
    kᵖ_wind = Cᵖ_wind * r * (1 + r)^GT / ((1 + r)^GT - 1)
    kᵖ_pv = Cᵖ_pv * r * (1 + r)^GT / ((1 + r)^GT - 1)
    α(n) = (n <= GT) ? 1 / (1 + r)^(n - 1) : 0
    @expression(m, eCG_INV[y = y_list], sum(α(epoch2year_start[y-n+1]) * sum(kᵖ_ng_cc_ccs[epoch2year_start[n]] * p_new_cap_ng_cc_ccs[s, n] * S_base for s in s_list) for n in 1:y)
                                        + sum(α(epoch2year_start[y-n+1]) * sum(kᵖ_ng_ct[epoch2year_start[n]] * p_new_cap_ng_ct[s, n] * S_base for s in s_list) for n in 1:y)
                                        + sum(α(epoch2year_start[y-n+1]) * sum(kᵖ_wind[epoch2year_start[n]] * p_new_cap_wind[s, n] * S_base for s in s_list) for n in 1:y)
                                        + sum(α(epoch2year_start[y-n+1]) * sum(kᵖ_pv[epoch2year_start[n]] * p_new_cap_pv[s, n] * S_base for s in s_list) for n in 1:y))

    # Emission (CO2, SO2, NOx, and PM2.5) tracking and pricing (CO2)
    # Emission factor for new natural gas plants # Averaged historical emission by zone
    emission_factor_CO2_ng = 0.725072727
    emission_factor_SO2_ng = 0.001663782
    emission_factor_NOx_ng = 0.019212469
    emission_factor_PM25_ng = 0.0000251682
    # Pricing Externalities such as CO2 emission and air quality damages
    @expression(m, eCE_Emission[y = y_list], (1 / ((1 + r)^(epoch2year_end[y] - 1))) * sum(τₑ[e] * (sum(P_CO2 * G[i].CO2mtonperMWh * p[y, i, e, t] * S_base for i in i_list,t in t_list)
                                                                                            +
                                                                                            sum(P_CO2 * emission_factor_CO2_ng * p_new_ng_cc_ccs[y, s, e, t] * S_base for s in s_list,t in t_list)
                                                                                            +
                                                                                            sum(P_CO2 * emission_factor_CO2_ng * p_new_ng_ct[y, s, e, t] * S_base for s in s_list,t in t_list)
                                                                                            )
                                                                                    for e in e_list))
    @expression(m, eCE_AirQuality[y = y_list], (1 / ((1 + r)^(epoch2year_end[y] - 1))) * sum(τₑ[e] * (sum(G[i].ave_marg_damages_isrm_LePeule * 1.08 * k_scale_AQ * p[y, i, e, t] * S_base for i in i_list,t in t_list) # 1 USD 2016 = USD 1.08 2020
                                                                                        +
                                                                                        sum(22.415 * 1.08 * k_scale_AQ * p_new_ng_cc_ccs[y, s, e, t] * S_base for s in s_list,t in t_list) # 22.415 $/MWh ave_marg_damages_isrm_LePeule for new NG-CC-CCS units
                                                                                        +
                                                                                        sum(26.22 * 1.08 * k_scale_AQ * p_new_ng_ct[y, s, e, t] * S_base for s in s_list,t in t_list) # 26.22 $/MWh ave_marg_damages_isrm_LePeule for new NG-CT units
                                                                                        )
                                                                                for e in e_list))
    # Tracking Global (CO2) and Local Pollutants (SO2, NOx, PM2.5)
    @expression(m, eE_CO2[y = y_list], sum(τₑ[e] * (sum((G[i].CO2mtonperMWh * p[y, i, e, t] * S_base) for i in i_list,t in t_list) + emission_factor_CO2_ng * sum(p_new_ng_cc_ccs[r, s, e, t] * S_base for r in 1:y, s in s_list, t in t_list) + emission_factor_CO2_ng * sum(p_new_ng_ct[r, s, e, t] * S_base for r in 1:y, s in s_list, t in t_list)) for e in e_list))
    @expression(m, eE_SO2[y = y_list], sum(τₑ[e] * (sum((G[i].SO2mtonperMWh * p[y, i, e, t] * S_base) for i in i_list,t in t_list) + emission_factor_SO2_ng * sum(p_new_ng_cc_ccs[r, s, e, t] * S_base for r in 1:y, s in s_list, t in t_list) + emission_factor_SO2_ng * sum(p_new_ng_ct[r, s, e, t] * S_base for r in 1:y, s in s_list, t in t_list)) for e in e_list))
    @expression(m, eE_NOx[y = y_list], sum(τₑ[e] * (sum((G[i].NOxmtonperMWh * p[y, i, e, t] * S_base) for i in i_list,t in t_list) + emission_factor_NOx_ng * sum(p_new_ng_cc_ccs[r, s, e, t] * S_base for r in 1:y, s in s_list, t in t_list) + emission_factor_NOx_ng * sum(p_new_ng_ct[r, s, e, t] * S_base for r in 1:y, s in s_list, t in t_list)) for e in e_list))
    @expression(m, eE_PM25[y = y_list], sum(τₑ[e] * (sum((G[i].PM25mtonperMWh * p[y, i, e, t] * S_base) for i in i_list,t in t_list) + emission_factor_PM25_ng * sum(p_new_ng_cc_ccs[r, s, e, t] * S_base for r in 1:y, s in s_list, t in t_list) + emission_factor_PM25_ng * sum(p_new_ng_ct[r, s, e, t] * S_base for r in 1:y, s in s_list, t in t_list)) for e in e_list))

    # Set the objective
    if GSw_MultiObj
        @objective(m, Min, sum(eC[y] * ObjScale + eCL_INV[y] * ObjScale + eCG_INV[y] * ObjScale + (GSw_Battery ? eCBESS_INV[y] * ObjScale : 0) + (GSw_Emission ? eCE_Emission[y] * ObjScale : 0) + (GSw_AirQuality ? eCE_AirQuality[y] * ObjScale : 0) for y in y_list))
    else
        @objective(m, Min, sum(eC[y] * ObjScale + eCL_INV[y] * ObjScale + eCG_INV[y] * ObjScale + (GSw_Battery ? eCBESS_INV[y] * ObjScale : 0) for y in y_list))
    end

    # ======== Constraints ======== #
    # Generator constraints
    @variable(m, z[y in y_list, i in i_list], Bin)
    RUrate = 400 # MW/h # NOTE: You can update this data for new builds
    for y in y_list
        for e in e_list, t in t_list, i in i_list
            fix(z[y, i], 1.0)
            @constraint(m, p[y, i, e, t] >= 0 / S_base * z[y, i])
            @constraint(m, p[y, i, e, t] <= G[i].Pmax / S_base * z[y, i])
        end
        for e in e_list, t in t_list, s in s_list
            @constraint(m, p_new_ng_cc_ccs[y, s, e, t] >= 0)
            @constraint(m, p_new_ng_cc_ccs[y, s, e, t] <= sum(p_new_cap_ng_cc_ccs[s, r] / S_base for r ∈ 1:y))
            @constraint(m, p_new_ng_ct[y, s, e, t] >= 0)
            @constraint(m, p_new_ng_ct[y, s, e, t] <= sum(p_new_cap_ng_ct[s, r] / S_base for r ∈ 1:y))
        end
        if GSw_Ramping
            for e in e_list, t in t_list[2:end], i in i_list
                @constraint(m, (p[y, i, e, t] + rg[y, i, e, t] - p[y, i, e, t-1]) <= G[i].RUrate / S_base * z[y, i])
                @constraint(m, -(p[y, i, e, t] - rg[y, i, e, t-1] - p[y, i, e, t-1]) <= G[i].RDrate / S_base * z[y, i])
            end
            for e in e_list, t in t_list[2:end], s in s_list
                @constraint(m, (p_new_ng_cc_ccs[y, s, e, t] + rg_new_ng_cc_ccs[y, s, e, t] - p_new_ng_cc_ccs[y, s, e, t-1]) <= RUrate / S_base)
                @constraint(m, -(p_new_ng_cc_ccs[y, s, e, t] - rg_new_ng_cc_ccs[y, s, e, t-1] - p_new_ng_cc_ccs[y, s, e, t-1]) <= RUrate / S_base)
                @constraint(m, (p_new_ng_ct[y, s, e, t] + rg_new_ng_ct[y, s, e, t] - p_new_ng_ct[y, s, e, t-1]) <= RUrate / S_base)
                @constraint(m, -(p_new_ng_ct[y, s, e, t] - rg_new_ng_ct[y, s, e, t-1] - p_new_ng_ct[y, s, e, t-1]) <= RUrate / S_base)
            end
        end
    end

    # Exogenous Retirements
    if GSw_ExogenousRetirements
        # Retire all plants that are already retired 
        for i in i_list, y in y_list, e in e_list, t in t_list
            if (G[i].PlannedRetirementYear <= Y_start)
                fix(z[y, i], 0.0; force=true)
            end
        end
        # Retire units that are planned to retire within planing horizon on planned retirement date and turn always on
        # y_epoch = (y_annual ÷ 5  + 1) # y_annual = (G[i].PlannedRetirementYear - 2022) # ÷ calculates quotient
        for i in i_list, y in y_list, e in e_list, t in t_list
            if (G[i].PlannedRetirementYear > Y_start) & (G[i].PlannedRetirementYear <= Y_start + Y)
                if y < ((G[i].PlannedRetirementYear - (Y_start-1)) ÷ (Y_epochs+1)  + 1)
                    # pass
                else
                    fix(z[y, i], 0.0; force=true)
                end
            end
        end
    end

    # Reserve constraints
    for y in y_list
        for e in e_list, t in t_list, i in i_list
            @constraint(m, p[y, i, e, t] + rg[y, i, e, t] <= G[i].Pmax / S_base * z[y, i])
            if GSw_Ramping
                @constraint(m, rg[y, i, e, t] <= G[i].RUrate / S_base * δt * z[y, i])
            end
        end
        for e in e_list, t in t_list, s in s_list
            @constraint(m, p_new_ng_cc_ccs[y, s, e, t] + rg_new_ng_cc_ccs[y, s, e, t] <= sum(p_new_cap_ng_cc_ccs[s, r] / S_base for r ∈ 1:y))
            if GSw_Ramping
                @constraint(m, rg_new_ng_cc_ccs[y, s, e, t] <= RUrate / S_base * δt)
            end
            @constraint(m, p_new_ng_ct[y, s, e, t] + rg_new_ng_ct[y, s, e, t] <= sum(p_new_cap_ng_ct[s, r] / S_base for r ∈ 1:y))
            if GSw_Ramping
                @constraint(m, rg_new_ng_ct[y, s, e, t] <= RUrate / S_base * δt)
            end
        end
    end

    if GSw_Battery
        # BESS constraints
        @constraint(m, [s = s_list, y = y_list], Eᵐⁱⁿ[s, y] == 0.2 * Eᵐᵃˣ[s, y]) # 20%
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list[2:end]], SoC[y, s, e, t] == SoC[y, s, e, t-1] +  ηₛᶜʰ * ch[y, s, e, t] * Δt - dis[y, s, e, t] * Δt / ηₛᵈⁱˢ) # duration of each period = 1
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list[1]], SoC[y, s, e, t] == 0.80 * AE[s, y]) # Head boundary condition # SoC⁰ = 80% or 50% or 100% depending on the utilization
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list[end]], SoC[y, s, e, t] == 0.80 * AE[s, y]) # Tail boundary condition # SoC⁰ = 80% or 50% or 100% depending on the utilization
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], Eᵐⁱⁿ[s, y] <= SoC[y, s, e, t])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], SoC[y, s, e, t] <= AE[s, y])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], ch[y, s, e, t] <= APe[s, y])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], dis[y, s, e, t] <= APe[s, y])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], re[y, s, e, t] + dis[y, s, e, t] - ch[y, s, e, t] <= APe[s, y])
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], (re[y, s, e, t] - ch[y, s, e, t]) * δt <= ηₛᵈⁱˢ * (SoC[y, s, e, t] - Eᵐⁱⁿ[s,y]))

        # Battery storage planning
        CF = 1 - k_cal - k_cycle
        a(n) = (n <= ST) ? 1 * (CF)^(n - 1) : 0
        @constraint(m, [s = s_list, y = y_list], AE[s, y] == sum(a(epoch2year_start[y-t₀+1]) * Eᵐᵃˣ[s, t₀] for t₀ in 1:y))
        b(n) = (n <= ST) ? 1 : 0
        @constraint(m, [s = s_list, y = y_list], APe[s, y] == sum(b(epoch2year_start[y-t₀+1]) * Pᵐᵃˣ[s, t₀] for t₀ in 1:y))
        @constraint(m, [s = s_list, y = y_list], battery_size_hr * Pᵐᵃˣ[s, y] == Eᵐᵃˣ[s, y])
    end

    # Power balance constraints
    if ~GSw_ExtremeScenarios
        @info "Power balance when e_numx == 0: normal e_list = $e_list"
        @constraint(m, λ[y = y_list, s = s_list, e = e_list, t = t_list], p_new_ng_cc_ccs[y, s, e, t] + p_new_ng_ct[y, s, e, t] + sum(p[y, g, e, t] for g in buses[s].generator) + (wind_node[epoch2year_end[y], s, e, t] / S_base - ws[y, s, e, t]) + sum(p_new_cap_wind[s, r] * wind_node_new[epoch2year_end[y], s, e, t] / S_base + p_new_cap_pv[s, r] * pv_node_new[epoch2year_end[y], s, e, t] / S_base for r in 1:y) == load_node[epoch2year_end[y], s, e, t] / S_base + (GSw_DemandFlexibility ? Δd[y, s, e, t] : 0) - ls[y, s, e, t] + bus_out_power[y, s, e, t] + (GSw_Battery ? (ch[y, s, e, t] - dis[y, s, e, t]) : 0))
        # Undergeneration slack
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list],  ls[y, s, e, t] <= load_node[epoch2year_end[y], s, e, t] / S_base + (GSw_DemandFlexibility ? Δd[y, s, e, t] : 0))
    else
        @info "Power balance e_numx == $e_numx: normal e_list = $(e_list[1:e_num-e_numx])"
        @info "Power balance e_numx == $e_numx: extreme e_list = $(e_list[e_num-e_numx+1:end])"
        println("e_num: $e_num")
        println("e_numx: $e_numx")
        @constraint(m, λ[y = y_list, s = s_list, e = e_list[1:e_num-e_numx], t = t_list], p_new_ng_cc_ccs[y, s, e, t] + p_new_ng_ct[y, s, e, t] + sum(p[y, g, e, t] for g in buses[s].generator) + (wind_node[epoch2year_end[y], s, e, t] / S_base - ws[y, s, e, t]) + sum(p_new_cap_wind[s, r] * wind_node_new[epoch2year_end[y], s, e, t] / S_base + p_new_cap_pv[s, r] * pv_node_new[epoch2year_end[y], s, e, t] / S_base for r in 1:y) == load_node[epoch2year_end[y], s, e, t] / S_base + (GSw_DemandFlexibility ? Δd[y, s, e, t] : 0) - ls[y, s, e, t] + bus_out_power[y, s, e, t] + (GSw_Battery ? (ch[y, s, e, t] - dis[y, s, e, t]) : 0))
        @constraint(m, λᴱ[y = y_list, s = s_list, e = e_list[e_num-e_numx+1:end], t = t_list], p_new_ng_cc_ccs[y, s, e, t] + p_new_ng_ct[y, s, e, t] + sum(p[y, g, e, t] for g in buses[s].generator) + (wind_node[epoch2year_end[y], s, e, t] / S_base - ws[y, s, e, t]) + sum(p_new_cap_wind[s, r] * wind_node_new[epoch2year_end[y], s, e, t] / S_base + p_new_cap_pv[s, r] * pv_node_new[epoch2year_end[y], s, e, t] / S_base for r in 1:y) == load_node[epoch2year_end[y], s, e, t] / S_base + (GSw_DemandFlexibility ? Δd[y, s, e, t] : 0) - lsx[y, s, e, t] + bus_out_power[y, s, e, t] + (GSw_Battery ? (ch[y, s, e, t] - dis[y, s, e, t]) : 0)) 
        # Undergeneration slack
        @constraint(m, [y = y_list, s = s_list, e = e_list[1:e_num-e_numx], t = t_list],  ls[y, s, e, t] <= load_node[epoch2year_end[y], s, e, t] / S_base + (GSw_DemandFlexibility ? Δd[y, s, e, t] : 0))
        @constraint(m, [y = y_list, s = s_list, e = e_list[e_num-e_numx+1:end], t = t_list], lsx[y, s, e, t] <= load_node[epoch2year_end[y], s, e, t] / S_base + (GSw_DemandFlexibility ? Δd[y, s, e, t] : 0))
    end

    @expression(m, p_new_cap_wind_annual[s=s_list,y=y_list], sum(p_new_cap_wind[s, r] for r in 1:y))
    @expression(m, p_new_cap_pv_annual[s=s_list,y=y_list], sum(p_new_cap_pv[s, r] for r in 1:y))

    if GSw_DemandFlexibility
        @constraint(m, [y = y_list, s = s_list, e = e_list], sum(Δd[y, s, e, t] for t in t_list) == 0)
    end
    # Overganeration slack # Renewables Curtailment <= 50% * Renewable Generation
    @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list],  (1/α⁺) * ws[y, s, e, t] <= (wind_node[epoch2year_end[y], s, e, t] / S_base) + sum(p_new_cap_wind[s, r] * wind_node_new[epoch2year_end[y], s, e, t] / S_base + p_new_cap_pv[s, r] * pv_node_new[epoch2year_end[y], s, e, t] / S_base for r in 1:y))

    # Reserve requirements
    if GSw_Reserve
        for y in y_list, e in e_list, t in t_list
            @constraint(m, sum(rg[y, i, e, t] for i in i_list) + sum(rg_new_ng_cc_ccs[y, s, e, t] for s in s_list) + sum(rg_new_ng_ct[y, s, e, t] for s in s_list) + (GSw_Battery ? sum(re[y, s, e, t] for s in s_list; init = 0) : 0) >= 0.03 * sum(load_node[epoch2year_end[y], :, e, t] / S_base) + 0.05 * sum(wind_node[epoch2year_end[y], :, e, t] / S_base) + 0.05 * sum(p_new_cap_wind[s, r] * wind_node_new[epoch2year_end[y], s, e, t] / S_base + p_new_cap_pv[s, r] * pv_node_new[epoch2year_end[y], s, e, t] / S_base for r in 1:y, s in s_list))
        end
    else
        @constraint(m, [y = y_list, i = i_list, e = e_list, t = t_list], rg[y, i, e, t] == 0) # Reserve provided by a generator
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], re[y, s, e, t] == 0) # Reserve provided by a BESS
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], rg_new_ng_cc_ccs[y, s, e, t]== 0) # Reserve provided by a new NG-CC-CCS
        @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], rg_new_ng_ct[y, s, e, t] == 0) # Reserve provided by a new NG-CT
    end

    # Renewable Portfolio Standards (RPS)
    # ME 100% by 2050, VT 75% by 2032, NH 25.2% by 2025, MA 41.1% by 2030, CT 44% by 2030, RI 38.5% by 2035 # https://ieeexplore.ieee.org/document/9428503
    # ME 80% by 2030, NH 25.2% by 2025, VT 75% by 2032, MA 80% by clean energy 2030, CT 48% by 2030, RI 38.5% by 2035 # https://www.spglobal.com/marketintelligence/en/news-insights/research/new-england-renewable-policies-to-drive-12500-mw-of-renewable-capacity-by-2030
    # zone_zonenum = Dict("ME" => 1, "NH" => 2, "VT" => 3, "WCMASS" => 4, "WCMA" => 4, "NEMA" => 5, "NEMA/BOST" => 5, "CT" => 6, "RI" => 7,"SEMASS" => 8, "SEMA" => 8)
    # RPS by years # RPS = Dict("ME" => (2030, 0.80), "NH" => (2025, 0.252), "VT" => (2032, 0.75),"MA" => (2030, 0.411),"CT" => (2030, 0.48),"RI" => (2035, 0.385)) # RPS by epochs, assuming 2023 as year 1
    RPS = Dict("ME" => (2, 0.80), "NH" => (1, 0.252), "VT" => (2, 0.75),"MA" => (2, 0.35),"CT" => (2, 0.48),"RI" => (3, 0.385)) # 35% by 2030 and an additional 1% each year after
    state_zone = Dict("ME" => [1], "NH" => [2], "VT" => [3],"MA" => [4, 5, 8],"CT" => [6],"RI" => [7])

    if GSw_RPS
        for state in keys(RPS)
                @constraint(m, [y = RPS[state][1]:y_list[end], _s=state_zone[state]], sum(sum((GSw_OSW ? sum(new_DC_flow_osw[y, _l, e, t] for _l ∈ intersect(l_list_candidate_new, buses[s].inlist); init = 0) : 0) + wind_node[epoch2year_end[y], s, e, t] / S_base + sum(p_new_cap_wind[s, r] * wind_node_new[epoch2year_end[y], s, e, t] / S_base + p_new_cap_pv[s, r] * pv_node_new[epoch2year_end[y], s, e, t] / S_base for r in 1:y) - ws[y, s, e, t] for e in e_list, t in t_list) for s in _s) + (GSw_RPSNC ? sum(rps_nc[y, s] for s in _s) : 0) >= RPS[state][2] * sum(sum(load_node[epoch2year_end[y], s, e, t] / S_base + (GSw_DemandFlexibility ? Δd[y, s, e, t] : 0) for e in e_list, t in t_list) for s in _s))
        end
    end

    # Line flows and transmission capacity (for existing lines)
    @constraint(m, [y = y_list, l = l_list_existing, e = e_list, t = t_list], DC_flow[y, l, e, t] * lines[l].x == (θ[y, lines[l].from_node, e, t] - θ[y, lines[l].to_node, e, t]))
    @constraint(m, [y = y_list, l = l_list_existing, e = e_list, t = t_list], -lines[l].s_max / S_base <= DC_flow[y, l, e, t] <= lines[l].s_max / S_base)
    @constraint(m, [y = y_list, s = s_list, e = e_list, t = t_list], -π <= θ[y, s, e, t] <= π)
    @constraint(m, [y = y_list, s = slack_bus, e = e_list, t = t_list], θ[y, s, e, t] == 0)
    for s in s_list
        if buses[s].inlist == [] && buses[s].outlist == []
            @constraint(m, [y = y_list, e = e_list, t = t_list], bus_out_power[y, s, e, t] == 0)
        elseif buses[s].inlist != [] && buses[s].outlist == []
            @constraint(m, [y = y_list, e = e_list, t = t_list], bus_out_power[y, s, e, t] == sum(-DC_flow[y, k, e, t] for k in intersect(buses[s].inlist, l_list_existing); init=0) + sum(-new_DC_flow[y, k, e, t] for k in intersect(buses[s].inlist, l_list_candidate_upgrade); init=0) + (GSw_OSW ? sum(-new_DC_flow_osw[y, k, e, t] for k in intersect(buses[s].inlist, l_list_candidate_new); init=0) : 0))
        elseif buses[s].inlist == [] && buses[s].outlist != []
            @constraint(m, [y = y_list, e = e_list, t = t_list], bus_out_power[y, s, e, t] == sum(DC_flow[y, k, e, t] for k in intersect(buses[s].outlist, l_list_existing); init=0) + sum(new_DC_flow[y, k, e, t] for k in intersect(buses[s].outlist, l_list_candidate_upgrade); init=0) + (GSw_OSW ? sum(new_DC_flow_osw[y, k, e, t] for k in intersect(buses[s].outlist, l_list_candidate_new); init=0) : 0))
        elseif buses[s].inlist != [] && buses[s].outlist != []
            @constraint(m, [y = y_list, e = e_list, t = t_list], bus_out_power[y, s, e, t] == sum(-DC_flow[y, k, e, t] for k in intersect(buses[s].inlist, l_list_existing); init=0) + sum(DC_flow[y, k, e, t] for k in intersect(buses[s].outlist, l_list_existing); init=0) + sum(-new_DC_flow[y, k, e, t] for k in intersect(buses[s].inlist, l_list_candidate_upgrade); init=0) + sum(new_DC_flow[y, k, e, t] for k in intersect(buses[s].outlist, l_list_candidate_upgrade); init=0) + (GSw_OSW ? (sum(-new_DC_flow_osw[y, k, e, t] for k in intersect(buses[s].inlist, l_list_candidate_new); init=0) + sum(new_DC_flow_osw[y, k, e, t] for k in intersect(buses[s].outlist, l_list_candidate_new); init=0)) : 0))
        end
    end

    # Line planning and operation
    # Allow only one time upgrade of a candidate line over the planning horizon
    @constraint(m, [l = l_list_candidate_upgrade], 0 <= sum(Iᴸ[l, r] for r in y_list) <= 1)

    @constraint(m,[l = l_list_candidate_upgrade, y = y_list], Zᴸ[l, y] == sum(Iᴸ[l, r] for r in 1:y)) # No delay
    # @constraint(m,[l = l_list_candidate_upgrade, y = y_list], Zᴸ[l, y] == sum(Iᴸ[l, r] for r in 1:y-1)) # Delay of 1 year
    @constraint(m,[y = y_list[2:end], l = l_list_candidate_upgrade], Zᴸ[l, y] >= Zᴸ[l, y-1])

    # # Assuming all 12 lines as candidate upgrades
    # Munoz et al.: An ideal choice of Mˡ is equal to the maximum angle difference times the susceptance of candidate line l, implemented below as: π / 2 * 1 / lines[l].x # Additional ref. : https://dl.acm.org/doi/10.1145/3396851.3397688
    @constraint(m, [y = y_list, l = l_list_candidate_upgrade, e = e_list, t = t_list], - (π / 2 * 1 / lines[l].x * (1 - Zᴸ[l, y])) <= new_DC_flow[y, l, e, t] * lines[l].x - (θ[y, lines[l].from_node, e, t] - θ[y, lines[l].to_node, e, t]))
    @constraint(m, [y = y_list, l = l_list_candidate_upgrade, e = e_list, t = t_list],  new_DC_flow[y, l, e, t] * lines[l].x - (θ[y, lines[l].from_node, e, t] - θ[y, lines[l].to_node, e, t]) <= (π / 2 * 1 / lines[l].x * (1 - Zᴸ[l, y])))
    @constraint(m, [y = y_list, l = l_list_candidate_upgrade, e = e_list, t = t_list], -lines[l].s_max / S_base * Zᴸ[l, y] <= new_DC_flow[y, l, e, t]) # NOTE: Use Sᴸ[l, y] instead of lines[l].s_max
    @constraint(m, [y = y_list, l = l_list_candidate_upgrade, e = e_list, t = t_list], new_DC_flow[y, l, e, t] <= lines[l].s_max / S_base * Zᴸ[l, y]) # NOTE: Use Sᴸ[l, y] instead of lines[l].s_max

    # New line builds for Offshore wind integration
    if GSw_OSW
        @constraint(m, [l = l_list_candidate_new], sum(Iᴸ_osw[l, c, r] for r in y_list, c in c_list) <= 1)

        @constraint(m,[l = l_list_candidate_new, y = y_list, c = c_list], Zᴸ_osw[l, c, y] == sum(Iᴸ_osw[l, c, r] for r in 1:y)) # No delay
        # @constraint(m,[l = l_list_candidate_new,y = y_list], Zᴸ_osw[l, c, y] == sum(Iᴸ_osw[l, c, r] for r in 1:y-1)) # Delay of 1 year
        @constraint(m,[y = y_list[2:end],c = c_list, l = l_list_candidate_new], Zᴸ_osw[l, c, y] >= Zᴸ_osw[l, c, y-1])

        # Munoz et al.: An ideal choice of Mˡ is equal to the maximum angle difference times the susceptance of candidate line l, implemented below as: π / 2 * 1 / lines[l].x # Additional ref. : https://dl.acm.org/doi/10.1145/3396851.3397688
        @constraint(m, [y = y_list, l = l_list_candidate_new, c = c_list[1], e = e_list, t = t_list], - (π / 2 * 1 / lines[l].x * (1 - Zᴸ_osw[l, 1, y])) <= new_DC_flow_osw[y, l, e, t] * lines[l].x - (θ[y, lines[l].from_node, e, t] - θ[y, lines[l].to_node, e, t]))
        @constraint(m, [y = y_list, l = l_list_candidate_new, c = c_list[1], e = e_list, t = t_list],  new_DC_flow_osw[y, l, e, t] * lines[l].x - (θ[y, lines[l].from_node, e, t] - θ[y, lines[l].to_node, e, t]) <= (π / 2 * 1 / lines[l].x * (1 - Zᴸ_osw[l, 1, y])))
        @constraint(m, [y = y_list, l = l_list_candidate_new, e = e_list, t = t_list], - (400 / S_base * Zᴸ_osw[l, 1, y] + 1400 / S_base * Zᴸ_osw[l, 2, y] + 2200 / S_base * Zᴸ_osw[l, 3, y]) <= new_DC_flow_osw[y, l, e, t])
        @constraint(m, [y = y_list, l = l_list_candidate_new, e = e_list, t = t_list], new_DC_flow_osw[y, l, e, t] <= (400 / S_base * Zᴸ_osw[l, 1, y] + 1400 / S_base * Zᴸ_osw[l, 2, y] + 2200 / S_base * Zᴸ_osw[l, 3, y]))
        
        # Build at least one line on the online year and let the model choose which configuration to build and how and where to connect to land.
        for s in s_list_new
            if ~isempty(intersect(union(buses[s].inlist, buses[s].outlist), l_list_candidate_new))
                    @constraint(m, sum(Iᴸ_osw[k, c, wind_online_year[s]] for k in intersect(union(buses[s].inlist, buses[s].outlist), l_list_candidate_new), c in c_list) >= 1)
            end
        end
    end
    return m
end

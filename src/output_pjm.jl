function save_results(m, network, t_num)
    G = network.generators
    n_generators = network.n_generators
    n_buses = network.n_buses
    n_lines = network.n_lines

    t_list = collect(1:t_num)
    s_list = collect(1:n_buses)
    e_list = collect(1:e_num)

    function _reshape_line_ylet_el_yt(mat4d)
        mat2d = zeros(n_lines*e_num, length(y_list)*t_num)
        for e=1:e_num, l=1:n_lines
            for y=y_list, t=1:t_num
                mat2d[l+n_lines*(e-1), t+t_num*(y-1)] = mat4d[y, l, e, t]
            end
        end
        return mat2d
    end

    function _reshape_node_yset_es_yt(mat4d)
        mat2d = zeros(n_buses*e_num, length(y_list)*t_num)
        for e=1:e_num, s=1:n_buses
            for y=y_list, t=1:t_num
                mat2d[s+n_buses*(e-1), t+t_num*(y-1)] = mat4d[y, s, e, t]
            end
        end
        return mat2d
    end

    function _reshape_line_lcy_l_cy(mat3d)
        mat2d = zeros(n_lines, length(y_list)*3)
        for l=1:n_lines
            for y=y_list, c=1:3
                mat2d[l, c+3*(y-1)] = mat3d[l, c, y]
            end
        end
        return mat2d
    end

    Pg = value.(m[:p]) * S_base
    Pg2d = zeros(n_generators*e_num, length(y_list)*t_num)
    for e=1:e_num, i=1:n_generators
        for y=y_list, t=1:t_num
            Pg2d[i+n_generators*(e-1), t+t_num*(y-1)] = Pg[y, i, e, t]
        end
    end
    CSV.write(joinpath(output_dir, "power_output_existing_generators.csv"), Tables.table(Pg2d * S_base), header=false)
    CSV.write(joinpath(output_dir, "power_output_new_ng_cc_ccs.csv"), Tables.table(_reshape_node_yset_es_yt(value.(m[:p_new_ng_cc_ccs]) * S_base)), header=false)
    CSV.write(joinpath(output_dir, "power_output_new_ng_ct.csv"), Tables.table(_reshape_node_yset_es_yt(value.(m[:p_new_ng_ct]) * S_base)), header=false)
    CSV.write(joinpath(output_dir, "load_node_out.csv"), Tables.table(_reshape_node_yset_es_yt(load_node_out) * S_base), header=false)
    CSV.write(joinpath(output_dir, "wind_node_out.csv"), Tables.table(_reshape_node_yset_es_yt(wind_node_out) * S_base), header=false)
    CSV.write(joinpath(output_dir, "pv_node_out.csv"), Tables.table(_reshape_node_yset_es_yt(pv_node_out) * S_base), header=false)
    CSV.write(joinpath(output_dir, "wind_node_new_out_raw.csv"), Tables.table(_reshape_node_yset_es_yt(wind_node_new_out) * S_base), header=false)
    CSV.write(joinpath(output_dir, "wind_node_new_out.csv"), Tables.table(_reshape_node_yset_es_yt(wind_node_new_out)  * S_base .* repeat(repeat(value.(m[:p_new_cap_wind_annual])', inner=(t_num,1))',e_num)), header=false)
    CSV.write(joinpath(output_dir, "pv_node_new_out_raw.csv"), Tables.table(_reshape_node_yset_es_yt(pv_node_new_out) * S_base), header=false)
    CSV.write(joinpath(output_dir, "pv_node_new_out.csv"), Tables.table(_reshape_node_yset_es_yt(pv_node_new_out) * S_base .* repeat(repeat(value.(m[:p_new_cap_pv_annual])', inner=(t_num,1))',e_num)), header=false)

    Pg2d_node = zeros(n_buses*e_num, length(y_list)*t_num)
    for e=1:e_num, s=1:n_buses
        for y=y_list, t=1:t_num
            Pg2d_node[s+n_buses*(e-1), t+t_num*(y-1)] = isempty(network.buses[s].generator) ? 0 : sum(Pg[y, i, e, t] for i in network.buses[s].generator)
        end
    end
    CSV.write(joinpath(output_dir, "power_output_existing_generators_node.csv"), Tables.table(Pg2d_node * S_base), header=false)

    CSV.write(joinpath(output_dir, "load_curtailments_normal.csv"), Tables.table(_reshape_node_yset_es_yt(value.(m[:ls]) * S_base)), header=false)
    CSV.write(joinpath(output_dir, "load_curtailments_extreme.csv"), Tables.table(_reshape_node_yset_es_yt(value.(m[:lsx]) * S_base)), header=false)
    CSV.write(joinpath(output_dir, "load_curtailments.csv"), Tables.table(_reshape_node_yset_es_yt(value.(m[:ls]) * S_base .+ value.(m[:lsx]) * S_base)), header=false)

    ws = value.(m[:ws]) * S_base
    CSV.write(joinpath(output_dir, "wind_spillage.csv"), Tables.table(_reshape_node_yset_es_yt(ws)), header=false)

    bus_out_value = value.(m[:bus_out_power]) * S_base
    CSV.write(joinpath(output_dir, "bus_out_power.csv"), Tables.table(_reshape_node_yset_es_yt(bus_out_value)), header=false)
    if GSw_Battery
        CSV.write(joinpath(output_dir, "ch_node.csv"), Tables.table(_reshape_node_yset_es_yt(value.(m[:ch] ./ ηₛᶜʰ) * S_base)), header=false)
        CSV.write(joinpath(output_dir, "dis_node.csv"), Tables.table(_reshape_node_yset_es_yt(value.(m[:dis] ./ ηₛᵈⁱˢ) * S_base)), header=false)
        CSV.write(joinpath(output_dir, "new_build_gen_cap_bess_power.csv"), Tables.table(value.(m[:Pᵐᵃˣ]) * S_base); header=false)
        CSV.write(joinpath(output_dir, "new_build_gen_cap_bess_energy.csv"), Tables.table(value.(m[:Eᵐᵃˣ]) * S_base); header=false)
        CSV.write(joinpath(output_dir, "annual_cost_bess_investment.csv"), Tables.table(value.(m[:eCBESS_INV])); header=false)
    end
    DC_flow_value = value.(m[:DC_flow]) * S_base
    CSV.write(joinpath(output_dir, "line_flow.csv"), Tables.table(_reshape_line_ylet_el_yt(DC_flow_value)), header=false)

    θ_value = value.(m[:θ])
    CSV.write(joinpath(output_dir, "theta.csv"), Tables.table(_reshape_node_yset_es_yt(θ_value)), header=false)

    new_DC_flow_value = value.(m[:new_DC_flow]) * S_base
    CSV.write(joinpath(output_dir, "upgraded_line_flow.csv"), Tables.table(_reshape_line_ylet_el_yt(new_DC_flow_value)), header=false)

    C_value = value.(m[:eC]) .+ value.(m[:eCG_INV]) .+ value.(m[:eCL_INV]) .+ (GSw_Emission ? value.(m[:eCE_Emission]) : 0) .+ (GSw_AirQuality ? value.(m[:eCE_AirQuality]) : 0) .+ (GSw_Battery ? value.(m[:eCBESS_INV]) : 0)
    CSV.write(joinpath(output_dir, "annual_cost_total.csv"), Tables.table(C_value); header=false)
    CSV.write(joinpath(output_dir, "annual_cost_operation.csv"), Tables.table(value.(m[:eC])); header=false)
    CSV.write(joinpath(output_dir, "annual_cost_gen_investment.csv"), Tables.table(value.(m[:eCG_INV])); header=false)
    CSV.write(joinpath(output_dir, "annual_cost_line_investment.csv"), Tables.table(value.(m[:eCL_INV])); header=false)

    CSV.write(joinpath(output_dir, "annual_cost_operation_over_generation.csv"), Tables.table(value.(m[:eCOverGenPenalty])); header=false)
    CSV.write(joinpath(output_dir, "annual_cost_operation_under_generation.csv"), Tables.table(value.(m[:eCUnderGenPenalty])); header=false)
    CSV.write(joinpath(output_dir, "annual_cost_operation_rps_noncompliance.csv"), Tables.table(value.(m[:eCRPSNCPenalty])); header=false)

    CSV.write(joinpath(output_dir, "new_build_gen_cap_ng_cc_ccs.csv"), Tables.table(value.(m[:p_new_cap_ng_cc_ccs]) * S_base); header=false)
    CSV.write(joinpath(output_dir, "new_build_gen_cap_ng_ct.csv"), Tables.table(value.(m[:p_new_cap_ng_ct]) * S_base); header=false)
    CSV.write(joinpath(output_dir, "new_build_gen_cap_wind.csv"), Tables.table(value.(m[:p_new_cap_wind]) * S_base); header=false)
    CSV.write(joinpath(output_dir, "new_build_gen_cap_pv.csv"), Tables.table(value.(m[:p_new_cap_pv]) * S_base); header=false)

    new_line_decision = round.(Int, value.(m[:Iᴸ]))
    CSV.write(joinpath(output_dir, "line_upgrade_decision.csv"), Tables.table(new_line_decision); header=false)
    CSV.write(joinpath(output_dir, "line_upgrade_decision_size.csv"), Tables.table(value.(m[:Sᴸ]) * S_base); header=false)

    if GSw_DemandFlexibility
        flexible_demand = value.(m[:Δd]) * S_base
        CSV.write(joinpath(output_dir, "flexible_load.csv"), Tables.table(_reshape_node_yset_es_yt(flexible_demand)), header=false)
        CSV.write(joinpath(output_dir, "annual_cost_operation_demand_flexibility.csv"), Tables.table(value.(m[:eCDemandFlexibility])); header=false)
    end

    CSV.write(joinpath(output_dir, "Emissions_CO2.csv"), Tables.table(value.(m[:eE_CO2])); header=false)
    CSV.write(joinpath(output_dir, "Emissions_SO2.csv"), Tables.table(value.(m[:eE_SO2])); header=false)
    CSV.write(joinpath(output_dir, "Emissions_NOx.csv"), Tables.table(value.(m[:eE_NOx])); header=false)
    CSV.write(joinpath(output_dir, "Emissions_PM25.csv"), Tables.table(value.(m[:eE_PM25])); header=false)
    CSV.write(joinpath(output_dir, "annual_cost_externalities_emission.csv"), Tables.table(value.(m[:eCE_Emission])); header=false)
    CSV.write(joinpath(output_dir, "annual_cost_externalities_airquality.csv"), Tables.table(value.(m[:eCE_AirQuality])); header=false)
    CSV.write(joinpath(output_dir, "annual_cost_externalities.csv"), Tables.table(value.(m[:eCE_Emission]) .+ value.(m[:eCE_AirQuality])); header=false)

    if GSw_OSW
        CSV.write(joinpath(output_dir, "new_line_flow.csv"), Tables.table(_reshape_line_ylet_el_yt(value.(m[:new_DC_flow_osw]) * S_base)), header=false)
        CSV.write(joinpath(output_dir, "new_line_decision.csv"), Tables.table(_reshape_line_lcy_l_cy(round.(Int, value.(m[:Iᴸ_osw])))); header=false)
    end
end

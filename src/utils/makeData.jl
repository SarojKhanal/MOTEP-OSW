# Set up environment from .toml files
import Pkg
Pkg.activate(".")
Pkg.instantiate()

# Load necessary packages
using DataFrames, CSV, XLSX, Query

# construct power plant database from EIA 2015 data
# relevant technologies
techn = Dict("Petroleum Liquids" => "Oil",
    "Natural Gas Steam Turbine" => "Gas_Steam",
    "Nuclear" => "Nuclear",
    "Conventional Steam Coal" => "Coal",
    "Natural Gas Fired Combined Cycle" => "Gas_CCGT",
    "Natural Gas Fired Combustion Turbine" => "Gas_GT",
    "Natural Gas Internal Combustion Engine" => "Gas_IC",
    "Other Natural Gas" => "Gas_Other")

isone_states = ["CT", "MA", "ME", "NH", "RI", "VT"] # there are some small plants outside of these states

# Load data
plants = DataFrame(XLSX.readtable("data/eia860/2015/2___Plant_Y2015.xlsx", "Plant", first_row=2))
generators = DataFrame(XLSX.readtable("data/eia860/2015/3_1_Generator_Y2015.xlsx", "Operable", first_row=2))
generators = generators[1:size(generators)[1]-1, :] # skip last entry

plants_isone = plants[∈(isone_states).(plants[:, Symbol("State")]), :]
generators_isone = generators[∈(isone_states).(generators[:, Symbol("State")]), :]


# subset ISONE
plants_isone = DataFrame(
    @from i in plants_isone begin
        @where i[Symbol("Balancing Authority Code")] == "ISNE"
        @select i
        @collect
    end)



generators_isone = generators_isone[in.(generators_isone[:, "Technology"], [Set(keys(techn))]), :]
generators_isone[:, "Type"] = missings(String, nrow(generators_isone))
for j in 1:nrow(generators_isone)
    generators_isone[j, "Type"] = techn[generators_isone[j, "Technology"]]
end

generators_isone[:, "Minimum Load (MW)"] = replace(generators_isone[:, "Minimum Load (MW)"], missing => 0.0)
generators_isone[:, "Minimum Load (MW)"] = replace(generators_isone[:, "Minimum Load (MW)"], " " => 0.0)
gdf = groupby(generators_isone, [:"Utility ID", :"Plant Code", :"State", :"Type", "Utility Name", "Plant Name"])
#combine(gdf, :"County" => sum)
generators_clustered_isone = combine(gdf, :"Nameplate Capacity (MW)" => sum, :"Minimum Load (MW)" => sum)
generators_clustered_isone = innerjoin(generators_clustered_isone, plants_isone[:, ["Plant Code", "Zip", "County", "Latitude", "Longitude"]], on=:"Plant Code")

# add standard operational characteristics
# https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2019/Sep/IRENA_Flexibility_in_CPPs_2019.pdf
# Techology                 Min Load    Avg RampRate    Min DownTime    Min UpTime
# Hard Coal, Average plant, 25–40%      1.5–4%          48h             48h
# CCGT, Average plant       40–50%      2–4%            4h              2h
# OCGT, Average plant       40–50%      2–4%            10-30min        30-60min

# from generator_list.xlsx  Min Load    Avg RampRate        Min DownTime    Min UpTime
# Nuclear                   80%         120                 48h             24h
# Coal                      40%         120                 8-60h           8-24h
# Oil                       35%         400                 
# Gas                       20%         400                 

# average linear cost (at 80% PMAX) and SU by Technology
# Nuclear:      11,     SU: 900,000
# Coal:         25,     SU:  12,586    
# Oil:          207,    SU:  94,146
# Gas:          56,     SU:  14,238
# Gas (Peaker): 343,    SU:   0     


# index, type, Zone, Pmax, Pmin, c2, c1, c0, SUcost, SDcost, RUrate, RDrate, UPtime, DNtime
generators_clustered_isone[:, "index"] = 1:nrow(generators_clustered_isone)
generators_clustered_isone[:, "type"] = generators_clustered_isone[:, "Type"]
for j in 1:nrow(generators_clustered_isone)
    if startswith(generators_clustered_isone[j, "Type"], "Gas_")
        generators_clustered_isone[j, "type"] = "Gas"
    end
end

# pmin and pmax, ramprate
generators_clustered_isone[:, "Pmax"] = generators_clustered_isone[:, "Nameplate Capacity (MW)_sum"]
generators_clustered_isone[:, "Pmin"] = generators_clustered_isone[:, "Minimum Load (MW)_sum"]
generators_clustered_isone[!, "RUrate"] .= 0.0
generators_clustered_isone[!, "RDrate"] .= 0.0

for j in 1:nrow(generators_clustered_isone)
    if generators_clustered_isone[j, "type"] == "Gas"
        generators_clustered_isone[j, "Pmin"] = 0.2 * generators_clustered_isone[j, "Pmax"]
        generators_clustered_isone[j, "RUrate"] = 400.0
        generators_clustered_isone[j, "RDrate"] = 400.0
    elseif generators_clustered_isone[j, "type"] == "Oil"
        generators_clustered_isone[j, "Pmin"] = 0.35 * generators_clustered_isone[j, "Pmax"]
        generators_clustered_isone[j, "RUrate"] = 400.0
        generators_clustered_isone[j, "RDrate"] = 400.0
    elseif generators_clustered_isone[j, "type"] == "Coal"
        generators_clustered_isone[j, "Pmin"] = 0.4 * generators_clustered_isone[j, "Pmax"]
        generators_clustered_isone[j, "RUrate"] = 120.0
        generators_clustered_isone[j, "RDrate"] = 120.0
    elseif generators_clustered_isone[j, "type"] == "Nuclear"
        generators_clustered_isone[j, "Pmin"] = 0.8 * generators_clustered_isone[j, "Pmax"]
        generators_clustered_isone[j, "RUrate"] = 120.0
        generators_clustered_isone[j, "RDrate"] = 120.0
    end
end

# only linear cost, UPtime	DNtime
generators_clustered_isone[!, "c0"] .= 0.0
generators_clustered_isone[!, "c1"] .= 0.0
generators_clustered_isone[!, "c2"] .= 0.0
generators_clustered_isone[!, "SUcost"] .= 0.0
generators_clustered_isone[!, "SDcost"] .= 0.0
generators_clustered_isone[!, "UPtime"] .= 0.0
generators_clustered_isone[!, "DNtime"] .= 0.0
for j in 1:nrow(generators_clustered_isone)
    if generators_clustered_isone[j, "Type"] == "Nuclear"
        generators_clustered_isone[j, "c1"] = 11.0
        generators_clustered_isone[j, "SUcost"] = 900000.0
        generators_clustered_isone[j, "UPtime"] = 24.0
        generators_clustered_isone[j, "DNtime"] = 48.0

    elseif generators_clustered_isone[j, "Type"] == "Coal"
        generators_clustered_isone[j, "c1"] = 25.0
        generators_clustered_isone[j, "SUcost"] = 12586.0
        generators_clustered_isone[j, "UPtime"] = 16.0 # median
        generators_clustered_isone[j, "DNtime"] = 16.0 # median

    elseif generators_clustered_isone[j, "Type"] == "Oil"
        generators_clustered_isone[j, "c1"] = 207.0
        generators_clustered_isone[j, "SUcost"] = 94146.0
        generators_clustered_isone[j, "UPtime"] = 16.0 # median
        generators_clustered_isone[j, "DNtime"] = 11.0 # median

    elseif generators_clustered_isone[j, "Type"] == "Gas_GT"
        generators_clustered_isone[j, "c1"] = 343.0
        generators_clustered_isone[j, "SUcost"] = 0.0
        generators_clustered_isone[j, "UPtime"] = 1.0
        generators_clustered_isone[j, "DNtime"] = 1.0

    else # all other gas including CCGT
        generators_clustered_isone[j, "c1"] = 56.0
        generators_clustered_isone[j, "SUcost"] = 14238.0
        generators_clustered_isone[j, "UPtime"] = 8.0 # median
        generators_clustered_isone[j, "DNtime"] = 6.0 # median
    end
end

# add a small cost adder to ensure that units differ slightly
generators_clustered_isone[:, "c1"] += generators_clustered_isone[:, "Pmax"] / 1000

# West-Center Massachusetts (WCMASS): Berkshire, Franklin, Hampshire, Hampden, Worcester
# Rest: SEMASS
generators_clustered_isone[:, "Zone"] = generators_clustered_isone[:, "State"]
for j in 1:nrow(generators_clustered_isone)
    generators_clustered_isone[j, "County"] = uppercase(generators_clustered_isone[j, "County"])
    if (generators_clustered_isone[j, "State"] == "MA") & (generators_clustered_isone[j, "County"] in Set(["BERKSHIRE", "FRANKLIN", "HAMPSHIRE", "HAMPDEN", "WORCESTER"]))
        generators_clustered_isone[j, "Zone"] = "WCMASS"
    elseif (generators_clustered_isone[j, "State"] == "MA") & (generators_clustered_isone[j, "County"] ∉ Set(["BERKSHIRE", "FRANKLIN", "HAMPSHIRE", "HAMPDEN", "WORCESTER"]))
        generators_clustered_isone[j, "Zone"] = "SEMASS"
    end
end



CSV.write("./data/generators_clustered_isone.csv", generators_clustered_isone)



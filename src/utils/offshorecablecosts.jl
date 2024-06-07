using Plots

# Constants
k_scale_OffshoreLineCost = 0.7340821889389222 # Derived scaling factor

# Function to calculate Cˡ_400
function calculate_Cˡ_400(Lˡ)
    return k_scale_OffshoreLineCost * (0.007225 .* (Lˡ .* 1.609).^2 + 0.767 .* (Lˡ .* 1.609) .+ 32.8125) .* 1e6 .* (1.5285 * 1.09) ./ (0.4 * 1000)
end

# Function to calculate Cˡ_1400
function calculate_Cˡ_1400(Lˡ)
    return k_scale_OffshoreLineCost * (1.36 .* (Lˡ .* 1.609) .+ 366.78) .* 1e6 .* (1.5285 * 1.09) ./ (1.4 * 1000)
end

# Function to calculate Cˡ_2200
function calculate_Cˡ_2200(Lˡ)
    return k_scale_OffshoreLineCost * (1.8 .* (Lˡ .* 1.609) .+ 562.08) .* 1e6 .* (1.5285 * 1.09) ./ (2.2 * 1000)
end

# Generate Lˡ values from 1.38 to 373.8
Lˡ_range = 1.38:0.1:373.8

# Calculate costs for each Lˡ value
Cˡ_400_values = calculate_Cˡ_400.(Lˡ_range)
Cˡ_1400_values = calculate_Cˡ_1400.(Lˡ_range)
Cˡ_2200_values = calculate_Cˡ_2200.(Lˡ_range)

# Plot
plot(Lˡ_range, Cˡ_400_values, label="Cˡ_400", title="Offshore Line Costs", xlabel="Distance (Lˡ) in miles", ylabel="Cost in USD", legend=:topright)
plot!(Lˡ_range, Cˡ_1400_values, label="Cˡ_1400")
plot!(Lˡ_range, Cˡ_2200_values, label="Cˡ_2200")


Lˡ = 37.2823 # Evaluate at reference point
k_scale_OffshoreLineCost = 0.7340821889389221 ## https://guidetoanoffshorewindfarm.com/wind-farm-costs $ 1 2020 USD = $1.01 USD, Average exchange rate in 2019: 1.2772 USD, £/MW: 320000, 37.2823 miles (60 km)
Cˡ_400 = k_scale_OffshoreLineCost * (0.007225 .* (Lˡ .* 1.609).^2 + 0.767 .* (Lˡ .* 1.609) .+ 32.8125) .* 1e6 .* (1.5285 * 1.09) ./(0.4 * 1000 * Lˡ)
Cˡ_1400 = k_scale_OffshoreLineCost * (1.36 .* (Lˡ .* 1.609) .+ 366.78) .* 1e6 .* (1.5285 * 1.09) ./ (1.4 * 1000 * Lˡ)
Cˡ_2200 = k_scale_OffshoreLineCost * (1.8 .* (Lˡ .* 1.609) .+ 562.08) .* 1e6 .* (1.5285 * 1.09) ./ (2.2 * 1000 * Lˡ)
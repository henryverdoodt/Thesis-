## Master Thesis: Analyse the impact of climate change on a G&TEP of EU Power System
# Author: Henry Verdoodt
# Last update: April 4, 2023

## Step 0: Activate environment - ensure consistency accross computers
#= using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
Pkg.instantiate() # If a Manifest.toml file exist in the current project, download all the packages declared in that manifest. Else, resolve a set of feasible packages from the Project.toml files and install them.
Pkg.add("InteractiveUtils")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Distributions")
Pkg.add("Gurobi")
Pkg.add("Images")
Pkg.add("JuMP")
Pkg.add("Plots")
Pkg.add("PolyChaos")
Pkg.add("YAML")
Pkg.add("StatsPlots")
Pkg.add("Images")  =#

#############################################################################################################################################################
#############################################################################################################################################################
###############################                                           ENTSO                                               ###############################
#############################################################################################################################################################
#############################################################################################################################################################



using CSV, YAML, JuMP, DataFrames, Distributions, Gurobi, Images, Plots, PolyChaos, InteractiveUtils, StatsPlots, Images


# Read the CSV file into a DataFrame
demand = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/DEMAND/PECD-country-demand_national_estimates-2025.csv", DataFrame)
solar = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/SOLAR/PECD-2021.3-country-LFSolarPV-2025.csv", DataFrame)
windon = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDON/PECD-2021.3-country-Onshore-2025.csv", DataFrame)
windoff = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDOFF/PECD-2021.3-country-Offshore-2025.csv", DataFrame)

# Get the unique countries in the 'country' column
countries_demand = unique(demand.country) # ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "ME", "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UA", "UK"]
countries_solar= unique(solar.country)    # ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "ME", "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UA", "UK"]
countries_windon = unique(windon.country) # MT is in solar & demand but not in windon => # ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "ME", "MK", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UA", "UK"]
countries_windoff= unique(windoff.country) # ["BE", "DE", "DK", "FI", "FR", "IE", "IT", "NL", "PT", "SE", "UK"]

# filter the data for a specific country and select only the relevant columns
y = 2016.0
countries = ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO"]

dem = DataFrame(col1 = Float64[])
for c in intersect(countries, countries_demand)
    demand1 = dropmissing(demand[demand.country .== c, [:year, :dem_MW]], disallowmissing=true)
    demand2 = dropmissing(demand1[demand1.year .== 2016.0, [:dem_MW]], disallowmissing=true)
    rename!(demand2, :dem_MW => c)
    if isempty(dem)
        dem = hcat(demand2)
    else
        dem = hcat(dem, demand2)
    end
end

sol = DataFrame(col1 = Float64[])
for c in intersect(countries, countries_solar)
    solar1 = dropmissing(solar[solar.country .== c, [:year, :cf]], disallowmissing=true)
    solar2 = dropmissing(solar1[solar1.year .== 2016.0, [:cf]], disallowmissing=true)
    rename!(solar2, :cf => c)
    if isempty(sol)
        sol = hcat(solar2)
    else
        sol = hcat(sol, solar2)
    end
end

won = DataFrame(col1 = Float64[])
for c in intersect(countries, countries_windon)
    windon1 = dropmissing(windon[windon.country .== c, [:year, :cf]], disallowmissing=true)
    windon2 = dropmissing(windon1[windon1.year .== 2016.0, [:cf]], disallowmissing=true)
    rename!(windon2, :cf => c)
    if isempty(won)
        won = hcat(windon2)
    else
        won = hcat(won, windon2)
    end
end


woff = DataFrame(col1 = Float64[])
for c in intersect(countries, countries_windoff)
    windoff1 = dropmissing(windoff[windoff.country .== c, [:year, :cf]], disallowmissing=true)
    windoff2 = dropmissing(windoff1[windoff1.year .== 2016.0, [:cf]], disallowmissing=true)
    rename!(windoff2, :cf => c)
    if isempty(woff)
        woff = hcat(windoff2)
    else
        woff = hcat(woff, windoff2)
    end
end


#############################################################################################################################################################
#############################################################################################################################################################
###############################                                        COPERNICUS                                             ###############################
#############################################################################################################################################################
#############################################################################################################################################################


using CSV, YAML, JuMP, DataFrames, Distributions, Gurobi, Images, Plots, PolyChaos, InteractiveUtils, StatsPlots, Images


# Read the CSV file into a DataFrame
demand = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/DEMAND/PECD-country-demand_national_estimates-2025.csv", DataFrame)
solar = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/SOLAR/PECD-2021.3-country-LFSolarPV-2025.csv", DataFrame)
windon = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDON/PECD-2021.3-country-Onshore-2025.csv", DataFrame)
windoff = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDOFF/PECD-2021.3-country-Offshore-2025.csv", DataFrame)
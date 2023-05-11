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
demand_entso = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/DEMAND/PECD-country-demand_national_estimates-2030.csv", DataFrame)
demand_PT_IT = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/Demand Data/new2_PT_IT.csv", delim=';',header=true, DataFrame)
solar_entso = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/SOLAR/PECD-2021.3-country-LFSolarPV-2030.csv", DataFrame)
windon_entso = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDON/PECD-2021.3-country-Onshore-2030.csv", DataFrame)
windoff_entso = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDOFF/PECD-2021.3-country-Offshore-2030.csv", DataFrame)

new_demand_PT_IT = transform(demand_PT_IT, :PT => ByRow(x -> parse(Float64, replace(x, "," => "."))) => :PT,
                                  :IT => ByRow(x -> parse(Float64, replace(x, "," => "."))) => :IT)


# Available data for following countries in data demand, solar, windon, windoff
countries_demand_entso = unique(demand_entso.country) # ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "ME", "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UA", "UK"]
countries_solar_entso= unique(solar_entso.country)    # ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "ME", "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UA", "UK"]
countries_windon_entso = unique(windon_entso.country) # MT is in solar & demand but not in windon => # ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "ME", "MK", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UA", "UK"]
countries_windoff_entso= unique(windoff_entso.country) # ["BE", "DE", "DK", "FI", "FR", "IE", "IT", "NL", "PT", "SE", "UK"]


# Functions to reformat ENTSO DATA
#= 
function reformat_entso_demand(demand::DataFrame, countries::Vector{String} , countries_demand::Vector{String3}, y::Float64, aggregate_3h::Bool)
    dem = DataFrame(col1 = Float64[])
    for c in intersect(countries, countries_demand)
        demand1 = dropmissing(demand[demand.country .== c, [:year, :dem_MW]], disallowmissing=true)
        demand2 = dropmissing(demand1[demand1.year .== y, [:dem_MW]], disallowmissing=true)
        rename!(demand2, :dem_MW => c)
        dem = isempty(dem) ? hcat(demand2) : hcat(dem, demand2)
    end
    if aggregate_3h
        nrows_new = Int(ceil(size(dem, 1)/3))
        dem_3h = DataFrame()
        for col in names(dem)
            col_data = []
            for i in 1:nrows_new
                row_start = (i-1)*3 + 1
                row_end = min(i*3, size(dem, 1))
                push!(col_data, mean(dem[row_start:row_end, col]))
            end
            dem_3h[!, col] = col_data
        end
        return dem_3h
    else
        return dem
    end
end
 =#

function reformat_entso_demand(data::DataFrame, countries::Vector{String}, countries_dem::Vector{String3}, years::Vector{Float64}, weights::Vector{Float64}, aggregate_3h::Bool)
    d = DataFrame()
    for c in intersect(countries, countries_dem)
        dem_values = []
        for i in 1:length(years)
            year = years[i]
            weight = weights[i]
            data1 = dropmissing(data[data.country .== c, [:year, :dem_MW]], disallowmissing=true)
            data2 = dropmissing(data1[data1.year .== year, [:dem_MW]], disallowmissing=true)
            rename!(data2, :dem_MW => c)
            push!(dem_values, data2[!, c] * weight)
        end
        sum_vectors(x, y) = x .+ y
        sum_of_vectors = reduce(sum_vectors, dem_values)
        sov = DataFrame(c => sum_of_vectors)
        d = isempty(d) ? hcat(sov) : hcat(d, sov)
    end
    if aggregate_3h
        nrows_new = Int(ceil(size(d, 1)/3))
        dem_3h = DataFrame()
        for col in names(d)
            col_data = []
            for i in 1:nrows_new
                row_start = (i-1)*3 + 1
                row_end = min(i*3, size(d, 1))
                push!(col_data, mean(d[row_start:row_end, col]))
            end
            dem_3h[!, col] = col_data
        end
        return dem_3h
    else
        return d   
    end
end

function reformat_demand_PT_IT(data::DataFrame, aggregate_3h::Bool)
    d = data
    if aggregate_3h
        nrows_new = Int(ceil(size(d, 1)/3))
        dem_3h = DataFrame()
        for col in names(d)
            col_data = []
            for i in 1:nrows_new
                row_start = (i-1)*3 + 1
                row_end = min(i*3, size(d, 1))
                push!(col_data, mean(d[row_start:row_end, col]))
            end
            dem_3h[!, col] = col_data
        end
        return dem_3h
    else
        return d   
    end
end

#= 
function reformat_entso_solar(data::DataFrame, countries::Vector{String} , countries_cf::Vector{String3}, y::Float64, aggregate_3h::Bool)
    d = DataFrame(col1 = Float64[])
    for c in intersect(countries, countries_cf)
        data1 = dropmissing(data[data.country .== c, [:year, :cf]], disallowmissing=true)
        data2 = dropmissing(data1[data1.year .== y, [:cf]], disallowmissing=true)
        rename!(data2, :cf => c)
        d = isempty(d) ? hcat(data2) : hcat(d, data2)
    end
    if aggregate_3h
        nrows_new = Int(ceil(size(d, 1)/3))
        sol_3h = DataFrame()
        for col in names(d)
            col_data = []
            for i in 1:nrows_new
                row_start = (i-1)*3 + 1
                row_end = min(i*3, size(d, 1))
                push!(col_data, mean(d[row_start:row_end, col]))
            end
            sol_3h[!, col] = col_data
        end
        return sol_3h
    else
        return d   
    end
end
 =#
function reformat_entso_solar(data::DataFrame, countries::Vector{String}, countries_cf::Vector{String3}, years::Vector{Float64}, weights::Vector{Float64}, aggregate_3h::Bool)
    d = DataFrame()
    for c in intersect(countries, countries_cf)
        cf_values = []
        for i in 1:length(years)
            year = years[i]
            weight = weights[i]
            data1 = dropmissing(data[data.country .== c, [:year, :cf]], disallowmissing=true)
            data2 = dropmissing(data1[data1.year .== year, [:cf]], disallowmissing=true)
            rename!(data2, :cf => c)
            push!(cf_values, data2[!, c] * weight)
        end
        sum_vectors(x, y) = x .+ y
        sum_of_vectors = reduce(sum_vectors, cf_values)
        sov = DataFrame(c => sum_of_vectors)
        d = isempty(d) ? hcat(sov) : hcat(d, sov)
    end
    if aggregate_3h
        nrows_new = Int(ceil(size(d, 1)/3))
        sol_3h = DataFrame()
        for col in names(d)
            col_data = []
            for i in 1:nrows_new
                row_start = (i-1)*3 + 1
                row_end = min(i*3, size(d, 1))
                push!(col_data, mean(d[row_start:row_end, col]))
            end
            sol_3h[!, col] = col_data
        end
        return sol_3h
    else
        return d   
    end
end

#countries = ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO", "CH", "FI", "IE", "IT", "AT", "PT", "SE"]   #["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO"]
#reformat_entso_solar(solar_entso, countries, countries_solar_entso, [1995.0, 2008.0, 2009.0], [0.233, 0.367, 0.4], false)

#= 
function reformat_entso_windon(windon::DataFrame, countries::Vector{String} , countries_windon::Vector{String3}, y::Float64, aggregate_3h::Bool)
    won = DataFrame(col1 = Float64[])
    for c in intersect(countries, countries_windon)
        windon1 = dropmissing(windon[windon.country .== c, [:year, :cf]], disallowmissing=true)
        windon2 = dropmissing(windon1[windon1.year .== y, [:cf]], disallowmissing=true)
        rename!(windon2, :cf => c)
        won = isempty(won) ? hcat(windon2) : hcat(won, windon2)
    end
    if aggregate_3h
        nrows_new = Int(ceil(size(won, 1)/3))
        won_3h = DataFrame()
        for col in names(won)
            col_data = []
            for i in 1:nrows_new
                row_start = (i-1)*3 + 1
                row_end = min(i*3, size(won, 1))
                push!(col_data, mean(won[row_start:row_end, col]))
            end
            won_3h[!, col] = col_data
        end
        return won_3h
    else
        return won
    end
end
 =#

function reformat_entso_windon(data::DataFrame, countries::Vector{String}, countries_windon::Vector{String3}, years::Vector{Float64}, weights::Vector{Float64}, aggregate_3h::Bool)
    d = DataFrame()
    for c in intersect(countries, countries_windon)
        cf_values = []
        for i in 1:length(years)
            year = years[i]
            weight = weights[i]
            data1 = dropmissing(data[data.country .== c, [:year, :cf]], disallowmissing=true)
            data2 = dropmissing(data1[data1.year .== year, [:cf]], disallowmissing=true)
            rename!(data2, :cf => c)
            push!(cf_values, data2[!, c] * weight)
        end
        sum_vectors(x, y) = x .+ y
        sum_of_vectors = reduce(sum_vectors, cf_values)
        sov = DataFrame(c => sum_of_vectors)
        d = isempty(d) ? hcat(sov) : hcat(d, sov)
    end
    if aggregate_3h
        nrows_new = Int(ceil(size(d, 1)/3))
        won_3h = DataFrame()
        for col in names(d)
            col_data = []
            for i in 1:nrows_new
                row_start = (i-1)*3 + 1
                row_end = min(i*3, size(d, 1))
                push!(col_data, mean(d[row_start:row_end, col]))
            end
            won_3h[!, col] = col_data
        end
        return won_3h
    else
        return d   
    end
end
#= 
function reformat_entso_windoff(windoff::DataFrame, countries::Vector{String} , countries_windoff::Vector{String3}, y::Float64, aggregate_3h::Bool)
    woff = DataFrame(col1 = Float64[])
    for c in intersect(countries, countries_windoff)
        windoff1 = dropmissing(windoff[windoff.country .== c, [:year, :cf]], disallowmissing=true)
        windoff2 = dropmissing(windoff1[windoff1.year .== y, [:cf]], disallowmissing=true)
        rename!(windoff2, :cf => c)
        woff = isempty(woff) ? hcat(windoff2) : hcat(woff, windoff2)
    end
    if aggregate_3h        
        nrows_new = Int(ceil(size(woff, 1)/3))
        woff_3h = DataFrame()
        for col in names(woff)
            col_data = []
            for i in 1:nrows_new
                row_start = (i-1)*3 + 1
                row_end = min(i*3, size(woff, 1))
                push!(col_data, mean(woff[row_start:row_end, col]))
            end
            woff_3h[!, col] = col_data
        end
        return woff_3h
    else
        return woff
    end
end
 =#
function reformat_entso_windoff(data::DataFrame, countries::Vector{String}, countries_windoff::Vector{String3}, years::Vector{Float64}, weights::Vector{Float64}, aggregate_3h::Bool)
    d = DataFrame()
    for c in intersect(countries, countries_windoff)
        cf_values = []
        for i in 1:length(years)
            year = years[i]
            weight = weights[i]
            data1 = dropmissing(data[data.country .== c, [:year, :cf]], disallowmissing=true)
            data2 = dropmissing(data1[data1.year .== year, [:cf]], disallowmissing=true)
            rename!(data2, :cf => c)
            push!(cf_values, data2[!, c] * weight)
        end
        sum_vectors(x, y) = x .+ y
        sum_of_vectors = reduce(sum_vectors, cf_values)
        sov = DataFrame(c => sum_of_vectors)
        d = isempty(d) ? hcat(sov) : hcat(d, sov)
    end
    if aggregate_3h
        nrows_new = Int(ceil(size(d, 1)/3))
        woff_3h = DataFrame()
        for col in names(d)
            col_data = []
            for i in 1:nrows_new
                row_start = (i-1)*3 + 1
                row_end = min(i*3, size(d, 1))
                push!(col_data, mean(d[row_start:row_end, col]))
            end
            woff_3h[!, col] = col_data
        end
        return woff_3h
    else
        return d   
    end
end


# Parameter for function to reformat data
#countries = ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO", "CH", "FI", "IE", "IT", "AT", "PT", "SE"]   #["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO"]
#y = 1982.0    # demand_entso: 1982.0 - 2016.0    and   solar,windon,windoff_entso: 1982.0 - 2019.0
#aggregate_3h = true

#= 
# Call reformated ENTSO DATA
reformat_entso_demand(demand_entso, countries, countries_demand_entso, y, aggregate_3h)
reformat_entso_solar(solar_entso, countries, countries_solar_entso, y, aggregate_3h)
reformat_entso_windon(windon_entso, countries, countries_windon_entso, y, aggregate_3h)
reformat_entso_windoff(windoff_entso, countries, countries_windoff_entso, y, aggregate_3h)

 =#

#############################################################################################################################################################
#############################################################################################################################################################
###############################                                        COPERNICUS                                             ###############################
#############################################################################################################################################################
#############################################################################################################################################################

using CSV, YAML, JuMP, DataFrames, Distributions, Gurobi, Images, Plots, PolyChaos, InteractiveUtils, StatsPlots, Images, Dates

# Read the CSV file into a DataFrame
solar_coper_cnrm = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/SOLAR/SOLAR_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)
windon_coper_cnrm = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDON/WINDON_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)
windoff_coper_cnrm = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDOFF/WINDOFF_MARITIME_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)

# Available data for following countries in solar, windon, windoff
countries_solar_coper_cnrm = setdiff(names(solar_coper_cnrm), ["Date"])
countries_windon_coper_cnrm = setdiff(names(windon_coper_cnrm), ["Date"])
countries_windoff_coper_cnrm = setdiff(names(windoff_coper_cnrm), ["Date"])

# Functions to reformat COPERNICUS DATA
function reformat_coper_solar(solar::DataFrame, countries::Vector{String} , countries_solar::Vector{String}, year::Int64)
    df_countries = select(solar, [:Date, Symbol.(intersect(countries, countries_solar))]...)
    df_year = df_countries[findall(x -> year == parse(Int, split(x, "-")[1]), df_countries.Date), :]
    sol = dropmissing(select(df_year, Not(:Date)), disallowmissing=true)

    replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v)
    df2 = map(replace_nan, eachcol(sol))

    sol2 = DataFrame()
    for (i, col) in enumerate(names(sol))
        sol2[!, col] = df2[i]
    end
    return sol2
end


function reformat_coper_windon(windon::DataFrame, countries::Vector{String} , countries_windon::Vector{String}, year::Int64)
    df_countries = select(windon, [:Date, Symbol.(intersect(countries, countries_windon))]...)
    df_year = df_countries[findall(x -> year == parse(Int, split(x, "-")[1]), df_countries.Date), :]
    won = dropmissing(select(df_year, Not(:Date)), disallowmissing=true)

    replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v)
    df2 = map(replace_nan, eachcol(won))

    won2 = DataFrame()
    for (i, col) in enumerate(names(won))
        won2[!, col] = df2[i]
    end
    return won2
end

function reformat_coper_windoff(windoff::DataFrame, countries::Vector{String} , countries_windoff::Vector{String}, year::Int64)
    df_countries = select(windoff, [:Date, Symbol.(intersect(countries, countries_windoff))]...)
    df_year = df_countries[findall(x -> year == parse(Int, split(x, "-")[1]), df_countries.Date), :]
    woff = dropmissing(select(df_year, Not(:Date)), disallowmissing=true)

    replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v)
    df2 = map(replace_nan, eachcol(woff))

    woff2 = DataFrame()
    for (i, col) in enumerate(names(woff))
        woff2[!, col] = df2[i]
    end
    return woff2
end


# Parameter for function to reformat data
#year = 2030    # 1952-2100
#countries = ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO", "CH", "FI", "IE", "IT", "AT", "PT", "SE"]     #["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO"]

#= 
# Call reformated COPERNICUS DATA
reformat_coper_solar(solar_coper_cnrm, countries, countries_solar_coper_cnrm, year)
reformat_coper_windon(windon_coper_cnrm, countries, countries_windon_coper_cnrm, year)
reformat_coper_windoff(windoff_coper_cnrm, countries, countries_windoff_coper_cnrm, year)
 =#

#############################################################################################################################################################
#############################################################################################################################################################


using CSV, FileIO

# Read the CSV file
df = CSV.File("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/Demand Data/PT_IT.csv") |> DataFrame

# Replace commas with dots
for col in names(df)
    if eltype(df[!, col]) == String
        df[!, col] = replace.(df[!, col], "," => ".")
    end
end

# Write the modified DataFrame back to a CSV file
CSV.write("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/Demand Data/new4_PT_IT.csv", df)
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
solar_coper_CNRM = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/SOLAR/SOLAR_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)
windon_coper_CNRM = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDON/WINDON_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)
windoff_coper_CNRM = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDOFF/WINDOFF_MARITIME_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)

solar_coper_EARTH = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/SOLAR/SOLAR_CL_CF_3H_RCP8_HIRHAM_EARTH_1951_2100.csv", DataFrame)
windon_coper_EARTH = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDON/WINDON_CL_CF_3H_RCP8_HIRHAM_EARTH_1951_2100.csv", DataFrame)
windoff_coper_EARTH = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDOFF/WINDOFF_MARITIME_CL_CF_3H_RCP8_HIRHAM_EARTH_1951_2100.csv", DataFrame)

solar_coper_HadGEM = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/SOLAR/SOLAR_CL_CF_3H_RCP8_RegCM_HadGEM_1971_2099.csv", DataFrame)
windon_coper_HadGEM = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDON/WINDON_CL_CF_3H_RCP8_RegCM_HadGEM_1971_2099.csv", DataFrame)
windoff_coper_HadGEM = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDOFF/WINDOFF_MARITIME_CL_CF_3H_RCP8_RegCM_HadGEM_1971_2099.csv", DataFrame)

# Available data for following countries in solar, windon, windoff
countries_solar_coper_CNRM = setdiff(names(solar_coper_CNRM), ["Date"])
countries_windon_coper_CNRM = setdiff(names(windon_coper_CNRM), ["Date"])
countries_windoff_coper_CNRM = setdiff(names(windoff_coper_CNRM), ["Date"])

countries_solar_coper_EARTH = setdiff(names(solar_coper_EARTH), ["Date"])
countries_windon_coper_EARTH = setdiff(names(windon_coper_EARTH), ["Date"])
countries_windoff_coper_EARTH = setdiff(names(windoff_coper_EARTH), ["Date"])

countries_solar_coper_HadGEM = setdiff(names(solar_coper_HadGEM), ["Date"])
countries_windon_coper_HadGEM = setdiff(names(windon_coper_HadGEM), ["Date"])
countries_windoff_coper_HadGEM = setdiff(names(windoff_coper_HadGEM), ["Date"])

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

function reformat_coper_solar_multiple_years(solar::DataFrame, countries::Vector{String}, countries_solar::Vector{String}, start_year::Int64, end_year::Int64)
    df_countries = select(solar, [:Date, Symbol.(intersect(countries, countries_solar))]...)
    df_years = df_countries[findall(x -> start_year <= parse(Int, split(x, "-")[1]) <= end_year, df_countries.Date), :]
    sol = dropmissing(select(df_years, Not(:Date)), disallowmissing=true)

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

function reformat_coper_windon_multiple_years(windon::DataFrame, countries::Vector{String}, countries_windon::Vector{String}, start_year::Int64, end_year::Int64)
    df_countries = select(windon, [:Date, Symbol.(intersect(countries, countries_windon))]...)
    df_years = df_countries[findall(x -> start_year <= parse(Int, split(x, "-")[1]) <= end_year, df_countries.Date), :]
    won = dropmissing(select(df_years, Not(:Date)), disallowmissing=true)

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

function reformat_coper_windoff_multiple_years(windoff::DataFrame, countries::Vector{String}, countries_windoff::Vector{String}, start_year::Int64, end_year::Int64)
    df_countries = select(windoff, [:Date, Symbol.(intersect(countries, countries_windoff))]...)
    df_years = df_countries[findall(x -> start_year <= parse(Int, split(x, "-")[1]) <= end_year, df_countries.Date), :]
    woff = dropmissing(select(df_years, Not(:Date)), disallowmissing=true)

    replace_nan(v) = map(x -> isnan(x) ? zero(x) : x, v)
    df2 = map(replace_nan, eachcol(woff))

    woff2 = DataFrame()
    for (i, col) in enumerate(names(woff))
        woff2[!, col] = df2[i]
    end
    return woff2
end


using DataFrames
using StatsPlots
using Statistics
using DataStructures
using Plots
using PrettyTables

countries = ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO", "CH", "FI", "IE", "IT", "AT", "PT", "SE"] 
wind_df1 = reformat_coper_windon_multiple_years(windon_coper_CNRM, countries, countries_windon_coper_CNRM, 1986, 2015)
solar_df1 = reformat_coper_solar_multiple_years(solar_coper_CNRM, countries, countries_solar_coper_CNRM, 1986, 2015)
wind_df2 = reformat_coper_windon_multiple_years(windon_coper_CNRM, countries, countries_windon_coper_CNRM, 2070, 2099)
solar_df2 = reformat_coper_solar_multiple_years(solar_coper_CNRM, countries, countries_solar_coper_CNRM, 2070, 2099)
wind_df3 = reformat_coper_windon_multiple_years(windon_coper_EARTH, countries, countries_windon_coper_EARTH, 2070, 2099)
solar_df3 = reformat_coper_solar_multiple_years(solar_coper_EARTH, countries, countries_solar_coper_EARTH, 2070, 2099)
wind_df4 = reformat_coper_windon_multiple_years(windon_coper_HadGEM, countries, countries_windon_coper_HadGEM, 2070, 2099)
solar_df4 = reformat_coper_solar_multiple_years(solar_coper_HadGEM, countries, countries_solar_coper_HadGEM, 2070, 2099)

North = [:NO, :SE, :FI]
South = [:PT, :ES, :IT]

function countmap(data)
    counts = Dict{eltype(data), Int}()
    for value in data
        counts[value] = get(counts, value, 0) + 1
    end
    return counts
end

function cumul_dunkelflaute_periods(original_dict::OrderedDict)
    modified_dict = OrderedDict()
    keys_list = collect(keys(original_dict))
    count = 2 
    for (key, value) in original_dict
        if count <= length(keys_list)
            modified_value = value + sum(original_dict[keys_list[j]] for j in (count):length(keys_list))
            modified_dict[key] = modified_value
            count += 1
        end
    end
    return modified_dict
end

function analyze_and_plot_dunkelflaute(solar_df1::DataFrame, wind_df1::DataFrame, solar_df2::DataFrame, wind_df2::DataFrame, solar_df3::DataFrame, wind_df3::DataFrame, solar_df4::DataFrame, wind_df4::DataFrame, cf_threshold::Float64, duration_threshold::Int64 , countries::Vector{Symbol}, numerical::Bool)
    # Identify dunkelflaute events
    dunkelflaute1_1 = (solar_df1[!, countries[1]] .< cf_threshold) .& (wind_df1[!, countries[1]] .< cf_threshold)
    dunkelflaute1_2 = (solar_df1[!, countries[2]] .< cf_threshold) .& (wind_df1[!, countries[2]] .< cf_threshold)
    dunkelflaute1_3 = (solar_df1[!, countries[3]] .< cf_threshold) .& (wind_df1[!, countries[3]] .< cf_threshold)
    dunkelflaute1 = dunkelflaute1_1 .& dunkelflaute1_2 .& dunkelflaute1_3

    dunkelflaute2_1 = (solar_df2[!, countries[1]] .< cf_threshold) .& (wind_df2[!, countries[1]] .< cf_threshold)
    dunkelflaute2_2 = (solar_df2[!, countries[2]] .< cf_threshold) .& (wind_df2[!, countries[2]] .< cf_threshold)
    dunkelflaute2_3 = (solar_df2[!, countries[3]] .< cf_threshold) .& (wind_df2[!, countries[3]] .< cf_threshold)
    dunkelflaute2 = dunkelflaute2_1 .& dunkelflaute2_2 .& dunkelflaute2_3

    dunkelflaute3_1 = (solar_df3[!, countries[1]] .< cf_threshold) .& (wind_df3[!, countries[1]] .< cf_threshold)
    dunkelflaute3_2 = (solar_df3[!, countries[2]] .< cf_threshold) .& (wind_df3[!, countries[2]] .< cf_threshold)
    dunkelflaute3_3 = (solar_df3[!, countries[3]] .< cf_threshold) .& (wind_df3[!, countries[3]] .< cf_threshold)
    dunkelflaute3 = dunkelflaute3_1 .& dunkelflaute3_2 .& dunkelflaute3_3

    dunkelflaute4_1 = (solar_df4[!, countries[1]] .< cf_threshold) .& (wind_df4[!, countries[1]] .< cf_threshold)
    dunkelflaute4_2 = (solar_df4[!, countries[2]] .< cf_threshold) .& (wind_df4[!, countries[2]] .< cf_threshold)
    dunkelflaute4_3 = (solar_df4[!, countries[3]] .< cf_threshold) .& (wind_df4[!, countries[3]] .< cf_threshold)
    dunkelflaute4 = dunkelflaute4_1 .& dunkelflaute4_2 .& dunkelflaute4_3
    
    #return dunkelflaute
    #dunkelflaute = (solar_df1[!, countries[1]] .< cf_threshold) .& (wind_df1[!, countries[1]] .< cf_threshold)

    # Initialize variables
    event_lengths1 = Int[]
    count1 = 0
    for value in dunkelflaute1  #column_data
        if value
            count1 += 1
        else
            if count1 >= duration_threshold
                push!(event_lengths1, count1)
            end
            count1 = 0
        end
    end
    # Initialize variables
    event_lengths2 = Int[]
    count2 = 0
    for value in dunkelflaute2  #column_data
        if value
            count2 += 1
        else
            if count2 >= duration_threshold
                push!(event_lengths2, count2)
            end
            count2 = 0
        end
    end
    # Initialize variables
    event_lengths3 = Int[]
    count3 = 0
    for value in dunkelflaute3  #column_data
        if value
            count3 += 1
        else
            if count3 >= duration_threshold
                push!(event_lengths3, count3)
            end
            count3 = 0
        end
    end
    # Initialize variables
    event_lengths4 = Int[]
    count4 = 0
    for value in dunkelflaute4  #column_data
        if value
            count4 += 1
        else
            if count4 >= duration_threshold
                push!(event_lengths4, count4)
            end
            count4 = 0
        end
    end
    # Calculate event frequencies
    counts1 = Dict{Int, Int}()
    for value in event_lengths1
        counts1[value] = get(counts1, value, 0) + 1
    end
    # Calculate event frequencies
    counts2 = Dict{Int, Int}()
    for value in event_lengths2
        counts2[value] = get(counts2, value, 0) + 1
    end
    # Calculate event frequencies
    counts3 = Dict{Int, Int}()
    for value in event_lengths3
        counts3[value] = get(counts3, value, 0) + 1
    end
    # Calculate event frequencies
    counts4 = Dict{Int, Int}()
    for value in event_lengths4
        counts4[value] = get(counts4, value, 0) + 1
    end
    event_frequencies1_cumul = cumul_dunkelflaute_periods(OrderedDict(sort(collect(counts1))))
    event_frequencies2_cumul = cumul_dunkelflaute_periods(OrderedDict(sort(collect(counts2))))
    event_frequencies3_cumul = cumul_dunkelflaute_periods(OrderedDict(sort(collect(counts3))))
    event_frequencies4_cumul = cumul_dunkelflaute_periods(OrderedDict(sort(collect(counts4))))

    event_frequencies1 = OrderedDict(sort(collect(counts1)))
    event_frequencies2 = OrderedDict(sort(collect(counts2)))
    event_frequencies3 = OrderedDict(sort(collect(counts3)))
    event_frequencies4 = OrderedDict(sort(collect(counts4)))

    if numerical
        result1 = sum([k * v for (k, v) in zip(keys(event_frequencies1) ./ 8, values(event_frequencies1))])
        result2 = sum([k * v for (k, v) in zip(keys(event_frequencies2) ./ 8, values(event_frequencies2))])
        result3 = sum([k * v for (k, v) in zip(keys(event_frequencies3) ./ 8, values(event_frequencies3))])
        result4 = sum([k * v for (k, v) in zip(keys(event_frequencies4) ./ 8, values(event_frequencies4))])
        # Generate some example data
        results = OrderedDict("Historical" => result1, "CNRM" => result2, "EARTH" => result3, "HadGEM" => result4)
        return results
    end
    # Generate the plot
    plot(collect(keys(event_frequencies1_cumul)) ./ 8, collect(values(event_frequencies1_cumul)), xlabel = "Duration (days)", ylabel = "Frequency", legend = true, label = "Historical (1986-2015)", title = "DunkelFlaute - North Europe (CF threshold: $(cf_threshold))")
    plot!(collect(keys(event_frequencies2_cumul)) ./ 8, collect(values(event_frequencies2_cumul)), xlabel = "Duration (days)", ylabel = "Frequency", legend = true, label = "CNRM (2070-2099)")
    plot!(collect(keys(event_frequencies3_cumul)) ./ 8, collect(values(event_frequencies3_cumul)), xlabel = "Duration (days)", ylabel = "Frequency", legend = true, label = "EARTH (2070-2099)")
    plot!(collect(keys(event_frequencies4_cumul)) ./ 8, collect(values(event_frequencies4_cumul)), xlabel = "Duration (days)", ylabel = "Frequency", legend = true, label = "HadGEM (2070-2099)") 

    #bar(collect(keys(event_frequencies)), collect(values(event_frequencies)), xlabel = "Number of Consecutive 'true'", ylabel = "Frequency", legend = false, title = "Consecutive 'true' Occurrences")
    #return plot1
end

analyze_and_plot_dunkelflaute(solar_df1, wind_df1, solar_df2, wind_df2, solar_df3, wind_df3, solar_df4, wind_df4, 0.2, 8, North, false)


results = OrderedDict("Historical" => Float64[], "CNRM" => Float64[], "EARTH" => Float64[], "HadGEM" => Float64[]) 
for cf in 0:0.02:0.98
    dict = analyze_and_plot_dunkelflaute(solar_df1, wind_df1, solar_df2, wind_df2, solar_df3, wind_df3, solar_df4, wind_df4, cf, 8, North, true)
    push!(results["Historical"], dict["Historical"])
    push!(results["CNRM"], dict["CNRM"])
    push!(results["EARTH"], dict["EARTH"])
    push!(results["HadGEM"], dict["HadGEM"])
end
x_values = collect(1:length(results["Historical"])) ./ length(results["Historical"])
plot(x_values, results["Historical"], xlabel = "CF Threshold", ylabel = "Persistence time (days)", legend = true, label = "Historical (1986-2015)", title = "DunkelFlaute - North Europe")
plot!(x_values, results["CNRM"], xlabel = "CF Threshold", ylabel = "Persistence time (days)", legend = true, label = "CNRM (2070-2099)")
plot!(x_values, results["EARTH"], xlabel = "CF Threshold", ylabel = "Persistence time (days)", legend = true, label = "EARTH (2070-2099)")
plot!(x_values, results["HadGEM"], xlabel = "CF Threshold", ylabel = "Persistence time (days)", legend = true, label = "HadGEM (2070-2099)") 








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
## Master Thesis: Analyse the impact of climate change on a G&TEP of EU Power System
# Author: Henry Verdoodt
# Last update: April 4, 2023

#############################################################################################################################################################
#############################################################################################################################################################
###############################                                        COPERNICUS                                             ###############################
#############################################################################################################################################################
#############################################################################################################################################################

using CSV, DataFrames, Plots, Statistics

# Load data from csv file
solar_CNRM = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/SOLAR/SOLAR_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)
windon_CNRM = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDON/WINDON_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)
windoff_CNRM = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDOFF/WINDOFF_MARITIME_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)

solar_EARTH = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/SOLAR/SOLAR_CL_CF_3H_RCP8_HIRHAM_EARTH_1951_2100.csv", DataFrame)
windon_EARTH = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDON/WINDON_CL_CF_3H_RCP8_HIRHAM_EARTH_1951_2100.csv", DataFrame)
windoff_EARTH = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDOFF/WINDOFF_MARITIME_CL_CF_3H_RCP8_HIRHAM_EARTH_1951_2100.csv", DataFrame)

solar_HadGEM = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/SOLAR/SOLAR_CL_CF_3H_RCP8_RegCM_HadGEM_1971_2100.csv", DataFrame)
windon_HadGEM = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDON/WINDON_CL_CF_3H_RCP8_RegCM_HadGEM_1971_2100.csv", DataFrame)
windoff_HadGEM = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDOFF/WINDOFF_MARITIME_CL_CF_3H_RCP8_RegCM_HadGEM_1971_2100.csv", DataFrame)

CF_CNRM = [solar_CNRM, windon_CNRM, windoff_CNRM]
CF_CNRM_name = ["solar_CNRM", "windon_CNRM", "windoff_CNRM"]

CF_EARTH = [solar_EARTH, windon_EARTH, windoff_EARTH]
CF_EARTH_name = ["solar_EARTH", "windon_EARTH", "windoff_EARTH"]

CF_HadGEM = [solar_HadGEM, windon_HadGEM, windoff_HadGEM]
CF_HadGEM_name = ["solar_HadGEM", "windon_HadGEM", "windoff_HadGEM"]

summer = ["06", "07", "08"]; autumn = ["09", "10", "11"]; winter = ["12", "01", "02"]; spring = ["03", "04", "05"]
seasons = [summer, autumn, winter, spring]
seasons_name = ["summer", "autumn", "winter", "spring"]

North_EU = ["NO", "SE", "FI"]; South_EU = ["PT", "ES", "IT"]
region = [North_EU, South_EU]
region_name = ["North_EU", "South_EU"]



#############################################################################################################################################################


function csv_boxplot_CF_Coper(data::DataFrame, data_name::String, BeginYear::Int64, EndYear::Int64, Season::Vector{String}, season_name::String, Yearly::Bool, region::Vector{String},region_name::String)
    rows1 = findall(row -> (BeginYear <= parse(Int, split(row, "-")[1]) <= EndYear), data[!, :Date])  
    data_year = data[rows1, :]
    if Yearly
        data_year_region = filter(!isnan, vcat(data_year[!, Symbol(region[1])], data_year[!, Symbol(region[2])], data_year[!, Symbol(region[3])]))
        output_file = "/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/$(data_name)_$(BeginYear)_$(EndYear)_year_$(region_name).csv"
        CSV.write(output_file, DataFrame(data=data_year_region))  
        #return data_year_region
    else
        data_season = data_year[findall(row -> split(row, "-")[2] in Set(Season), data_year[!, :Date]), :]
        data_season_region = filter(!isnan, vcat(data_season[!, Symbol(region[1])], data_season[!, Symbol(region[2])], data_season[!, Symbol(region[3])]))
        output_file = "/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/$(data_name)_$(BeginYear)_$(EndYear)_$(season_name)_$(region_name).csv"   
        CSV.write(output_file, DataFrame(data=data_season_region))
        #return data_season_region
    end
end

function loop_csv_boxplot_CF_Coper(CF_CM::Vector{DataFrame}, CF_CM_name::Vector{String},seasons::Vector{Vector{String}}, seasons_name::Vector{String}, region::Vector{Vector{String}}, region_name::Vector{String})
    count_cf = 0
    for cf in CF_CM
        count_cf +=1
        #println(CF_CM_name[count_cf])
        count_r = 0
        for r in region
            count_r +=1
            #println(region_name[count_r])
            count_s =0
            for s in seasons
                count_s +=1
                #println(seasons_name[count_s])
                csv_boxplot_CF_Coper(cf, CF_CM_name[count_cf], 2071, 2100, s, seasons_name[count_s], false, r, region_name[count_r])
            end
            csv_boxplot_CF_Coper(cf, CF_CM_name[count_cf], 2071, 2100, summer, "summer", true, r, region_name[count_r])
        end
    end
end

#csv_boxplot_CF_Coper(solar_CNRM, "solar_CNRM", 2071, 2100, summer, "summer", false, North_EU, "South_EU")
loop_csv_boxplot_CF_Coper(CF_CM, CF_CM_name, seasons, seasons_name, region, region_name)







#############################################################################################################################################################
#############################################################################################################################################################
###############################                                             ENTSO                                             ###############################
#############################################################################################################################################################
#############################################################################################################################################################

using CSV, DataFrames, Plots, Statistics, Dates


# Read the CSV file into a DataFrame
solar_ENTSO_25 = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/SOLAR/PECD-2021.3-country-LFSolarPV-2025.csv", DataFrame)
windon_ENTSO_25 = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDON/PECD-2021.3-country-Onshore-2025.csv", DataFrame)
windoff_ENTSO_25 = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDOFF/PECD-2021.3-country-Offshore-2025.csv", DataFrame)

solar_ENTSO_30 = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/SOLAR/PECD-2021.3-country-LFSolarPV-2030.csv", DataFrame)
windon_ENTSO_30 = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDON/PECD-2021.3-country-Onshore-2030.csv", DataFrame)
windoff_ENTSO_30 = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDOFF/PECD-2021.3-country-Offshore-2030.csv", DataFrame)

CF_ENTSO_25 = [solar_ENTSO_25, windon_ENTSO_25, windoff_ENTSO_25]
CF_ENTSO_25_name = ["solar_ENTSO_25", "windon_ENTSO_25", "windoff_ENTSO_25"]

CF_ENTSO_30 = [solar_ENTSO_30, windon_ENTSO_30, windoff_ENTSO_30]
CF_ENTSO_30_name = ["solar_ENTSO_30", "windon_ENTSO_30", "windoff_ENTSO_30"]


##### PARAMETERS #####
south_EU = ["Portugal", "Spain", "Italy"] 
north_EU = ["Norway", "Sweden", "Finland"]

reg = [south_EU, north_EU]
reg_name = ["South_EU", "North_EU"]
seas_name = ["summer", "autumn", "winter", "spring"]


countries = Dict("Albania" => "AL","Austria" => "AT", "Bosnia and Herzegovina" => "BA","Belgium" => "BE","Bulgaria" => "BG","Switzerland" => "CH","Cyprus" => "CY","Czech Republic" => "CZ","Germany" => "DE","Denmark" => "DK",
                        "Estonia" => "EE","Spain" => "ES", "Finland" => "FI", "France" => "FR", "Greece" => "GR", "Croatia" => "HR", "Hungary" => "HU", "Ireland" => "IE", "Italy" => "IT", "Lithuania" => "LT",
                        "Luxembourg" => "LU", "Latvia" => "LV", "Montenegro" => "ME", "North Macedonia" => "MK", "Malta" => "MT", "Netherlands" => "NL", "Norway" => "NO", "Poland" => "PL", "Portugal" => "PT",
                        "Romania" => "RO","Serbia" => "RS", "Sweden" => "SE", "Slovenia" => "SI", "Slovakia" => "SK", "Turkey" => "TR", "Ukraine" => "UA", "United Kingdom" => "UK")



#############################################################################################################################################################


function csv_boxplot_CF_ENTSO(data::DataFrame, data_name::String, countries::Dict{String, String}, BeginYear::Float64, EndYear::Float64, region::Vector{String}, region_name::String, Yearly::Bool, Season::String)
    data_region= dropmissing(data[(data.country .== countries[region[1]]) .| (data.country .== countries[region[2]]) .| (data.country .== countries[region[3]]) , [:year, :month, :day, :hour, :cf]], disallowmissing=true)
    data_region_annual= dropmissing(data_region[data_region.year .>= BeginYear .&& data_region.year .<= EndYear, [:month, :day, :hour, :cf]], disallowmissing=true)
    if Yearly
        output_file = "/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/$(data_name)_$(convert(Int, BeginYear))_$(convert(Int, EndYear))_year_$(region_name).csv"
        CSV.write(output_file, DataFrame(data=data_region_annual[!,:cf])) 
        return
        #return data_region_annual[!,:cf]
    end
    if Season == "winter"
        data_region_season= dropmissing(data_region_annual[(data_region_annual.month .== 12.0) .| (data_region_annual.month .== 1.0) .| (data_region_annual.month .== 2.0), [:day, :hour, :cf]], disallowmissing=true)
        output_file = "/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/$(data_name)_$(convert(Int, BeginYear))_$(convert(Int, EndYear))_$(Season)_$(region_name).csv"
        CSV.write(output_file, DataFrame(data=data_region_season[!,:cf])) 
        return
        #return data_region_season[!,:cf]
    end
    if Season == "summer"
        data_region_season= dropmissing(data_region_annual[data_region_annual.month .>= 6.0 .&& data_region_annual.month .<= 8.0, [:day, :hour, :cf]], disallowmissing=true)
        output_file = "/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/$(data_name)_$(convert(Int, BeginYear))_$(convert(Int, EndYear))_$(Season)_$(region_name).csv"
        CSV.write(output_file, DataFrame(data=data_region_season[!,:cf])) 
        return
        #return data_region_season[!,:cf]
    end
    if Season == "spring"
        data_region_season= dropmissing(data_region_annual[data_region_annual.month .>= 3.0 .&& data_region_annual.month .<= 5.0, [:day, :hour, :cf]], disallowmissing=true)
        output_file = "/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/$(data_name)_$(convert(Int, BeginYear))_$(convert(Int, EndYear))_$(Season)_$(region_name).csv"
        CSV.write(output_file, DataFrame(data=data_region_season[!,:cf])) 
        return
        #return data_region_season[!,:cf]
    end
    if Season == "autumn"
        data_region_season= dropmissing(data_region_annual[data_region_annual.month .>= 9.0 .&& data_region_annual.month .<= 11.0, [:day, :hour, :cf]], disallowmissing=true)
        output_file = "/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/$(data_name)_$(convert(Int, BeginYear))_$(convert(Int, EndYear))_$(Season)_$(region_name).csv"
        CSV.write(output_file, DataFrame(data=data_region_season[!,:cf])) 
        return
        #return data_region_season[!,:cf]
    else
        return "Season not known"
    end
end

#csv_boxplot_CF_ENTSO(windoff_ENTSO_25, "windoff_ENTSO_25", countries, 1986.0, 2015.0, south_EU, "south_EU", false, "winter")


function loop_csv_boxplot_CF_ENTSO(CF_ENTSO::Vector{DataFrame}, CF_ENTSO_name::Vector{String},seasons::Vector{String}, reg::Vector{Vector{String}}, reg_name::Vector{String})
    count_cf = 0
    for cf in CF_ENTSO
        count_cf +=1
        #println(CF_ENTSO_name[count_cf])
        count_r = 0
        for r in reg
            count_r +=1
            #println(reg_name[count_r])
            for s in seasons
                #println(s)
                csv_boxplot_CF_ENTSO(cf, CF_ENTSO_name[count_cf], countries, 1986.0, 2015.0, r, reg_name[count_r], false, s)
            end
            csv_boxplot_CF_ENTSO(cf, CF_ENTSO_name[count_cf], countries, 1986.0, 2015.0, r, reg_name[count_r], true, "summer")
        end
    end
end

#loop_csv_boxplot_CF_ENTSO(CF_ENTSO_30, CF_ENTSO_30_name, seas_name, reg, reg_name)


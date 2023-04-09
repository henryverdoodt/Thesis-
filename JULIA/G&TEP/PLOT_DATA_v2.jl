## Master Thesis: Analyse the impact of climate change on a G&TEP of EU Power System
# Author: Henry Verdoodt
# Last update: April 4, 2023

using CSV, DataFrames, Plots, Statistics

#############################################################################################################################################################
###############################                                      PLOT OF DEMAND                                           ###############################
#############################################################################################################################################################

# Load data from csv file
df = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/DEMAND/ED_CL_E_DAY_RCP8_ALADIN_CNRM_1971_2100.csv", DataFrame)

##### PARAMETERS #####
count = "Norway" # Name of Country

# Dictionaries with countries name, abreviation, number FOR DEMAND
countries_dict = Dict("Austria" => "AT", "Bosnia and Herzegovina" => "BA", "Belgium" => "BE", "Bulgaria" => "BG", "Switzerland" => "CH", "Czech Republic" => "CZ", "Germany" => "DE",
                      "Denmark" => "DK", "Estonia" => "EE",  "Greece" => "EL", "Spain" => "ES", "Finland" => "FI", "France" => "FR", "Croatia" => "HR",  "Hungary" => "HU",
                      "Ireland" => "IE", "Italy" => "IT", "Lithuania" => "LT", "Luxembourg" => "LU", "Latvia" => "LV", "Montenegro" => "ME", "North Macedonia" => "MK",
                      "Netherlands" => "NL",  "Norway" => "NO", "Poland" => "PL", "Portugal" => "PT", "Romania" => "RO", "Serbia" => "RS", "Sweden" => "SE", "Slovenia" => "SI",
                      "Slovakia" => "SK", "United Kingdom" => "UK", "Cyprus" => "CY")			  
countries = ["AT","BA","BE","BG","CH","CZ","DE","DK","EE","EL","ES","FI","FR","HR","HU","IE","IT","LT","LU","LV","ME","MK","NL","NO","PL","PT","RO","RS","SE","SI","SK","UK","CY"]
num_countries = Dict(zip(countries, 2:length(countries)+1))
num_countries[countries_dict[count]]

# select the year(s) you want to plot
year1 = ["1971"]
year2 = ["2030"]
year3 = ["2100"]

# filter the dataframe to only include rows for the selected year(s) FOR DEMAND
rows1 = findall(row -> split(row, "/")[1] in Set(year1), df[!, :Date])
df_year1 = df[rows1, :]
rows2 = findall(row -> split(row, "/")[1] in Set(year2), df[!, :Date])
df_year2 = df[rows2, :]
rows3 = findall(row -> split(row, "/")[1] in Set(year3), df[!, :Date])
df_year3 = df[rows3, :]

# extract the date column as an array of strings
dates = df_year1[!, :Date]

# extract the data for each country as an array of floats
data1 = Matrix(df_year1[:, num_countries[countries_dict[count]]:num_countries[countries_dict[count]]])
data2 = Matrix(df_year2[:, num_countries[countries_dict[count]]:num_countries[countries_dict[count]]])
data3 = Matrix(df_year3[:, num_countries[countries_dict[count]]:num_countries[countries_dict[count]]])

# Plot of Demand evolution
locations = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
p1 = plot(dates, data1, ylabel="Demand (MWh)", title= "Evolution of Demand of $count", xticks=(locations, month_names), label=year1[1], legend=:bottomright)
p2 = plot!(dates, data2, label=year2[1], legend=:bottomright)
p3 = plot!(dates, data3, label=year3[1], legend=:bottomright)

# UNCOMMENT TO SAVE the plot to the IMAGES file 
#savefig("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/DEMAND/IMAGES/Evolution_Demand_$count.png")




#############################################################################################################################################################
###############################                                      PLOT OF SOLAR                                            ###############################
#############################################################################################################################################################

using CSV, DataFrames, Plots, Statistics

# Load data from csv file
solar = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/SOLAR/SOLAR_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)

##### PARAMETERS #####
year = ["1972"] # 1952-2100
s = "winter" # plot for season s ("summer", "autumn", "winter", "spring")

# Dictionaries with countries name, abreviation, number FOR SOLAR
country_codes = Dict("Albania" => "AL", "Austria" => "AT", "Bosnia and Herzegovina" => "BA", "Belgium" => "BE", "Bulgaria" => "BG", "Switzerland" => "CH","Cyprus" => "CY","Czech Republic" => "CZ",
                     "Germany" => "DE","Denmark" => "DK", "Estonia" => "EE","Greece" => "EL", "Spain" => "ES", "Finland" => "FI", "France" => "FR", "Croatia" => "HR","Hungary" => "HU",
                     "Ireland" => "IE","Iceland" => "IS","Italy" => "IT", "Liechtenstein" => "LI", "Lithuania" => "LT", "Luxembourg" => "LU", "Latvia" => "LV","Montenegro" => "ME",
                     "North Macedonia" => "MK","Malta" => "MT","Netherlands" => "NL","Norway" => "NO", "Poland" => "PL","Portugal" => "PT","Romania" => "RO","Serbia" => "RS","Sweden" => "SE",
                     "Slovenia" => "SI", "Slovakia" => "SK","Turkey" => "TR","United Kingdom" => "UK")
countries2 = ["AL","AT","BA","BE","BG","CH","CY","CZ","DE","DK","EE","EL","ES","FI","FR","HR","HU","IE","IS","IT","LI","LT","LU","LV","ME","MK","MT","NL","NO","PL","PT","RO","RS","SE","SI","SK","TR","UK"]

num_countries2 = Dict(zip(countries2, 2:length(countries)+1))
num_countries2[country_codes[count]]

# Seasons and Timestep parameter
summer = ["06", "07", "08"]; autumn = ["09", "10", "11"]; winter = ["12", "01", "02"]; spring = ["03", "04", "05"]
t1 = ["01:30:00"]; t2 = ["04:30:00"]; t3 = ["07:30:00"]; t4 = ["10:30:00"]; t5 = ["13:30:00"]; t6 = ["16:30:00"]; t7 = ["19:30:00"]; t8 = ["22:30:00"]

#############################################################################################################################################################
# filter the dataframe to only include rows for the selected year(s) FOR SOLAR
rows1 = findall(row -> split(row, "-")[1] in Set(year), solar[!, :Date])
solar_year = solar[rows1, :]

solar_summer = solar_year[findall(row -> split(row, "-")[2] in Set(summer), solar_year[!, :Date]), :]
solar_autumn = solar_year[findall(row -> split(row, "-")[2] in Set(autumn), solar_year[!, :Date]), :]
solar_winter = solar_year[findall(row -> split(row, "-")[2] in Set(winter), solar_year[!, :Date]), :]
solar_spring = solar_year[findall(row -> split(row, "-")[2] in Set(spring), solar_year[!, :Date]), :]

seasons = Dict("summer" => solar_summer, "autumn" => solar_autumn, "winter" => solar_winter, "spring" => solar_spring)

function plot_solar_curve(season::AbstractString, df::DataFrame)
    t = [t1, t2, t3, t4, t5, t6, t7, t8]
    solar_means = Dict("ES" => Float64[], "NO" => Float64[], "BE" => Float64[], "BG" => Float64[])
    
    for ti in t
        subset = df[findall(row -> split(row, " ")[2] in Set(ti), df[!, :Date]), :]
        for key in keys(solar_means)
            push!(solar_means[key], mean(skipmissing(subset[!, Symbol(key)]))) 
        end
    end
    
    p1 = plot(title="$season Solar yield curve of year $(year[1])", xlabel="Time", ylabel="Solar_CF", legend=:bottomright, markershape=:circle)
    for (key, values) in solar_means
        plot!(p1, values, label=key)
    end
    return p1
end

plot_solar_curve(s, seasons[s])

# UNCOMMENT TO SAVE the plot to the IMAGES file 
#savefig("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/SOLAR/IMAGES/$season Solar yield curve of year $year.png")


using CSV, DataFrames, Plots, Statistics

# Load data from csv file
solar = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/SOLAR/SOLAR_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)

# Define the season and country as parameters
season = "winter"
country = "ES"

# Seasons and Timestep parameter
summer = ["06", "07", "08"]; autumn = ["09", "10", "11"]; winter = ["12", "01", "02"]; spring = ["03", "04", "05"]
t1 = ["01:30:00"]; t2 = ["04:30:00"]; t3 = ["07:30:00"]; t4 = ["10:30:00"]; t5 = ["13:30:00"]; t6 = ["16:30:00"]; t7 = ["19:30:00"]; t8 = ["22:30:00"]

# filter the dataframe to only include rows for the selected year(s) FOR SOLAR
year1 = "1972"; year2 = "2030"; year3 = "2100"
rows1 = findall(row -> split(row, "-")[1] in Set([year1, year2, year3]), solar[!, :Date])
solar_years = solar[rows1, :]

# Filter the dataframe based on season and country
solar_country = solar_years[findall(row -> split(row, "-")[2] in Set(eval(country)), solar_years[!, :Date]), :]
solar_seasons = Dict("summer" => solar_summer, "autumn" => solar_autumn, "winter" => solar_winter, "spring" => solar_spring)
df = solar_seasons[season][findall(row -> split(row, "-")[2] in Set(eval(season)), solar_seasons[season][!, :Date]), :]

function plot_solar_curve(country::AbstractString, season::AbstractString, df::DataFrame)
    t = [t1, t2, t3, t4, t5, t6, t7, t8]
    solar_means = Dict(year1 => Dict("ES" => Float64[], "NO" => Float64[], "BE" => Float64[], "BG" => Float64[]),
                       year2 => Dict("ES" => Float64[], "NO" => Float64[], "BE" => Float64[], "BG" => Float64[]),
                       year3 => Dict("ES" => Float64[], "NO" => Float64[], "BE" => Float64[], "BG" => Float64[]))

    for ti in t
        subset = df[findall(row -> split(row, " ")[2] in Set(ti), df[!, :Date]), :]
        for (year, solar_means_year) in solar_means
            for key in keys(solar_means_year)
                push!(solar_means[year][key], mean(skipmissing(subset[!, Symbol(key)]))) 
            end
        end
    end
    
    p1 = plot(title="$country $season Solar yield curve for $year1, $year2 and $year3", xlabel="Time", ylabel="Solar_CF", legend=:bottomright, markershape=:circle)
    for (year, solar_means_year) in solar_means
        for (key, values) in solar_means_year
            plot!(p1, values, label="$key $year")
        end
    end
    return p1
end

plot_solar_curve(country, season, solar_country)


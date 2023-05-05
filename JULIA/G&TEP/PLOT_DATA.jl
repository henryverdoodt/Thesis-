## Master Thesis: Analyse the impact of climate change on a G&TEP of EU Power System
# Author: Henry Verdoodt
# Last update: April 4, 2023
using Pkg
Pkg.add("PlotlyJS")


# Include Files
include("/Users/henryverdoodt/Documents/CODE/JULIA/G&TEP/REFORMAT_DATA.jl")

#############################################################################################################################################################
#############################################################################################################################################################
###############################                                        COPERNICUS                                             ###############################
#############################################################################################################################################################
#############################################################################################################################################################

using CSV, DataFrames, Plots, Statistics

#############################################################################################################################################################
###############################                                      PLOT OF DEMAND                                           ###############################
#############################################################################################################################################################

# Load data from csv file
df = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/DEMAND/ED_CL_E_DAY_RCP8_ALADIN_CNRM_1971_2100.csv", DataFrame)

##### PARAMETERS #####
count = "Belgium" # Name of Country

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
year = ["2100"] # 1952-2100
s = "winter" # plot for season s ("summer", "autumn", "winter", "spring")
count = "Belgium"

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
rows1 = findall(row -> begin
                      year = parse(Int, split(row, "-")[1])
                      2071 <= year <= 2100
                  end, solar[!, :Date])

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
    loc=[1,2,3,4,5,6,7,8]
    hours=["1:30h","4:30h","7:30h","10:30h", "13:30h","16:30h","19:30h","22:30h"]
    p1 = plot(title="$season Solar yield curve of year $(year[1])", xlabel="Time", ylabel="Solar_CF", xticks=(loc, hours), legend=:bottomright, markershape=:circle)
    for (key, values) in solar_means
        plot!(p1, values, label=key)
    end
    return p1
end

plot_solar_curve(s, seasons[s])

# UNCOMMENT TO SAVE the plot to the IMAGES file 
#savefig("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/SOLAR/IMAGES/$season Solar yield curve of year $year.png")

#############################################################################################################################################################
###############################                                      BOXPLOT SOLAR                                            ###############################
#############################################################################################################################################################
using CSV, DataFrames, Plots, Statistics

# Load data from csv file
solar = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/SOLAR/SOLAR_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)

summer = ["06", "07", "08"]; autumn = ["09", "10", "11"]; winter = ["12", "01", "02"]; spring = ["03", "04", "05"]

#############################################################################################################################################################
# filter the dataframe to only include rows for the selected year(s) FOR SOLAR
rows1 = findall(row -> begin
                      year = parse(Int, split(row, "-")[1])
                      2071 <= year <= 2100
                  end, solar[!, :Date])

solar_year = solar[rows1, :]

# Dataframe with all solar cf per country for specific years
solar_summer = solar_year[findall(row -> split(row, "-")[2] in Set(summer), solar_year[!, :Date]), :]
solar_autumn = solar_year[findall(row -> split(row, "-")[2] in Set(autumn), solar_year[!, :Date]), :]
solar_winter = solar_year[findall(row -> split(row, "-")[2] in Set(winter), solar_year[!, :Date]), :]
solar_spring = solar_year[findall(row -> split(row, "-")[2] in Set(spring), solar_year[!, :Date]), :]

# Vector with data of solar cf for South EU
solar_year_south = vcat(solar_year[!, :PT], solar_year[!, :ES], solar_year[!, :IT])
solar_summer_south = vcat(solar_summer[!, :PT], solar_summer[!, :ES], solar_summer[!, :IT])
solar_autumn_south = vcat(solar_autumn[!, :PT], solar_autumn[!, :ES], solar_autumn[!, :IT])
solar_winter_south = vcat(solar_winter[!, :PT], solar_winter[!, :ES], solar_winter[!, :IT])
solar_spring_south = vcat(solar_spring[!, :PT], solar_spring[!, :ES], solar_spring[!, :IT])

# Vector with data of solar cf for North EU
solar_year_north = filter(!isnan, vcat(solar_year[!, :NO], solar_year[!, :SE], solar_year[!, :FI]))
solar_summer_north = filter(!isnan, vcat(solar_summer[!, :NO], solar_summer[!, :SE], solar_summer[!, :FI]))
solar_autumn_north = filter(!isnan, vcat(solar_autumn[!, :NO], solar_autumn[!, :SE], solar_autumn[!, :FI]))
solar_winter_north = filter(!isnan, vcat(solar_winter[!, :NO], solar_winter[!, :SE], solar_winter[!, :FI]))
solar_spring_north = filter(!isnan, vcat(solar_spring[!, :NO], solar_spring[!, :SE], solar_spring[!, :FI]))

# Boxplot of solar cf for South EU
boxplot_solar_south_annual_coper = boxplot(solar_year_south, xlabel="ANNUAL", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_solar_south_summer_coper = boxplot(solar_summer_south, xlabel="SUMMER", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_solar_south_autumn_coper = boxplot(solar_autumn_south, xlabel="AUTUMN", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_solar_south_winter_coper = boxplot(solar_winter_south, xlabel="WINTER", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_solar_south_spring_coper = boxplot(solar_spring_south, xlabel="SPRING", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"

# Boxplot of solar cf for North EU
boxplot_solar_north_annual_coper = boxplot(solar_year_north, xlabel="ANNUAL", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_solar_north_summer_coper = boxplot(solar_summer_north, xlabel="SUMMER", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_solar_north_autumn_coper = boxplot(solar_autumn_north, xlabel="AUTUMN", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_solar_north_winter_coper = boxplot(solar_winter_north, xlabel="WINTER", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_solar_north_spring_coper = boxplot(solar_spring_north, xlabel="SPRING", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"

plot(boxplot_solar_south_annual_coper, boxplot_solar_south_summer_coper, boxplot_solar_south_autumn_coper, boxplot_solar_south_winter_coper, boxplot_solar_south_spring_coper, layout=(1,5))
plot(boxplot_solar_north_annual_coper, boxplot_solar_north_summer_coper, boxplot_solar_north_autumn_coper, boxplot_solar_north_winter_coper, boxplot_solar_north_spring_coper, layout=(1,5))



#############################################################################################################################################################
###############################                                      PLOT OF WIND                                             ###############################
#############################################################################################################################################################

using CSV, DataFrames, Plots, Statistics

# Load data from csv file
windon = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDON/WINDON_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)

##### PARAMETERS #####
year = ["2100"] # 1952-2100
s = "winter" # plot for season s ("summer", "autumn", "winter", "spring")

# Dictionaries with countries name, abreviation, number FOR WIND
country_codes = Dict("Albania" => "AL", "Austria" => "AT", "Bosnia and Herzegovina" => "BA", "Belgium" => "BE", "Bulgaria" => "BG", "Switzerland" => "CH","Cyprus" => "CY","Czech Republic" => "CZ",
                     "Germany" => "DE","Denmark" => "DK", "Estonia" => "EE","Greece" => "EL", "Spain" => "ES", "Finland" => "FI", "France" => "FR", "Croatia" => "HR","Hungary" => "HU",
                     "Ireland" => "IE","Iceland" => "IS","Italy" => "IT", "Liechtenstein" => "LI", "Lithuania" => "LT", "Luxembourg" => "LU", "Latvia" => "LV","Montenegro" => "ME",
                     "North Macedonia" => "MK","Malta" => "MT","Netherlands" => "NL","Norway" => "NO", "Poland" => "PL","Portugal" => "PT","Romania" => "RO","Serbia" => "RS","Sweden" => "SE",
                     "Slovenia" => "SI", "Slovakia" => "SK","Turkey" => "TR","United Kingdom" => "UK")
countries2 = ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "EL", "ES", "FI", "FR", "HR", "HU", "IE", "IS", "IT", "LI", "LT", "LU", "LV", "ME", "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UK"]

num_countries2 = Dict(zip(countries2, 2:length(countries)+1))
num_countries2[country_codes[count]]

# Seasons and Timestep parameter
summer = ["06", "07", "08"]; autumn = ["09", "10", "11"]; winter = ["12", "01", "02"]; spring = ["03", "04", "05"]
t1 = ["00:00:00"]; t2 = ["03:00:00"]; t3 = ["06:00:00"]; t4 = ["09:00:00"]; t5 = ["12:00:00"]; t6 = ["15:00:00"]; t7 = ["18:00:00"]; t8 = ["21:00:00"]

#############################################################################################################################################################
# filter the dataframe to only include rows for the selected year(s) FOR WIND
rows1 = findall(row -> split(row, "-")[1] in Set(year), windon[!, :Date])
windon_year = windon[rows1, :]

windon_summer = windon_year[findall(row -> split(row, "-")[2] in Set(summer), windon_year[!, :Date]), :]
windon_autumn = windon_year[findall(row -> split(row, "-")[2] in Set(autumn), windon_year[!, :Date]), :]
windon_winter = windon_year[findall(row -> split(row, "-")[2] in Set(winter), windon_year[!, :Date]), :]
windon_spring = windon_year[findall(row -> split(row, "-")[2] in Set(spring), windon_year[!, :Date]), :]

seasons = Dict("summer" => windon_summer, "autumn" => windon_autumn, "winter" => windon_winter, "spring" => windon_spring)

function plot_windon_curve(season::AbstractString, df::DataFrame)
    t = [t1, t2, t3, t4, t5, t6, t7, t8]
    windon_means = Dict("ES" => Float64[], "NO" => Float64[], "BE" => Float64[], "BG" => Float64[])
    
    for ti in t
        subset = df[findall(row -> split(row, " ")[2] in Set(ti), df[!, :Date]), :]
        for key in keys(windon_means)
            push!(windon_means[key], mean(skipmissing(subset[!, Symbol(key)]))) 
        end
    end
    loc=[1,2,3,4,5,6,7,8]
    hours=["0h","3h","6h","9h", "12h","15h","18h","21h"]
    p1 = plot(title="$season Windon yield curve of year $(year[1])", xlabel="Time", ylabel="Windon_CF", xticks=(loc, hours), legend=:bottomright, markershape=:circle)
    for (key, values) in windon_means
        plot!(p1, values, label=key)
    end
    return p1
end

plot_windon_curve(s, seasons[s])

# UNCOMMENT TO SAVE the plot to the IMAGES file 
#savefig("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDON/IMAGES/$season Windon yield curve of year $year.png")


#############################################################################################################################################################
###############################                                      BOXPLOT WINDON COPERNICUS                                ###############################
#############################################################################################################################################################
using CSV, DataFrames, Plots, Statistics

# Load data from csv file
windon = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDON/WINDON_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)

summer = ["06", "07", "08"]; autumn = ["09", "10", "11"]; winter = ["12", "01", "02"]; spring = ["03", "04", "05"]

#############################################################################################################################################################
# filter the dataframe to only include rows for the selected year(s) FOR WINDON
rows1 = findall(row -> begin
                      year = parse(Int, split(row, "-")[1])
                      2071 <= year <= 2100
                  end, windon[!, :Date])

windon_year = windon[rows1, :]

# Dataframe with all windon cf per country for specific years
windon_summer = windon_year[findall(row -> split(row, "-")[2] in Set(summer), windon_year[!, :Date]), :]
windon_autumn = windon_year[findall(row -> split(row, "-")[2] in Set(autumn), windon_year[!, :Date]), :]
windon_winter = windon_year[findall(row -> split(row, "-")[2] in Set(winter), windon_year[!, :Date]), :]
windon_spring = windon_year[findall(row -> split(row, "-")[2] in Set(spring), windon_year[!, :Date]), :]

# Vector with data of windon cf for South EU
windon_year_south = vcat(windon_year[!, :PT], windon_year[!, :ES], windon_year[!, :IT])
windon_summer_south = vcat(windon_summer[!, :PT], windon_summer[!, :ES], windon_summer[!, :IT])
windon_autumn_south = vcat(windon_autumn[!, :PT], windon_autumn[!, :ES], windon_autumn[!, :IT])
windon_winter_south = vcat(windon_winter[!, :PT], windon_winter[!, :ES], windon_winter[!, :IT])
windon_spring_south = vcat(windon_spring[!, :PT], windon_spring[!, :ES], windon_spring[!, :IT])

# Vector with data of windon cf for North EU
windon_year_north = filter(!isnan, vcat(windon_year[!, :NO], windon_year[!, :SE], windon_year[!, :FI]))
windon_summer_north = filter(!isnan, vcat(windon_summer[!, :NO], windon_summer[!, :SE], windon_summer[!, :FI]))
windon_autumn_north = filter(!isnan, vcat(windon_autumn[!, :NO], windon_autumn[!, :SE], windon_autumn[!, :FI]))
windon_winter_north = filter(!isnan, vcat(windon_winter[!, :NO], windon_winter[!, :SE], windon_winter[!, :FI]))
windon_spring_north = filter(!isnan, vcat(windon_spring[!, :NO], windon_spring[!, :SE], windon_spring[!, :FI]))

# Boxplot of windon cf for South EU
boxplot_windon_south_annual_coper = boxplot(windon_year_south, xlabel="ANNUAL", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windon_south_summer_coper = boxplot(windon_summer_south, xlabel="SUMMER", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windon_south_autumn_coper = boxplot(windon_autumn_south, xlabel="AUTUMN", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windon_south_winter_coper = boxplot(windon_winter_south, xlabel="WINTER", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windon_south_spring_coper = boxplot(windon_spring_south, xlabel="SPRING", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"

# Boxplot of windon cf for North EU
boxplot_windon_north_annual_coper = boxplot(windon_year_north, xlabel="ANNUAL", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windon_north_summer_coper = boxplot(windon_summer_north, xlabel="SUMMER", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windon_north_autumn_coper = boxplot(windon_autumn_north, xlabel="AUTUMN", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windon_north_winter_coper = boxplot(windon_winter_north, xlabel="WINTER", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windon_north_spring_coper = boxplot(windon_spring_north, xlabel="SPRING", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"

plot(boxplot_windon_south_annual_coper, boxplot_windon_south_summer_coper, boxplot_windon_south_autumn_coper, boxplot_windon_south_winter_coper, boxplot_windon_south_spring_coper, layout=(1,5))
plot(boxplot_windon_north_annual_coper, boxplot_windon_north_summer_coper, boxplot_windon_north_autumn_coper, boxplot_windon_north_winter_coper, boxplot_windon_north_spring_coper, layout=(1,5))





#############################################################################################################################################################
###############################                                      BOXPLOT WINDOFF COPERNICUS                                ###############################
#############################################################################################################################################################
using CSV, DataFrames, Plots, Statistics

# Load data from csv file
windoff = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/COPERNICUS/WINDOFF/WINDOFF_MARITIME_CL_CF_3H_RCP8_ALADIN_CNRM_1952_2100.csv", DataFrame)

summer = ["06", "07", "08"]; autumn = ["09", "10", "11"]; winter = ["12", "01", "02"]; spring = ["03", "04", "05"]

#############################################################################################################################################################
# filter the dataframe to only include rows for the selected year(s) FOR WINDOFF
rows1 = findall(row -> begin
                      year = parse(Int, split(row, "-")[1])
                      2071 <= year <= 2100
                  end, windoff[!, :Date])

windoff_year = windoff[rows1, :]

# Dataframe with all windoff cf per country for specific years
windoff_summer = windoff_year[findall(row -> split(row, "-")[2] in Set(summer), windoff_year[!, :Date]), :]
windoff_autumn = windoff_year[findall(row -> split(row, "-")[2] in Set(autumn), windoff_year[!, :Date]), :]
windoff_winter = windoff_year[findall(row -> split(row, "-")[2] in Set(winter), windoff_year[!, :Date]), :]
windoff_spring = windoff_year[findall(row -> split(row, "-")[2] in Set(spring), windoff_year[!, :Date]), :]

# Vector with data of windoff cf for South EU
windoff_year_south = vcat(windoff_year[!, :PT], windoff_year[!, :ES], windoff_year[!, :IT])
windoff_summer_south = vcat(windoff_summer[!, :PT], windoff_summer[!, :ES], windoff_summer[!, :IT])
windoff_autumn_south = vcat(windoff_autumn[!, :PT], windoff_autumn[!, :ES], windoff_autumn[!, :IT])
windoff_winter_south = vcat(windoff_winter[!, :PT], windoff_winter[!, :ES], windoff_winter[!, :IT])
windoff_spring_south = vcat(windoff_spring[!, :PT], windoff_spring[!, :ES], windoff_spring[!, :IT])

# Vector with data of windoff cf for North EU
windoff_year_north = filter(!isnan, vcat(windoff_year[!, :NO], windoff_year[!, :SE], windoff_year[!, :FI]))
windoff_summer_north = filter(!isnan, vcat(windoff_summer[!, :NO], windoff_summer[!, :SE], windoff_summer[!, :FI]))
windoff_autumn_north = filter(!isnan, vcat(windoff_autumn[!, :NO], windoff_autumn[!, :SE], windoff_autumn[!, :FI]))
windoff_winter_north = filter(!isnan, vcat(windoff_winter[!, :NO], windoff_winter[!, :SE], windoff_winter[!, :FI]))
windoff_spring_north = filter(!isnan, vcat(windoff_spring[!, :NO], windoff_spring[!, :SE], windoff_spring[!, :FI]))

# Boxplot of windoff cf for South EU
boxplot_windoff_south_annual_coper = boxplot(windoff_year_south, xlabel="ANNUAL", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"   # title="WINDOFF CF SOUTH EU PROJECTIONS"
boxplot_windoff_south_summer_coper = boxplot(windoff_summer_south, xlabel="SUMMER", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windoff_south_autumn_coper = boxplot(windoff_autumn_south, xlabel="AUTUMN", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windoff_south_winter_coper = boxplot(windoff_winter_south, xlabel="WINTER", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windoff_south_spring_coper = boxplot(windoff_spring_south, xlabel="SPRING", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"

# Boxplot of windoff cf for North EU
boxplot_windoff_north_annual_coper = boxplot(windoff_year_north, xlabel="ANNUAL", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windoff_north_summer_coper = boxplot(windoff_summer_north, xlabel="SUMMER", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windoff_north_autumn_coper = boxplot(windoff_autumn_north, xlabel="AUTUMN", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windoff_north_winter_coper = boxplot(windoff_winter_north, xlabel="WINTER", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"
boxplot_windoff_north_spring_coper = boxplot(windoff_spring_north, xlabel="SPRING", whisker_width=0.5, marker_size=3, color="red", label = "PROJ") # ylabel="Variable Value"

plot(boxplot_windoff_south_annual_coper, boxplot_windoff_south_summer_coper, boxplot_windoff_south_autumn_coper, boxplot_windoff_south_winter_coper, boxplot_windoff_south_spring_coper, layout=(1,5))
plot(boxplot_windoff_north_annual_coper, boxplot_windoff_north_summer_coper, boxplot_windoff_north_autumn_coper, boxplot_windoff_north_winter_coper, boxplot_windoff_north_spring_coper, layout=(1,5))



plot(boxplot_windoff_south_annual_ref, boxplot_windoff_south_annual_coper, boxplot_windoff_south_summer_ref, boxplot_windoff_south_summer_coper, boxplot_windoff_south_autumn_ref, boxplot_windoff_south_autumn_coper, boxplot_windoff_south_spring_ref, boxplot_windoff_south_spring_coper, boxplot_windoff_south_winter_ref,boxplot_windoff_south_winter_coper, layout=(1,10), size=(1200,500))











#############################################################################################################################################################
#############################################################################################################################################################
###############################                                             ENTSO                                             ###############################
#############################################################################################################################################################
#############################################################################################################################################################

#############################################################################################################################################################
###############################                                      PLOT OF DEMAND                                           ###############################
#############################################################################################################################################################


using CSV, DataFrames, Plots, Statistics, Dates


# Read the CSV file into a DataFrame
demand = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/DEMAND/PECD-country-demand_national_estimates-2025.csv", DataFrame)
solar = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/SOLAR/PECD-2021.3-country-LFSolarPV-2025.csv", DataFrame)
windon = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDON/PECD-2021.3-country-Onshore-2025.csv", DataFrame)
windoff = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDOFF/PECD-2021.3-country-Offshore-2025.csv", DataFrame)

##### PARAMETERS #####
count = "Norway" # Name of Country

countries = Dict("Albania" => "AL","Austria" => "AT", "Bosnia and Herzegovina" => "BA","Belgium" => "BE","Bulgaria" => "BG","Switzerland" => "CH","Cyprus" => "CY","Czech Republic" => "CZ","Germany" => "DE","Denmark" => "DK",
                        "Estonia" => "EE","Spain" => "ES", "Finland" => "FI", "France" => "FR", "Greece" => "GR", "Croatia" => "HR", "Hungary" => "HU", "Ireland" => "IE", "Italy" => "IT", "Lithuania" => "LT",
                        "Luxembourg" => "LU", "Latvia" => "LV", "Montenegro" => "ME", "North Macedonia" => "MK", "Malta" => "MT", "Netherlands" => "NL", "Norway" => "NO", "Poland" => "PL", "Portugal" => "PT",
                        "Romania" => "RO","Serbia" => "RS", "Sweden" => "SE", "Slovenia" => "SI", "Slovakia" => "SK", "Turkey" => "TR", "Ukraine" => "UA", "United Kingdom" => "UK")

# Get the unique countries in the 'country' column
countries_demand = unique(demand.country) # ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "ME", "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UA", "UK"]
countries_solar= unique(solar.country)    # ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "ME", "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UA", "UK"]
countries_windon = unique(windon.country) # MT is in solar & demand but not in windon => # ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "ME", "MK", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UA", "UK"]
countries_windoff= unique(windoff.country) # ["BE", "DE", "DK", "FI", "FR", "IE", "IT", "NL", "PT", "SE", "UK"]


# filter the data for a specific country and select only the relevant columns
es_demand = demand[demand.country .== countries[count], [:year, :month, :day, :hour, :dem_MW]]
es_demand = dropmissing(es_demand, disallowmissing=true)

# create a new column for the date/time
es_demand.date = DateTime.(es_demand.year, es_demand.month, es_demand.day, es_demand.hour, 0, 0)
date = es_demand.date

# Make a daily average of the dem_MW column
daily_avg = combine(groupby(es_demand, [:year, :month, :day]), :dem_MW => mean => :dem_MW)
daily_avg.dem_MW .= daily_avg.dem_MW .* 24

# Subset the data to years 1982, 2000, and 2016
subset1 = daily_avg[(daily_avg.year .== 1982.0), :]
subset2 = daily_avg[(daily_avg.year .== 2000.0), :]
subset3 = daily_avg[(daily_avg.year .== 2016.0), :]

plot(1:size(subset1, 1), subset1.dem_MW, xlabel="days", ylabel="Demand (MWh)", title= "Evolution of Demand of $count", label="1982",)
plot!(1:size(subset2, 1), subset2.dem_MW, xlabel="days", ylabel="Demand (MWh)", title= "Evolution of Demand of $count", label="2000",)
plot!(1:size(subset3, 1), subset3.dem_MW, xlabel="days", ylabel="Demand (MWh)", title= "Evolution of Demand of $count", label="2016",)

#############################################################################################################################################################

############################ DATA PARAMETERS ############################
countries = ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO", "CH", "FI", "IE", "IT", "AT", "PT", "SE"]   #["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO"]
y = 2009.0    # demand_entso: 1982.0 - 2016.0    and   solar,windon,windoff_entso: 1982.0 - 2019.0
year = 2030
representative_years = [1995.0, 2008.0, 2009.0]
weights = [0.233, 0.367, 0.4]
aggregate_3h = true

countries_demand = countries_demand_entso       
countries_solar = countries_solar_entso         # countries_solar_coper_cnrm
countries_windon = countries_windon_entso       # countries_windon_coper_cnrm
countries_windoff = countries_windoff_entso     # countries_windoff_coper_cnrm

countries_dem = ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO", "CH", "FI", "IE", "AT", "SE"]
countries_dem2 = ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO", "CH", "FI", "IE", "AT", "SE", "IT", "PT"]

dem = reformat_entso_demand(demand_entso, countries_dem, countries_demand_entso, representative_years, weights, aggregate_3h)
dem_PT_IT = reformat_demand_PT_IT(demand_PT_IT, aggregate_3h)
dem_RP = hcat(dem, dem_PT_IT)

dem_1982 = reformat_entso_demand(demand_entso, countries_dem2, countries_demand_entso, [1982.0, 1982.0], [0.5,0.5], aggregate_3h)

plot(dem_1982, title="ENTSO demand 1982", label=names(dem_1982))

X = :PT        # ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO", "CH", "FI", "IE", "IT", "AT", "PT", "SE"]
plot(dem_1982[!, X], label="1982") 
plot!(dem_RP[!, X], label="RP")

sol_RP = reformat_entso_solar(solar_entso, countries, countries_solar_entso, representative_years, weights, aggregate_3h)
sol_1982 = reformat_entso_solar(solar_entso, countries, countries_solar_entso, [1982.0], [1.0], aggregate_3h)

won_RP = reformat_entso_windon(windon_entso, countries, countries_windon_entso, representative_years, weights, aggregate_3h) 
won_1982 = reformat_entso_windon(windon_entso, countries, countries_windon_entso, [1982.0], [1.0], aggregate_3h) 

woff_RP = reformat_entso_windoff(windoff_entso, countries, countries_windoff_entso, representative_years, weights, aggregate_3h)
woff_1982 = reformat_entso_windoff(windoff_entso, countries, countries_windoff_entso, [1982.0], [1.0], aggregate_3h)


X = :BE        # ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO", "CH", "FI", "IE", "IT", "AT", "PT", "SE"]
plot(won_1982[!, X], label="1982") 
plot!(won_RP[!, X], label="RP")


#############################################################################################################################################################
###############################                                      PLOT OF SOLAR                                            ###############################
#############################################################################################################################################################

using CSV, DataFrames, Plots, Statistics, Dates


# Read the CSV file into a DataFrame
solar = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/SOLAR/PECD-2021.3-country-LFSolarPV-2025.csv", DataFrame)

##### PARAMETERS #####
count = "Spain" # Name of Country 
season = "summer"
y1 = 1982.0
y2 = 2000.0
y3 = 2016.0


countries = Dict("Albania" => "AL","Austria" => "AT", "Bosnia and Herzegovina" => "BA","Belgium" => "BE","Bulgaria" => "BG","Switzerland" => "CH","Cyprus" => "CY","Czech Republic" => "CZ","Germany" => "DE","Denmark" => "DK",
                        "Estonia" => "EE","Spain" => "ES", "Finland" => "FI", "France" => "FR", "Greece" => "GR", "Croatia" => "HR", "Hungary" => "HU", "Ireland" => "IE", "Italy" => "IT", "Lithuania" => "LT",
                        "Luxembourg" => "LU", "Latvia" => "LV", "Montenegro" => "ME", "North Macedonia" => "MK", "Malta" => "MT", "Netherlands" => "NL", "Norway" => "NO", "Poland" => "PL", "Portugal" => "PT",
                        "Romania" => "RO","Serbia" => "RS", "Sweden" => "SE", "Slovenia" => "SI", "Slovakia" => "SK", "Turkey" => "TR", "Ukraine" => "UA", "United Kingdom" => "UK")

seasons = Dict("summer" => [6.0, 7.0, 8.0], "autumn" => [9.0, 10.0, 11.0], "winter" => [12.0, 1.0, 2.0], "spring" => [3.0, 4.0, 5.0])

# filter the data for a specific country and select only the relevant columns
solar= dropmissing(solar[(solar.country .== countries[count]) , [:year, :month, :day, :hour, :cf]], disallowmissing=true)

# Filter the DataFrame to keep only the summer days of year y
df_subset1 = filter(row -> row.year == y1 && row.month in seasons[season], solar)
df_subset2 = filter(row -> row.year == y2 && row.month in seasons[season], solar)
df_subset3 = filter(row -> row.year == y3 && row.month in seasons[season], solar)

# group the rows by hour
hourly_avg1 = combine(groupby(df_subset1, :hour), :cf => mean => :cf)
hourly_avg2 = combine(groupby(df_subset2, :hour), :cf => mean => :cf)
hourly_avg3 = combine(groupby(df_subset3, :hour), :cf => mean => :cf)

# create a new column for the date/time
#solar.date = DateTime.(solar.year, solar.month, solar.day, solar.hour, 0, 0)
#date = solar.date

# Plot
plot(1:size(hourly_avg1, 1), hourly_avg1.cf, xlabel="hours", ylabel="Capacity Factor", title= "Evolution of Capacity Factor of $count for $season", label="1982",)
plot!(1:size(hourly_avg2, 1), hourly_avg2.cf, xlabel="hours", ylabel="Capacity Factor", title= "Evolution of Capacity Factor of $count for $season", label="2000",)
plot!(1:size(hourly_avg3, 1), hourly_avg3.cf, xlabel="hours", ylabel="Capacity Factor", title= "Evolution of Capacity Factor of $count for $season", label="2016",)


#############################################################################################################################################################
###############################                                      BOXPLOT SOLAR                                            ###############################
#############################################################################################################################################################
using CSV, DataFrames, Plots, Statistics, Dates


# Read the CSV file into a DataFrame
solar = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/SOLAR/PECD-2021.3-country-LFSolarPV-2025.csv", DataFrame)

##### PARAMETERS #####
count_south = ["Portugal", "Spain", "Italy"] 
count_north = ["Norway", "Sweden", "Finland"] 
season = "summer"


countries = Dict("Albania" => "AL","Austria" => "AT", "Bosnia and Herzegovina" => "BA","Belgium" => "BE","Bulgaria" => "BG","Switzerland" => "CH","Cyprus" => "CY","Czech Republic" => "CZ","Germany" => "DE","Denmark" => "DK",
                        "Estonia" => "EE","Spain" => "ES", "Finland" => "FI", "France" => "FR", "Greece" => "GR", "Croatia" => "HR", "Hungary" => "HU", "Ireland" => "IE", "Italy" => "IT", "Lithuania" => "LT",
                        "Luxembourg" => "LU", "Latvia" => "LV", "Montenegro" => "ME", "North Macedonia" => "MK", "Malta" => "MT", "Netherlands" => "NL", "Norway" => "NO", "Poland" => "PL", "Portugal" => "PT",
                        "Romania" => "RO","Serbia" => "RS", "Sweden" => "SE", "Slovenia" => "SI", "Slovakia" => "SK", "Turkey" => "TR", "Ukraine" => "UA", "United Kingdom" => "UK")

seasons = Dict("summer" => [6.0, 7.0, 8.0], "autumn" => [9.0, 10.0, 11.0], "winter" => [12.0, 1.0, 2.0], "spring" => [3.0, 4.0, 5.0])

# filter the data for a specific country and select only the relevant columns
solar_south= dropmissing(solar[(solar.country .== countries[count_south[1]]) .| (solar.country .== countries[count_south[2]]) .| (solar.country .== countries[count_south[3]]) , [:year, :month, :day, :hour, :cf]], disallowmissing=true)
solar_south_annual= dropmissing(solar_south[solar_south.year .>= 1986.0 .&& solar_south.year .<= 2015.0, [:month, :day, :hour, :cf]], disallowmissing=true)
solar_south_summer= dropmissing(solar_south_annual[solar_south_annual.month .>= 6.0 .&& solar_south_annual.month .<= 8.0, [:day, :hour, :cf]], disallowmissing=true)
solar_south_autumn= dropmissing(solar_south_annual[solar_south_annual.month .>= 9.0 .&& solar_south_annual.month .<= 11.0, [:day, :hour, :cf]], disallowmissing=true)
solar_south_spring= dropmissing(solar_south_annual[solar_south_annual.month .>= 3.0 .&& solar_south_annual.month .<= 5.0, [:day, :hour, :cf]], disallowmissing=true)
solar_south_winter= dropmissing(solar_south_annual[(solar_south_annual.month .== 12.0) .| (solar_south_annual.month .== 1.0) .| (solar_south_annual.month .== 2.0), [:day, :hour, :cf]], disallowmissing=true)

boxplot_solar_south_annual_ref = boxplot(solar_south_annual[!,:cf], xlabel="ANNUAL", title="                   SOLAR CF SOUTH EU", whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"
boxplot_solar_south_summer_ref = boxplot(solar_south_summer[!,:cf], xlabel="SUMMER",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="SOLAR CF SOUTH EU SUMMER",
boxplot_solar_south_autumn_ref = boxplot(solar_south_autumn[!,:cf], xlabel="AUTUMN",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="SOLAR CF SOUTH EU AUTUMN",
boxplot_solar_south_spring_ref = boxplot(solar_south_spring[!,:cf], xlabel="SPRING",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="SOLAR CF SOUTH EU SPRING",
boxplot_solar_south_winter_ref = boxplot(solar_south_winter[!,:cf], xlabel="WINTER",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="SOLAR CF SOUTH EU WINTER",


# filter the data for a specific country and select only the relevant columns
solar_north= dropmissing(solar[(solar.country .== countries[count_north[1]]) .| (solar.country .== countries[count_north[2]]) .| (solar.country .== countries[count_north[3]]) , [:year, :month, :day, :hour, :cf]], disallowmissing=true)
solar_north_annual= dropmissing(solar_north[solar_north.year .>= 1986.0 .&& solar_north.year .<= 2015.0, [:month, :day, :hour, :cf]], disallowmissing=true)
solar_north_summer= dropmissing(solar_north_annual[solar_north_annual.month .>= 6.0 .&& solar_north_annual.month .<= 8.0, [:day, :hour, :cf]], disallowmissing=true)
solar_north_autumn= dropmissing(solar_north_annual[solar_north_annual.month .>= 9.0 .&& solar_north_annual.month .<= 11.0, [:day, :hour, :cf]], disallowmissing=true)
solar_north_spring= dropmissing(solar_north_annual[solar_north_annual.month .>= 3.0 .&& solar_north_annual.month .<= 5.0, [:day, :hour, :cf]], disallowmissing=true)
solar_north_winter= dropmissing(solar_north_annual[(solar_north_annual.month .== 12.0) .| (solar_north_annual.month .== 1.0) .| (solar_north_annual.month .== 2.0), [:day, :hour, :cf]], disallowmissing=true)

boxplot_solar_north_annual_ref = boxplot(solar_north_annual[!,:cf], xlabel="ANNUAL", title="                   SOLAR CF NORTH EU", whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"
boxplot_solar_north_summer_ref = boxplot(solar_north_summer[!,:cf], xlabel="SUMMER",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="SOLAR CF NORTH EU SUMMER",
boxplot_solar_north_autumn_ref = boxplot(solar_north_autumn[!,:cf], xlabel="AUTUMN",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="SOLAR CF NORTH EU AUTUMN",
boxplot_solar_north_spring_ref = boxplot(solar_north_spring[!,:cf], xlabel="SPRING",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="SOLAR CF NORTH EU SPRING",
boxplot_solar_north_winter_ref = boxplot(solar_north_winter[!,:cf], xlabel="WINTER",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="SOLAR CF NORTH EU WINTER",

# Combine the box plots into a composite plot
plot(boxplot_solar_south_annual_ref, boxplot_solar_south_summer_ref, boxplot_solar_south_autumn_ref, boxplot_solar_south_spring_ref, boxplot_solar_south_winter_ref, layout=(1,5))
plot(boxplot_solar_north_annual_ref, boxplot_solar_north_summer_ref, boxplot_solar_north_autumn_ref, boxplot_solar_north_spring_ref, boxplot_solar_north_winter_ref, layout=(1,5))




#############################################################################################################################################################
###############################                                      PLOT OF WINDON                                           ###############################
#############################################################################################################################################################

using CSV, DataFrames, Plots, Statistics, Dates


# Read the CSV file into a DataFrame
windon = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDON/PECD-2021.3-country-Onshore-2025.csv", DataFrame)

##### PARAMETERS #####
count = "Norway" # Name of Country
season = "summer"
y1 = 1982.0
y2 = 2000.0
y3 = 2016.0


countries = Dict("Albania" => "AL","Austria" => "AT", "Bosnia and Herzegovina" => "BA","Belgium" => "BE","Bulgaria" => "BG","Switzerland" => "CH","Cyprus" => "CY","Czech Republic" => "CZ","Germany" => "DE","Denmark" => "DK",
                        "Estonia" => "EE","Spain" => "ES", "Finland" => "FI", "France" => "FR", "Greece" => "GR", "Croatia" => "HR", "Hungary" => "HU", "Ireland" => "IE", "Italy" => "IT", "Lithuania" => "LT",
                        "Luxembourg" => "LU", "Latvia" => "LV", "Montenegro" => "ME", "North Macedonia" => "MK", "Netherlands" => "NL", "Norway" => "NO", "Poland" => "PL", "Portugal" => "PT",
                        "Romania" => "RO","Serbia" => "RS", "Sweden" => "SE", "Slovenia" => "SI", "Slovakia" => "SK", "Turkey" => "TR", "Ukraine" => "UA", "United Kingdom" => "UK")

seasons = Dict("summer" => [6.0, 7.0, 8.0], "autumn" => [9.0, 10.0, 11.0], "winter" => [12.0, 1.0, 2.0], "spring" => [3.0, 4.0, 5.0])

# filter the data for a specific country and select only the relevant columns
windon= dropmissing(windon[windon.country .== countries[count], [:year, :month, :day, :hour, :cf]], disallowmissing=true)

# Filter the DataFrame to keep only the summer days of year y
df_subset1 = filter(row -> row.year == y1 && row.month in seasons[season], windon)
df_subset2 = filter(row -> row.year == y2 && row.month in seasons[season], windon)
df_subset3 = filter(row -> row.year == y3 && row.month in seasons[season], windon)

# group the rows by hour
hourly_avg1 = combine(groupby(df_subset1, :hour), :cf => mean => :cf)
hourly_avg2 = combine(groupby(df_subset2, :hour), :cf => mean => :cf)
hourly_avg3 = combine(groupby(df_subset3, :hour), :cf => mean => :cf)

# group the rows by day
daily_avg1 = combine(groupby(df_subset1, [:month, :day]), :cf => mean => :cf)
daily_avg2 = combine(groupby(df_subset2, [:month, :day]), :cf => mean => :cf)
daily_avg3 = combine(groupby(df_subset3, [:month, :day]), :cf => mean => :cf)

# Plot Daily CF dots for one season
locations = [16,46,77]
month_num = seasons[season]
p1 = plot(1:size(daily_avg1, 1), daily_avg1.cf, xlabel="$season months", ylabel="Capacity Factor", title= "WindOn Capacity Factor - $season months $count", xticks=(locations, month_num), label=y1,)
plot!(1:size(daily_avg2, 1), daily_avg2.cf, xlabel="$season months", ylabel="Capacity Factor", title= "WindOn Capacity Factor - $season months $count", label=y2,)
plot!(1:size(daily_avg3, 1), daily_avg3.cf, xlabel="$season months", ylabel="Capacity Factor", title= "WindOn Capacity Factor - $season months $count", label=y3,)

# Plot Hourly CF dots for one day of one season
locations = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
hours = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
p2 = plot(1:size(hourly_avg1, 1), hourly_avg1.cf, xlabel="$season day", ylabel="Capacity Factor", title= "WindOn Capacity Factor - $season day $count", xticks=(locations, hours), label=y1,)
plot!(1:size(hourly_avg2, 1), hourly_avg2.cf, xlabel="$season day", ylabel="Capacity Factor", title= "WindOn Capacity Factor - $season day $count", label=y2,)
plot!(1:size(hourly_avg3, 1), hourly_avg3.cf, xlabel="$season day", ylabel="Capacity Factor", title= "WindOn Capacity Factor - $season day $count", label=y3,)

# Plot Hourly CF dots for one season
locations = [372,1116,1860]
month_num = seasons[season]
p3 = plot(1:size(df_subset1, 1), df_subset1.cf, xlabel="$season months", ylabel="Capacity Factor", title= "WindOn Capacity Factor - $season months $count", xticks=(locations, month_num), label=y1,)
plot!(1:size(df_subset2, 1), df_subset2.cf, xlabel="$season months", ylabel="Capacity Factor", title= "WindOn Capacity Factor - $season months $count", label=y2,)
plot!(1:size(df_subset3, 1), df_subset3.cf, xlabel="$season months", ylabel="Capacity Factor", title= "WindOn Capacity Factor - $season months $count", label=y3,)

plot(p1, p2, layout=(1,2), size=(1600,600))


#############################################################################################################################################################
###############################                                      BOXPLOT WINDON                                           ###############################
#############################################################################################################################################################
using CSV, DataFrames, Plots, Statistics, Dates


# Read the CSV file into a DataFrame
windon = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDON/PECD-2021.3-country-Onshore-2025.csv", DataFrame)

##### PARAMETERS #####
count = ["Portugal", "Spain", "Italy"] # Name of Country 
count_north = ["Norway", "Sweden", "Finland"] 
season = "summer"


countries = Dict("Albania" => "AL","Austria" => "AT", "Bosnia and Herzegovina" => "BA","Belgium" => "BE","Bulgaria" => "BG","Switzerland" => "CH","Cyprus" => "CY","Czech Republic" => "CZ","Germany" => "DE","Denmark" => "DK",
                        "Estonia" => "EE","Spain" => "ES", "Finland" => "FI", "France" => "FR", "Greece" => "GR", "Croatia" => "HR", "Hungary" => "HU", "Ireland" => "IE", "Italy" => "IT", "Lithuania" => "LT",
                        "Luxembourg" => "LU", "Latvia" => "LV", "Montenegro" => "ME", "North Macedonia" => "MK", "Malta" => "MT", "Netherlands" => "NL", "Norway" => "NO", "Poland" => "PL", "Portugal" => "PT",
                        "Romania" => "RO","Serbia" => "RS", "Sweden" => "SE", "Slovenia" => "SI", "Slovakia" => "SK", "Turkey" => "TR", "Ukraine" => "UA", "United Kingdom" => "UK")

seasons = Dict("summer" => [6.0, 7.0, 8.0], "autumn" => [9.0, 10.0, 11.0], "winter" => [12.0, 1.0, 2.0], "spring" => [3.0, 4.0, 5.0])

# filter the data for a specific country and select only the relevant columns
windon_south= dropmissing(windon[(windon.country .== countries[count[1]]) .| (windon.country .== countries[count[2]]) .| (windon.country .== countries[count[3]]) , [:year, :month, :day, :hour, :cf]], disallowmissing=true)
windon_south_annual= dropmissing(windon_south[windon_south.year .>= 1986.0 .&& windon_south.year .<= 2015.0, [:month, :day, :hour, :cf]], disallowmissing=true)
windon_south_summer= dropmissing(windon_south_annual[windon_south_annual.month .>= 6.0 .&& windon_south_annual.month .<= 8.0, [:day, :hour, :cf]], disallowmissing=true)
windon_south_autumn= dropmissing(windon_south_annual[windon_south_annual.month .>= 9.0 .&& windon_south_annual.month .<= 11.0, [:day, :hour, :cf]], disallowmissing=true)
windon_south_spring= dropmissing(windon_south_annual[windon_south_annual.month .>= 3.0 .&& windon_south_annual.month .<= 5.0, [:day, :hour, :cf]], disallowmissing=true)
windon_south_winter= dropmissing(windon_south_annual[(windon_south_annual.month .== 12.0) .| (windon_south_annual.month .== 1.0) .| (windon_south_annual.month .== 2.0), [:day, :hour, :cf]], disallowmissing=true)

windon_south_annual[!,:cf]

boxplot_windon_south_annual_ref = boxplot(windon_south_annual[!,:cf], xlabel="ANNUAL", title="                   WINDON CF SOUTH EU", whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"
boxplot_windon_south_summer_ref = boxplot(windon_south_summer[!,:cf], xlabel="SUMMER",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDON CF SOUTH EU SUMMER",
boxplot_windon_south_autumn_ref = boxplot(windon_south_autumn[!,:cf], xlabel="AUTUMN",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDON CF SOUTH EU AUTUMN",
boxplot_windon_south_spring_ref = boxplot(windon_south_spring[!,:cf], xlabel="SPRING",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDON CF SOUTH EU SPRING",
boxplot_windon_south_winter_ref = boxplot(windon_south_winter[!,:cf], xlabel="WINTER",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDON CF SOUTH EU WINTER",


# filter the data for a specific country and select only the relevant columns
windon_north= dropmissing(windon[(windon.country .== countries[count_north[1]]) .| (windon.country .== countries[count_north[2]]) .| (windon.country .== countries[count_north[3]]) , [:year, :month, :day, :hour, :cf]], disallowmissing=true)
windon_north_annual= dropmissing(windon_north[windon_north.year .>= 1986.0 .&& windon_north.year .<= 2015.0, [:month, :day, :hour, :cf]], disallowmissing=true)
windon_north_summer= dropmissing(windon_north_annual[windon_north_annual.month .>= 6.0 .&& windon_north_annual.month .<= 8.0, [:day, :hour, :cf]], disallowmissing=true)
windon_north_autumn= dropmissing(windon_north_annual[windon_north_annual.month .>= 9.0 .&& windon_north_annual.month .<= 11.0, [:day, :hour, :cf]], disallowmissing=true)
windon_north_spring= dropmissing(windon_north_annual[windon_north_annual.month .>= 3.0 .&& windon_north_annual.month .<= 5.0, [:day, :hour, :cf]], disallowmissing=true)
windon_north_winter= dropmissing(windon_north_annual[(windon_north_annual.month .== 12.0) .| (windon_north_annual.month .== 1.0) .| (windon_north_annual.month .== 2.0), [:day, :hour, :cf]], disallowmissing=true)

boxplot_windon_north_annual_ref = boxplot(windon_north_annual[!,:cf], xlabel="ANNUAL", title="                   WINDON CF NORTH EU", whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"
boxplot_windon_north_summer_ref = boxplot(windon_north_summer[!,:cf], xlabel="SUMMER",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDON CF NORTH EU SUMMER",
boxplot_windon_north_autumn_ref = boxplot(windon_north_autumn[!,:cf], xlabel="AUTUMN",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDON CF NORTH EU AUTUMN",
boxplot_windon_north_spring_ref = boxplot(windon_north_spring[!,:cf], xlabel="SPRING",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDON CF NORTH EU SPRING",
boxplot_windon_north_winter_ref = boxplot(windon_north_winter[!,:cf], xlabel="WINTER",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDON CF NORTH EU WINTER",

# Combine the box plots into a composite plot
plot(boxplot_windon_south_annual_ref, boxplot_windon_south_summer_ref, boxplot_windon_south_autumn_ref, boxplot_windon_south_spring_ref, boxplot_windon_south_winter_ref, layout=(1,5))
plot(boxplot_windon_north_annual_ref, boxplot_windon_north_summer_ref, boxplot_windon_north_autumn_ref, boxplot_windon_north_spring_ref, boxplot_windon_north_winter_ref, layout=(1,5))



#############################################################################################################################################################
###############################                                      PLOT OF WINDOFF                                          ###############################
#############################################################################################################################################################

using CSV, DataFrames, Plots, Statistics, Dates


# Read the CSV file into a DataFrame
windoff = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDOFF/PECD-2021.3-country-Offshore-2025.csv", DataFrame)

##### PARAMETERS #####
count = "Belgium" # Name of Country
season = "winter"
y1 = 1982.0
y2 = 2000.0
y3 = 2016.0


countries = Dict("Belgium" => "BE", "Germany" => "DE", "Denmark" => "DK", "Finland" => "FI", "France" => "FR", "Ireland" => "IE", "Italy" => "IT", "Netherlands" => "NL", "Portugal" => "PT", "Sweden" => "SE", "United Kingdom" => "UK")

seasons = Dict("summer" => [6.0, 7.0, 8.0], "autumn" => [9.0, 10.0, 11.0], "winter" => [12.0, 1.0, 2.0], "spring" => [3.0, 4.0, 5.0])

# filter the data for a specific country and select only the relevant columns
windoff= dropmissing(windoff[windoff.country .== countries[count], [:year, :month, :day, :hour, :cf]], disallowmissing=true)

# Filter the DataFrame to keep only the season days of year y
df_subset1 = filter(row -> row.year == y1 && row.month in seasons[season], windoff)
df_subset2 = filter(row -> row.year == y2 && row.month in seasons[season], windoff)
df_subset3 = filter(row -> row.year == y3 && row.month in seasons[season], windoff)

# group the rows by hour
hourly_avg1 = combine(groupby(df_subset1, :hour), :cf => mean => :cf)
hourly_avg2 = combine(groupby(df_subset2, :hour), :cf => mean => :cf)
hourly_avg3 = combine(groupby(df_subset3, :hour), :cf => mean => :cf)

# group the rows by day
daily_avg1 = combine(groupby(df_subset1, [:month, :day]), :cf => mean => :cf)
daily_avg2 = combine(groupby(df_subset2, [:month, :day]), :cf => mean => :cf)
daily_avg3 = combine(groupby(df_subset3, [:month, :day]), :cf => mean => :cf)

# Plot Daily CF dots for one season
locations = [16,46,77]
month_num = seasons[season]
p1 = plot(1:size(daily_avg1, 1), daily_avg1.cf, xlabel="$season months", ylabel="Capacity Factor", title= "Seasonal Evolution of Windoff Capacity Factor of $count for $season", xticks=(locations, month_num), label=y1,)
plot!(1:size(daily_avg2, 1), daily_avg2.cf, xlabel="$season months", ylabel="Capacity Factor", title= "Seasonal Evolution of Windoff Capacity Factor of $count for $season", label=y2,)
plot!(1:size(daily_avg3, 1), daily_avg3.cf, xlabel="$season months", ylabel="Capacity Factor", title= "Seasonal Evolution of Windoff Capacity Factor of $count for $season", label=y3,)

# Plot Hourly CF dots for one day of one season
locations = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
hours = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
p2 = plot(1:size(hourly_avg1, 1), hourly_avg1.cf, xlabel="$season day", ylabel="Capacity Factor", title= "Daily Evolution of Windoff Capacity Factor of $count for $season", xticks=(locations, hours), label=y1,)
plot!(1:size(hourly_avg2, 1), hourly_avg2.cf, xlabel="$season day", ylabel="Capacity Factor", title= "Daily Evolution of Windoff Capacity Factor of $count for $season", label=y2,)
plot!(1:size(hourly_avg3, 1), hourly_avg3.cf, xlabel="$season day", ylabel="Capacity Factor", title= "Daily Evolution of Windoff Capacity Factor of $count for $season", label=y3,)

# Plot Hourly CF dots for one season
locations = [372,1116,1860]
month_num = seasons[season]
p3 = plot(1:size(df_subset1, 1), df_subset1.cf, xlabel="$season months", ylabel="Capacity Factor", title= "Seasonal Evolution of Windoff Capacity Factor of $count for $season", xticks=(locations, month_num), label=y1,)
plot!(1:size(df_subset2, 1), df_subset2.cf, xlabel="$season months", ylabel="Capacity Factor", title= "Seasonal Evolution of Windoff Capacity Factor of $count for $season", label=y2,)
plot!(1:size(df_subset3, 1), df_subset3.cf, xlabel="$season months", ylabel="Capacity Factor", title= "Seasonal Evolution of Windoff Capacity Factor of $count for $season", label=y3,)

plot(p1, p2, layout=(1,2), size=(1600,600))


#############################################################################################################################################################
###############################                                      BOXPLOT WINDOFF                                          ###############################
#############################################################################################################################################################
using CSV, DataFrames, Plots, Statistics, Dates


# Read the CSV file into a DataFrame
windoff = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDOFF/PECD-2021.3-country-Offshore-2025.csv", DataFrame)

##### PARAMETERS #####
count = ["Portugal", "Spain", "Italy"] # Name of Country 
count_north = ["Norway", "Sweden", "Finland"] 
season = "summer"


countries = Dict("Albania" => "AL","Austria" => "AT", "Bosnia and Herzegovina" => "BA","Belgium" => "BE","Bulgaria" => "BG","Switzerland" => "CH","Cyprus" => "CY","Czech Republic" => "CZ","Germany" => "DE","Denmark" => "DK",
                        "Estonia" => "EE","Spain" => "ES", "Finland" => "FI", "France" => "FR", "Greece" => "GR", "Croatia" => "HR", "Hungary" => "HU", "Ireland" => "IE", "Italy" => "IT", "Lithuania" => "LT",
                        "Luxembourg" => "LU", "Latvia" => "LV", "Montenegro" => "ME", "North Macedonia" => "MK", "Malta" => "MT", "Netherlands" => "NL", "Norway" => "NO", "Poland" => "PL", "Portugal" => "PT",
                        "Romania" => "RO","Serbia" => "RS", "Sweden" => "SE", "Slovenia" => "SI", "Slovakia" => "SK", "Turkey" => "TR", "Ukraine" => "UA", "United Kingdom" => "UK")

seasons = Dict("summer" => [6.0, 7.0, 8.0], "autumn" => [9.0, 10.0, 11.0], "winter" => [12.0, 1.0, 2.0], "spring" => [3.0, 4.0, 5.0])

# filter the data for a specific country and select only the relevant columns
windoff_south= dropmissing(windoff[(windoff.country .== countries[count[1]]) .| (windoff.country .== countries[count[2]]) .| (windoff.country .== countries[count[3]]) , [:year, :month, :day, :hour, :cf]], disallowmissing=true)
windoff_south_annual= dropmissing(windoff_south[windoff_south.year .>= 1986.0 .&& windoff_south.year .<= 2015.0, [:month, :day, :hour, :cf]], disallowmissing=true)
windoff_south_summer= dropmissing(windoff_south_annual[windoff_south_annual.month .>= 6.0 .&& windoff_south_annual.month .<= 8.0, [:day, :hour, :cf]], disallowmissing=true)
windoff_south_autumn= dropmissing(windoff_south_annual[windoff_south_annual.month .>= 9.0 .&& windoff_south_annual.month .<= 11.0, [:day, :hour, :cf]], disallowmissing=true)
windoff_south_spring= dropmissing(windoff_south_annual[windoff_south_annual.month .>= 3.0 .&& windoff_south_annual.month .<= 5.0, [:day, :hour, :cf]], disallowmissing=true)
windoff_south_winter= dropmissing(windoff_south_annual[(windoff_south_annual.month .== 12.0) .| (windoff_south_annual.month .== 1.0) .| (windoff_south_annual.month .== 2.0), [:day, :hour, :cf]], disallowmissing=true)

boxplot_windoff_south_annual_ref = boxplot(windoff_south_annual[!,:cf], xlabel="ANNUAL", title="                   WINDOFF CF SOUTH EU", whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"
boxplot_windoff_south_summer_ref = boxplot(windoff_south_summer[!,:cf], xlabel="SUMMER",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDOFF CF SOUTH EU SUMMER",
boxplot_windoff_south_autumn_ref = boxplot(windoff_south_autumn[!,:cf], xlabel="AUTUMN",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDOFF CF SOUTH EU AUTUMN",
boxplot_windoff_south_spring_ref = boxplot(windoff_south_spring[!,:cf], xlabel="SPRING",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDOFF CF SOUTH EU SPRING",
boxplot_windoff_south_winter_ref = boxplot(windoff_south_winter[!,:cf], xlabel="WINTER",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDOFF CF SOUTH EU WINTER",


# filter the data for a specific country and select only the relevant columns
windoff_north= dropmissing(windoff[(windoff.country .== countries[count_north[1]]) .| (windoff.country .== countries[count_north[2]]) .| (windoff.country .== countries[count_north[3]]) , [:year, :month, :day, :hour, :cf]], disallowmissing=true)
windoff_north_annual= dropmissing(windoff_north[windoff_north.year .>= 1986.0 .&& windoff_north.year .<= 2015.0, [:month, :day, :hour, :cf]], disallowmissing=true)
windoff_north_summer= dropmissing(windoff_north_annual[windoff_north_annual.month .>= 6.0 .&& windoff_north_annual.month .<= 8.0, [:day, :hour, :cf]], disallowmissing=true)
windoff_north_autumn= dropmissing(windoff_north_annual[windoff_north_annual.month .>= 9.0 .&& windoff_north_annual.month .<= 11.0, [:day, :hour, :cf]], disallowmissing=true)
windoff_north_spring= dropmissing(windoff_north_annual[windoff_north_annual.month .>= 3.0 .&& windoff_north_annual.month .<= 5.0, [:day, :hour, :cf]], disallowmissing=true)
windoff_north_winter= dropmissing(windoff_north_annual[(windoff_north_annual.month .== 12.0) .| (windoff_north_annual.month .== 1.0) .| (windoff_north_annual.month .== 2.0), [:day, :hour, :cf]], disallowmissing=true)

boxplot_windoff_north_annual_ref = boxplot(windoff_north_annual[!,:cf], xlabel="ANNUAL", title="                   WINDOFF CF NORTH EU", whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"
boxplot_windoff_north_summer_ref = boxplot(windoff_north_summer[!,:cf], xlabel="SUMMER",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDOFF CF NORTH EU SUMMER",
boxplot_windoff_north_autumn_ref = boxplot(windoff_north_autumn[!,:cf], xlabel="AUTUMN",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDOFF CF NORTH EU AUTUMN",
boxplot_windoff_north_spring_ref = boxplot(windoff_north_spring[!,:cf], xlabel="SPRING",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"    title="WINDOFF CF NORTH EU SPRING",
boxplot_windoff_north_winter_ref = boxplot(windoff_north_winter[!,:cf], xlabel="WINTER",  whisker_width=0.5, marker_size=3, label = "REF") # ylabel="Variable Value"   

# Combine the box plots into a composite plot
plot(boxplot_windoff_south_annual_ref, boxplot_windoff_south_summer_ref, boxplot_windoff_south_autumn_ref, boxplot_windoff_south_spring_ref, boxplot_windoff_south_winter_ref, layout=(1,5))
plot(boxplot_windoff_north_annual_ref, boxplot_windoff_north_summer_ref, boxplot_windoff_north_autumn_ref, boxplot_windoff_north_spring_ref, boxplot_windoff_north_winter_ref, layout=(1,5))



#################################################################################################################################################################################################
#####################################                                                           FINAL BOXPLOTS                                              ##################################### 
#################################################################################################################################################################################################

### SOUTH EU ###
# Solar
plot(boxplot_solar_south_annual_ref, boxplot_solar_south_annual_coper, boxplot_solar_south_summer_ref, boxplot_solar_south_summer_coper, boxplot_solar_south_autumn_ref, boxplot_solar_south_autumn_coper, boxplot_solar_south_spring_ref, boxplot_solar_south_spring_coper, boxplot_solar_south_winter_ref,boxplot_solar_south_winter_coper, layout=(1,10), size=(1200,500))

# Windon
plot(boxplot_windon_south_annual_ref, boxplot_windon_south_annual_coper, boxplot_windon_south_summer_ref, boxplot_windon_south_summer_coper, boxplot_windon_south_autumn_ref, boxplot_windon_south_autumn_coper, boxplot_windon_south_spring_ref, boxplot_windon_south_spring_coper, boxplot_windon_south_winter_ref,boxplot_windon_south_winter_coper, layout=(1,10), size=(1200,500))

# Windoff
plot(boxplot_windoff_south_annual_ref, boxplot_windoff_south_annual_coper, boxplot_windoff_south_summer_ref, boxplot_windoff_south_summer_coper, boxplot_windoff_south_autumn_ref, boxplot_windoff_south_autumn_coper, boxplot_windoff_south_spring_ref, boxplot_windoff_south_spring_coper, boxplot_windoff_south_winter_ref,boxplot_windoff_south_winter_coper, layout=(1,10), size=(1200,500))


### NORTH EU ###
# Solar
plot(boxplot_solar_north_annual_ref, boxplot_solar_north_annual_coper, boxplot_solar_north_summer_ref, boxplot_solar_north_summer_coper, boxplot_solar_north_autumn_ref, boxplot_solar_north_autumn_coper, boxplot_solar_north_spring_ref, boxplot_solar_north_spring_coper, boxplot_solar_north_winter_ref,boxplot_solar_north_winter_coper, layout=(1,10), size=(1200,500))

# Windon
plot(boxplot_windon_north_annual_ref, boxplot_windon_north_annual_coper, boxplot_windon_north_summer_ref, boxplot_windon_north_summer_coper, boxplot_windon_north_autumn_ref, boxplot_windon_north_autumn_coper, boxplot_windon_north_spring_ref, boxplot_windon_north_spring_coper, boxplot_windon_north_winter_ref,boxplot_windon_north_winter_coper, layout=(1,10), size=(1200,500))

# Windoff
plot(boxplot_windoff_north_annual_ref, boxplot_windoff_north_annual_coper, boxplot_windoff_north_summer_ref, boxplot_windoff_north_summer_coper, boxplot_windoff_north_autumn_ref, boxplot_windoff_north_autumn_coper, boxplot_windoff_north_spring_ref, boxplot_windoff_north_spring_coper, boxplot_windoff_north_winter_ref,boxplot_windoff_north_winter_coper, layout=(1,10), size=(1200,500))

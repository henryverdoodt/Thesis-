## Master Thesis: Analyse the impact of climate change on a G&TEP of EU Power System
# Author: Henry Verdoodt
# Last update: April 20, 2023

#= 
## Step 0: Activate environment - ensure consistency accross computers
using Pkg
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
Pkg.add("Images") 
 =#

 using CSV, YAML, JuMP, DataFrames, Distributions, Gurobi, Images, Plots, PolyChaos, InteractiveUtils
 using StatsPlots, Images

 # Start the timer
 t_start = time()
 
 # Include Files
 include("/Users/henryverdoodt/Documents/CODE/JULIA/G&TEP/REFORMAT_DATA.jl")
 include("/Users/henryverdoodt/Documents/CODE/JULIA/G&TEP/PLOT_MODEL.jl")
 
 
 ############################ DATA PARAMETERS ############################
 countries = ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO", "CH", "FI", "IE", "IT", "AT", "PT", "SE"]   #["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO"]

 # demand_entso: 1982.0 - 2016.0    and   solar,windon,windoff_entso: 1982.0 - 2019.0
 year = 1986  # year weather data
 representative_years = [1995.0, 2008.0, 2009.0] # year demand 
 weights = [0.233, 0.367, 0.4] # weights demand years
 aggregate_3h = true # aggregate 1h timestep data to 3h timestep data
 
 # Countries where demand/solar/windon/windoff data exist
 countries_demand = countries_demand_entso       
 countries_solar = countries_solar_coper_CNRM         # countries_solar_coper_cnrm        # countries_solar_entso 
 countries_windon = countries_windon_coper_CNRM       # countries_windon_coper_cnrm       # countries_windon_entso 
 countries_windoff = countries_windoff_coper_CNRM     # countries_windoff_coper_cnrm      # countries_windoff_entso
 
 #Countries from domain where demand data exist (not IT and PT)
 countries_dem = ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO", "CH", "FI", "IE", "AT", "SE"]
 
 # Read the CSV file into a DataFrame
 
 ## REF CASE DEMAND
 dem = reformat_entso_demand(demand_entso, countries_dem, countries_demand_entso, representative_years, weights, aggregate_3h)
 dem_PT_IT = reformat_demand_PT_IT(new_demand_PT_IT, aggregate_3h)
 dem = hcat(dem, dem_PT_IT)

 #=
 ## REF CASE CF
 sol = reformat_entso_solar(solar_entso, countries, countries_solar_entso, representative_years, weights, aggregate_3h)
 won = reformat_entso_windon(windon_entso, countries, countries_windon_entso, representative_years, weights, aggregate_3h) 
 woff = reformat_entso_windoff(windoff_entso, countries, countries_windoff_entso, representative_years, weights, aggregate_3h)
 =#
 
 # PROJECTIONS CF
 sol = reformat_coper_solar(solar_coper_CNRM, countries, countries_solar_coper_CNRM, year)          # reformat_coper_solar(solar_coper_cnrm, countries, countries_solar_coper_cnrm, year)           # reformat_entso_solar(solar_entso, countries, countries_solar_entso, y, aggregate_3h)
 won = reformat_coper_windon(windon_coper_CNRM, countries, countries_windon_coper_CNRM, year)       # reformat_coper_windon(windon_coper_cnrm, countries, countries_windon_coper_cnrm, year)        # reformat_entso_windon(windon_entso, countries, countries_windon_entso, y, aggregate_3h)
 woff = reformat_coper_windoff(windoff_coper_CNRM, countries, countries_windoff_coper_CNRM, year)   #reformat_coper_windoff(windoff_coper_cnrm, countries, countries_windoff_coper_cnrm, year)      # reformat_entso_windoff(windoff_entso, countries, countries_windoff_entso, y, aggregate_3h)
 
 
 # DATA GENERATION AND TRANSMISSION
 data = YAML.load_file(joinpath(@__DIR__, "overview_data.yaml"))
 network = YAML.load_file(joinpath(@__DIR__, "network_west_europe.yaml"))
 
 println("DATA Done")
 
 
 ## Step 2: create model & pass data to model
 m = Model(optimizer_with_attributes(Gurobi.Optimizer))
 
 
 function timestep(aggregate::Bool)
     if aggregate
         return 2920
     else
         return 8760
     end
 end
 
 # Step 2a: create sets
 function define_sets!(m::Model, data::Dict, network::Dict)
     m.ext[:sets] = Dict()
     J = m.ext[:sets][:J] = 1:timestep(aggregate_3h) # Timesteps: 2920  or  8760
     I = m.ext[:sets][:I] = [id for id in keys(data["PowerSector"]) if id ∉ ["CCGT_old", "Biomass", "CCGT_new"]] # Generators per type 
     L_ac = m.ext[:sets][:AC_Lines] = [network["AC_Lines"]["AC_$(i)"]["Connection"] for i in 1:length(keys(network["AC_Lines"]))] # AC Lines
     L_dc = m.ext[:sets][:DC_Lines] = [network["DC_Lines"]["DC_$(i)"]["Connection"] for i in 1:length(keys(network["DC_Lines"]))] # DC Lines
     L = m.ext[:sets][:Lines] = union(m.ext[:sets][:AC_Lines],m.ext[:sets][:DC_Lines]) # Lines
     N = m.ext[:sets][:Nodes] =  [network["Nodes"]["Node_$(i)"] for i in 1:length(keys(network["Nodes"]))]# Nodes
 
     #S = m.ext[:sets][:S] = [s for s in keys(data["Storage_H2"])] # Storage per type
     S = m.ext[:sets][:S] = [s for s in keys(data["Storage"])] # Storage per type
     return m # return model
 end
 
 
 # Step 2b: add time series
 function process_time_series_data!(m::Model, dem::DataFrame, sol::DataFrame, won::DataFrame, woff::DataFrame)
     # extract the relevant sets
     I = m.ext[:sets][:I]
     J = m.ext[:sets][:J]
     N = m.ext[:sets][:Nodes]
     S = m.ext[:sets][:S]
 
     # create dictionary to store time series
     m.ext[:timeseries] = Dict()
     m.ext[:timeseries][:D] = Dict()
     m.ext[:timeseries][:AFS] = Dict()
     m.ext[:timeseries][:AFWON] = Dict()
     m.ext[:timeseries][:AFWOFF] = Dict()
 
     # TimeSeries
     for c in intersect(countries, countries_demand)
         m.ext[:timeseries][:D][Symbol(c)] = dem[!, Symbol(c)]
     end
     for c in intersect(countries, countries_solar)
         m.ext[:timeseries][:AFS][Symbol(c)] = sol[!, Symbol(c)]
     end
     for c in intersect(countries, countries_windon)
         m.ext[:timeseries][:AFWON][Symbol(c)] = won[!, Symbol(c)]
     end
     for c in intersect(countries, countries_windoff)
         m.ext[:timeseries][:AFWOFF][Symbol(c)] = woff[!, Symbol(c)]
     end
 
     # return model
     return m
 end
 
 # step 2c: process input parameters
 function process_parameters!(m::Model, data::Dict, network::Dict)
     # extract sets
     I = m.ext[:sets][:I]
     S = m.ext[:sets][:S]
     N = m.ext[:sets][:Nodes]
     IED = ["Nuclear", "SPP_lignite", "SPP_coal", "OCGT", "ICE", "Biofuel", "Hydro_RoR"]
 
     # Create parameter dictonary
     m.ext[:parameters] = Dict()
 
     # general parameters
     disc_r = m.ext[:parameters][:disc_r] = data["General"]["discount_rate"] # -, discount rate
     VOLL = m.ext[:parameters][:VOLL] = data["General"]["valueOfLostLoad"] # EUR/MWh, VOLL
     hurdle_cost = m.ext[:parameters][:hurdle_cost] = data["General"]["hurdle_cost"] # EUR/MWh, hurdle_cost
     alphaCO2 = m.ext[:parameters][:alphaCO2] = data["ETS"]["P_calibration"] # EUR/ton, alphaCO2
     #Bl_ac = m.ext[:parameters][:Bl_ac] = network["Bl_ac"]*10^(-6) # S/km, Susceptance of AC line 
     #Bl_dc =  m.ext[:parameters][:Bl_dc] = network["Bl_dc"]*10^(-6) # S/km, Susceptance of DC line  
     #Bl = m.ext[:parameters][:Bl] = cat(Bl_ac, Bl_dc, dims=1) # S/km, Susceptance of All line [USE cat() OR union() function??]
     cap_exist = m.ext[:parameters][:cap_exist] = Dict(n => Dict(i => network["Existing_Capacity"][n][i] for i in IED) for n in N)
     Emax_Storage = m.ext[:parameters][:Emax_Storage] = Dict(n => Dict(i => network["Storage_Existing_Capacity"][n][i]["Emax"] for i in S) for n in N)
     Pch_Storage = m.ext[:parameters][:Pch_Storage] = Dict(n => Dict(i => network["Storage_Existing_Capacity"][n][i]["P_ch"] for i in S) for n in N)
     Pdis_Storage = m.ext[:parameters][:Pdis_Storage] = Dict(n => Dict(i => network["Storage_Existing_Capacity"][n][i]["P_dis"] for i in S) for n in N)


     # parameters of generators per unit
     d = data["PowerSector"]
     #betaD = m.ext[:parameters][:betaD] = Dict(i => (d[SubString(i,1:length(i))]["fuelCosts"]/d[SubString(i,1:length(i))]["efficiency"]) for i in I)  # EUR/MWh, Cost of generation of unit i
     #deltaD = m.ext[:parameters][:deltaD] = Dict(i => (d[SubString(i,1:length(i))]["emissions"]/d[SubString(i,1:length(i))]["efficiency"]) for i in I) # ton/MWh, Emissions of generation of unit i
     #GMAX = m.ext[:parameters][:GMAX] = Dict(i => d[SubString(i,1:length(i))]["gMAX"] for i in I) # MW, Maximum Power Output of one generation unit
     VC = m.ext[:parameters][:VC] = Dict(i => ((d[SubString(i,1:length(i))]["fuelCosts"] + d[SubString(i,1:length(i))]["emissions"]*alphaCO2)/d[SubString(i,1:length(i))]["efficiency"]) for i in I) # EUR/MWh, Variable Cost of unit i
     #OC = m.ext[:parameters][:OC] = Dict(i => 10^3*d[SubString(i,1:length(i))]["OC"] for i in I) # EUR/MW, Overnight investment cost of a unit
     #LifeT = m.ext[:parameters][:LifeT] = Dict(i => d[SubString(i,1:length(i))]["Lifetime"] for i in I) # years, Lifetime of Asset
     IC = m.ext[:parameters][:IC] = Dict(i => 10^3*(d[SubString(i,1:length(i))]["OC"]*disc_r)/(1-(1/(1+disc_r)^(d[SubString(i,1:length(i))]["Lifetime"]))) for i in I) # EUR/MWy, Investment Cost in generation unit i
     AF = m.ext[:parameters][:AF] = Dict(i => d[SubString(i,1:length(i))]["AF"] for i in I) # Availability facor of a generator
 
    #  # parameters of Storage per unit
    #  dS = data["Storage_H2"]
    #  VCS = m.ext[:parameters][:VCS] = Dict(i => (dS[SubString(i,1:length(i))]["fuelCosts"]/dS[SubString(i,1:length(i))]["efficiency"]) for i in S)  # EUR/MWh, Cost of storage of unit i
    #  OCS = m.ext[:parameters][:OCS] = Dict(i => 10^3*dS[SubString(i,1:length(i))]["OC"] for i in S) # EUR/MW, Overnight investment cost of a unit
    #  LifeTS = m.ext[:parameters][:LifeTS] = Dict(i => dS[SubString(i,1:length(i))]["Lifetime"] for i in S) # years, Lifetime of Asset
    #  ICS = m.ext[:parameters][:ICS] = Dict(i => 10^3*(dS[SubString(i,1:length(i))]["OC"]*disc_r)/(1-(1/(1+disc_r)^(dS[SubString(i,1:length(i))]["Lifetime"]))) for i in S) # EUR/MWy, Investment Cost in generation unit i
    #  AFSt = m.ext[:parameters][:AFSt] = Dict(i => dS[SubString(i,1:length(i))]["AF"] for i in S) # Availability facor of Storage technology
     
     # parameters of Storage per unit
     dS = data["Storage"]
     VCS = m.ext[:parameters][:VCS] = Dict(i => (dS[SubString(i,1:length(i))]["VC"]) for i in S)  # EUR/MWh, Variable Cost of storage of unit i
     AFSt = m.ext[:parameters][:AFSt] = Dict(i => dS[SubString(i,1:length(i))]["AF"] for i in S) # Availability facor of Storage technology
     eff_ch = m.ext[:parameters][:eff_ch] = Dict(i => dS[SubString(i,1:length(i))]["eff_ch"] for i in S) # Availability facor of Storage technology
     eff_dis = m.ext[:parameters][:eff_dis] = Dict(i => dS[SubString(i,1:length(i))]["eff_dis"] for i in S) # Availability facor of Storage technology


     #FC = m.ext[:parameters][:FC] = Dict(i => ((d[SubString(i,1:length(i))]["fuelCosts"])/d[SubString(i,1:length(i))]["efficiency"]) for i in I) # EUR/MWh, Fuelcost Cost of unit i
     #CO2 = m.ext[:parameters][:CO2] = Dict(i => ((d[SubString(i,1:length(i))]["emissions"])/d[SubString(i,1:length(i))]["efficiency"]) for i in I) # ton/MWh, Variable Cost of unit i
 
     # parameters of transmission lines AC
     Fl_MAX_AC = m.ext[:parameters][:Fl_MAX_AC] = [network["AC_Lines"]["AC_$(i)"]["Capacity_2030"] for i in 1:length(keys(network["AC_Lines"]))] # MW, AC Lines
     #L_length_AC = m.ext[:parameters][:L_length_AC] = [network["AC_Lines"]["AC_$(i)"]["Length"] for i in 1:length(keys(network["AC_Lines"]))] # km, Length of AC transmission line
     #L_price_AC = m.ext[:parameters][:L_price_AC] = [network["AC_Lines"]["AC_$(i)"]["Price"] for i in 1:length(keys(network["AC_Lines"]))] # EUR/km.MW, Price of AC line per unit length
     #OC_var_AC = m.ext[:parameters][:OC_var_AC] = [network["AC_Lines"]["AC_$(i)"]["Price"] * network["AC_Lines"]["AC_$(i)"]["Length"] for i in 1:length(keys(network["AC_Lines"]))] # EUR/MW, Investment Cost for AC Transmission line
     #IC_var_AC = m.ext[:parameters][:IC_var_AC] = [(network["AC_Lines"]["AC_$(i)"]["Price"] * network["AC_Lines"]["AC_$(i)"]["Length"]*disc_r)/(1-(1/(1+disc_r)^(network["Lifetime_ac"]))) for i in 1:length(keys(network["AC_Lines"]))] # EUR/MWy, Investment Cost of AC line
 
     # parameters of transmission lines DC 
     Fl_MAX_DC = m.ext[:parameters][:Fl_MAX_DC] = [network["DC_Lines"]["DC_$(i)"]["Capacity_2030"] for i in 1:length(keys(network["DC_Lines"]))] # MW, DC Lines
     #L_length_DC = m.ext[:parameters][:L_length_DC] = [network["DC_Lines"]["DC_$(i)"]["Length"] for i in 1:length(keys(network["DC_Lines"]))] # km, Length of DC transmission line
     #L_price_DC = m.ext[:parameters][:L_price_DC] = [network["DC_Lines"]["DC_$(i)"]["Price"] for i in 1:length(keys(network["DC_Lines"]))] # EUR/km.MW, Price of DC line per unit length
     #OC_var_DC = m.ext[:parameters][:OC_var_DC] = [network["DC_Lines"]["DC_$(i)"]["Price"] * network["DC_Lines"]["DC_$(i)"]["Length"] for i in 1:length(keys(network["DC_Lines"]))] # EUR/MW, Investment Cost for DC Transmission line
     #IC_var_DC = m.ext[:parameters][:IC_var_DC] = [(network["DC_Lines"]["DC_$(i)"]["Price"] * network["DC_Lines"]["DC_$(i)"]["Length"]*disc_r)/(1-(1/(1+disc_r)^(network["Lifetime_dc"]))) for i in 1:length(keys(network["DC_Lines"]))] # EUR/MWy, Investment Cost of DC line
     Fl_MAX = m.ext[:parameters][:Fl_MAX] = cat(Fl_MAX_AC, Fl_MAX_DC, dims=1) # MW, All Lines, Maximum capacity [USE cat() OR union() function??]
     # return model
     return m
 end
 
 # call functions
 define_sets!(m, data, network)  
 println("SETS Done")
 process_time_series_data!(m, dem, sol, won, woff)   # process_time_series_data!(m, dem_3h, sol_3h, won_3h, woff_3h)  # process_time_series_data!(m, dem, sol, won, woff) 
 println("TIME SERIES Done")
 process_parameters!(m, data, network)
 println("PARAMETERS Done")
 
 function find_line_number(network::Dict, countries::Vector)
     for (name, line) in network
         if line["Connection"] == countries
             # Return the line number if the cities match
             return parse(Int, split(name, "_")[2])
         end
     end
     # Return Nothing if no match is found
     return nothing
 end
 
 
 ## Step 3: construct your model
 function build_model!(m::Model)
     # Create m.ext entries "variables", "expressions" and "constraints"
     m.ext[:variables] = Dict()
     m.ext[:expressions] = Dict()
     m.ext[:constraints] = Dict()
 
     # Extract sets
     J = m.ext[:sets][:J] # Timesteps
     I = m.ext[:sets][:I] # Generators per type
     I_CO2 = ["CCGT_new", "OCGT", "ICE"] # Generators units emitting CO2
     ID = ["Nuclear", "SPP_lignite", "SPP_coal", "OCGT", "ICE", "Biofuel", "Hydro_RoR"]  # ID = ["Nuclear", "SPP_lignite", "SPP_coal", "CCGT_old", "CCGT_new", "OCGT", "ICE", "Biomass"]; As we assume Greenfield method we only take the new tech. and dont take coal into account
     IV = ["Solar", "WindOnshore", "WindOffshore"]
     ITOT = ["Solar", "WindOnshore", "WindOffshore", "Nuclear", "SPP_lignite", "SPP_coal", "OCGT", "ICE", "Biofuel", "Hydro_RoR"]
     L_ac = m.ext[:sets][:AC_Lines] # AC Lines
     L_dc = m.ext[:sets][:DC_Lines] # DC Lines
     L = m.ext[:sets][:Lines] # All Lines
     N = m.ext[:sets][:Nodes] # Nodes
     S = m.ext[:sets][:S] # Storage
 
 
     # Extract time series data
     D = m.ext[:timeseries][:D]
     AFWON = m.ext[:timeseries][:AFWON]
     AFWOFF = m.ext[:timeseries][:AFWOFF]
     AFS = m.ext[:timeseries][:AFS]
 
     # Extract parameters
         ### general parameters
     VOLL = m.ext[:parameters][:VOLL] 
     hurdle_cost = m.ext[:parameters][:hurdle_cost] 
     #Bl_ac = m.ext[:parameters][:Bl_ac] 
     #Bl_dc =  m.ext[:parameters][:Bl_dc] 
     #Bl = m.ext[:parameters][:Bl]
     cap_exist = m.ext[:parameters][:cap_exist]
     
         ### parameters of generators per unit
     VC = m.ext[:parameters][:VC]
     IC = m.ext[:parameters][:IC]
     AF = m.ext[:parameters][:AF]
 
     #FC = m.ext[:parameters][:FC]
     #CO2 = m.ext[:parameters][:CO2]
 
         ### parameters of storage per unit
     VCS = m.ext[:parameters][:VCS]
     #ICS = m.ext[:parameters][:ICS]
     AFSt = m.ext[:parameters][:AFSt]
     Pch_Storage = m.ext[:parameters][:Pch_Storage]
     Pdis_Storage = m.ext[:parameters][:Pdis_Storage]
     Emax_Storage = m.ext[:parameters][:Emax_Storage]
     eff_ch = m.ext[:parameters][:eff_ch]
     eff_dis = m.ext[:parameters][:eff_dis]
     
 
         ### parameters of transmission lines
     #IC_var_AC = m.ext[:parameters][:IC_var_AC]
     #IC_var_DC = m.ext[:parameters][:IC_var_DC] 
     Fl_MAX_AC = m.ext[:parameters][:Fl_MAX_AC]
     Fl_MAX_DC = m.ext[:parameters][:Fl_MAX_DC]
    
 
     # create variables 
     g = m.ext[:variables][:g] = @variable(m, [i=I,j=J,n=N], lower_bound=0, base_name="generation") # Power produced by candidate generating unit i at time j in node n [MW]
     cap = m.ext[:variables][:cap] = @variable(m, [i=I,n=N], lower_bound=0, base_name="capacity of generating unit") # Capacity of candidate generating unit i in node n[MW]
     ens = m.ext[:variables][:ens] = @variable(m, [j=J,n=N], lower_bound=0, base_name="energy not served") # Energy Not Served at time j in node n [MWh] (OR Should I use Load shed of demand instead in MW?)
     pl = m.ext[:variables][:pl] = @variable(m, [j=J,l=L], lower_bound=0, base_name="power flow in transmission") # Power flow through Transmission Line l at time j [MW]
     #θ = m.ext[:variables][:θ] = @variable(m, [j=J,n=N], lower_bound=0, base_name="voltage angle") # Voltage angle at node n and time j [rad]
     #varlac = m.ext[:variables][:varlac] = @variable(m, [la=L_ac], lower_bound=0, base_name="capacity of ac line") # Capacity of AC line [MW]
     #varldc = m.ext[:variables][:varldc] = @variable(m, [ld=L_dc], lower_bound=0, base_name="capacity of dc line") # Capacity of DC line [MW]
 
    #  # variable for storage
    #  g_PtH = m.ext[:variables][:g_PtH] = @variable(m, [j=J,n=N], lower_bound=0, base_name="Power to Hydrogen") # Power used to produce Hydrogen at time j in node n [MW]
    #  g_OCGT = m.ext[:variables][:g_OCGT] = @variable(m, [j=J,n=N], lower_bound=0, base_name="Hydrogen to Power") # Hydrogen to generate Power at time j in node n [MW]
    #  cap_PtH = m.ext[:variables][:cap_PtH] = @variable(m, [n=N], lower_bound=0, base_name="capacity of PtH unit") # Capacity of PtH unit in node n[MW]
    #  cap_OCGT = m.ext[:variables][:cap_OCGT] = @variable(m, [n=N], lower_bound=0, base_name="capacity of OCGT unit") # Capacity of OCGT unit in node n[MW]
 
     # variable for storage
     c = m.ext[:variables][:c] = @variable(m, [i=S,j=J,n=N], lower_bound=0, upper_bound= Pch_Storage[n][i], base_name="charging")
     d = m.ext[:variables][:d] = @variable(m, [i=S,j=J,n=N], lower_bound=0, upper_bound= Pdis_Storage[n][i], base_name="discharging")
     e = m.ext[:variables][:e] = @variable(m, [i=S,j=J,n=N], lower_bound=0, upper_bound= Emax_Storage[n][i], base_name="SOC")

     #alphaCO2 = m.ext[:variables][:alphaCO2] = @variable(m, lower_bound=0, base_name="carbon price") # Carbon price [EUR]
 
     # Objective
     obj = m.ext[:objective] = @objective(m, Min, sum(IC[i]*cap[i,n] for i in IV, n in N) + sum(VC[i]*g[i,j,n] for i in ITOT, j in J, n in N) + sum(VCS[i]*(d[i,j,n]+c[i,j,n]) for i in S, j in J, n in N) + sum(VOLL*ens[j,n] for j in J, n in N)) + sum(hurdle_cost*pl[j,l] for j in J, l in L)
     
     #obj = m.ext[:objective] = @objective(m, Min, sum(IC[i]*cap[i,n] for i in IV, n in N) + sum(VC[i]*g[i,j,n] for i in ITOT, j in J, n in N) + sum(VOLL*ens[j,n] for j in J, n in N)) + sum(hurdle_cost*pl[j,l] for j in J, l in L)
     #obj = m.ext[:objective] = @objective(m, Min, sum(IC[i]*cap[i,n] for i in I, n in N) + sum(ICS["PtH"]*cap_PtH[n] for n in N) + sum(ICS["OCGT_H"]*cap_OCGT[n] for n in N) + sum(VC[i]*g[i,j,n] for i in I, j in J, n in N) + sum(VCS["PtH"]*g_PtH[j,n] for j in J, n in N) + sum(VCS["OCGT_H"]*g_OCGT[j,n] for j in J, n in N) + sum(VOLL*ens[j,n] for j in J, n in N)) + sum(hurdle_cost*pl[j,l] for j in J, l in L)
     #obj = m.ext[:objective] = @objective(m, Min, sum(IC[i]*cap[i,n] for i in I, n in N) + sum(ICS["PtH"]*cap_PtH[n] for n in N) + sum(ICS["OCGT_H"]*cap_OCGT[n] for n in N) + sum(IC_var_AC[find_line_number(network["AC_Lines"], la)]*varlac[la] for la in L_ac) + sum(IC_var_DC[find_line_number(network["DC_Lines"], ld)]*varldc[ld] for ld in L_dc) + sum(VC[i]*g[i,j,n] for i in I, j in J, n in N) + sum(VCS["PtH"]*g_PtH[j,n] for j in J, n in N) + sum(VCS["OCGT_H"]*g_OCGT[j,n] for j in J, n in N) + sum(VOLL*ens[j,n] for j in J, n in N))
     #obj = m.ext[:objective] = @objective(m, Min, sum(IC[i]*cap[i,n] for i in I, n in N)  + sum(IC_var_AC[find_line_number(network["AC_Lines"], la)]*varlac[la] for la in L_ac) + sum(IC_var_DC[find_line_number(network["DC_Lines"], ld)]*varldc[ld] for ld in L_dc) + sum((FC[i]*g[i,j,n] + CO2[i]*alphaCO2*g[i,j,n]) for i in I, j in J, n in N) + sum(VOLL*ens[j,n] for j in J, n in N))
 
      
     # Constraints
     con_MC = m.ext[:constraints][:con_MC] = @constraint(m, [j=J, n=N], sum(g[i,j,n] for i in ITOT) + sum(d[i,j,n] for i in S) - sum(c[i,j,n] for i in S) + sum(pl[j,l] for l in L if l[1] == n; init=0)  >= D[Symbol(n)][j] - ens[j,n] + sum(pl[j,l] for l in L if l[2] == n; init=0)) # Market Clearing constraint (if we assume curtailment of RES: replace == with >=) NOT OKAY SHOULD USE + P_RECEIVING - P_SENDING

     #con_MC = m.ext[:constraints][:con_MC] = @constraint(m, [j=J, n=N], sum(g[i,j,n] for i in ITOT) + sum(pl[j,l] for l in L if l[1] == n; init=0)  >= D[Symbol(n)][j] - ens[j,n] + sum(pl[j,l] for l in L if l[2] == n; init=0)) # Market Clearing constraint (if we assume curtailment of RES: replace == with >=) NOT OKAY SHOULD USE + P_RECEIVING - P_SENDING
     #con_MC = m.ext[:constraints][:con_MC] = @constraint(m, [j=J, n=N], sum(g[i,j,n] for i in I) + sum(pl[j,l] for l in L if l[1] == n; init=0) + sum(g_OCGT[j,n])  >= D[Symbol(n)][j] - ens[j,n] + sum(pl[j,l] for l in L if l[2] == n; init=0) + sum(g_PtH[j,n])) # Market Clearing constraint (if we assume curtailment of RES: replace == with >=) NOT OKAY SHOULD USE + P_RECEIVING - P_SENDING
     con_LoL = m.ext[:constraints][:con_LoL] = @constraint(m, [j=J,n=N], 0 <= ens[j,n] <= D[Symbol(n)][j]) # Loss of Load constraint
     #con_PFap_ac = m.ext[:constraints][:con_PFap_ac] = @constraint(m, [j=J,la=L_ac], pl[j,la] == Bl_ac * network["AC_Lines"]["AC_$(find_line_number(network["AC_Lines"], la))"]["Length"] * (θ[j,la[1]] - θ[j,la[2]]) * 400 * 400) # 'DC' Power Flow Approximation [MW] # θ should have a node as argument but l[x] is a node # TAKE LINE VOLTAGE VALUE AS 400kV
     #con_PFap_dc = m.ext[:constraints][:con_PFap_dc] = @constraint(m, [j=J,ld=L_dc], pl[j,ld] == Bl_dc * network["DC_Lines"]["DC_$(find_line_number(network["DC_Lines"], ld))"]["Length"] * (θ[j,ld[1]] - θ[j,ld[2]]) * 400 * 400) # 'DC' Power Flow Approximation [MW] # θ should have a node as argument but l[x] is a node 
     con_varlac1 = m.ext[:constraints][:con_varlac1] = @constraint(m, [j=J,la=L_ac], pl[j,la] <= Fl_MAX_AC[find_line_number(network["AC_Lines"], la)]) # Upper Limit Power Flow AC
     con_varlac2 = m.ext[:constraints][:con_varlac2] = @constraint(m, [j=J,la=L_ac], -Fl_MAX_AC[find_line_number(network["AC_Lines"], la)] <= pl[j,la]) # Lower Limit Power Flow AC
     con_varldc1 = m.ext[:constraints][:con_varldc1] = @constraint(m, [j=J,ld=L_dc], pl[j,ld] <= Fl_MAX_DC[find_line_number(network["DC_Lines"], ld)]) # Upper Limit Power Flow DC
     con_varldc2 = m.ext[:constraints][:con_varldc2] = @constraint(m, [j=J,ld=L_dc], -Fl_MAX_DC[find_line_number(network["DC_Lines"], ld)] <= pl[j,ld]) # Lower Limit Power Flow DC
     #con_θb = m.ext[:constraints][:con_θb] = @constraint(m, [j=J,n=N], -pi <= θ[j,n] <= pi) # Bound Voltage angles
     #con_θref = m.ext[:constraints][:con_θref] = @constraint(m, [j=J], θ[j,"BE"] == 0.0) # Voltage angle at ref node
     con_DGl = m.ext[:constraints][:con_DGl] = @constraint(m, [i=ID, j=J, n=N], g[i,j,n] <= AF[i]*cap_exist[n][i]) # Dispatchable generation limit
     #con_DGl = m.ext[:constraints][:con_DGl] = @constraint(m, [i=ID, j=J, n=N], g[i,j,n] <= AF[i]*cap[i,n]) # Dispatchable generation limit
     con_SGl = m.ext[:constraints][:con_SGl] = @constraint(m, [i=["Solar"], j=J, n=N], g[i,j,n] <= AFS[Symbol(n)][j]*AF[i]*cap[i,n]) # Solar generation limit
     con_WONGl = m.ext[:constraints][:con_WONGl] = @constraint(m, [i=["WindOnshore"], j=J, n=N], g[i,j,n] <= AFWON[Symbol(n)][j]*AF[i]*cap[i,n]) # Windon generation limit
     con_WOFFGl = m.ext[:constraints][:con_WOFFGl] = @constraint(m, [i=["WindOffshore"], j=J, n=intersect(countries, countries_windoff)], g[i,j,n] <= 1.2*AFWOFF[Symbol(n)][j]*AF[i]*cap[i,n]) # Windoff generation limit (assume hubheight at 150m --> factor 1.2)
     con_CWO_WOFF = m.ext[:constraints][:con_CWO_WOFF] = @constraint(m, [i=["WindOffshore"], j=J, n=setdiff(countries, countries_windoff)], g[i,j,n] == 0.0) # No wind offshore in these countries
 
     #con_alphaCO2 = m.ext[:constraints][:con_alphaCO2] = @constraint(m, [i=I_CO2, j=J, n=N], g[i,j,n] == 0.0)

     #con_ENS = m.ext[:constraints][:con_ENS] = @constraint(m, (1 - (sum(ens[j,n] for j in J, n in N)/sum(D[Symbol(n)][j] for n in N, j in J))) >= 0.90) # Reliability percentage of Power system
 
     #con_PtH = m.ext[:constraints][:con_PtH] = @constraint(m, [j=J, n=N], g_PtH[j,n] <= AFSt["PtH"]*cap_PtH[n]) # PtH generation limit
     #con_OCGT = m.ext[:constraints][:con_OCGT] = @constraint(m, [j=J, n=N], g_OCGT[j,n] <= AFSt["OCGT_H"]*cap_OCGT[n]) # OCGT generation limit

     #con_SB = m.ext[:constraints][:con_SB] = @constraint(m, [j=1:365, n=N], sum(g_PtH[8*(j-1)+i,n] for i in 1:8) == sum(g_OCGT[8*(j-1)+i,n] for i in 1:8)) # Daily Storage Balance in each country
     #con_SB = m.ext[:constraints][:con_SB] = @constraint(m, [j=1:365], sum(g_PtH[8*(j-1)+i,n] for i in 1:8, n in N) == sum(g_OCGT[8*(j-1)+i,n] for i in 1:8, n in N)) # Daily Storage Balance over whole EU 

     #con_SB = m.ext[:constraints][:con_SB] = @constraint(m, [n=N], sum(g_PtH[j,n] for j in 1:timestep(aggregate_3h)) == sum(g_OCGT[j,n] for j in 1:timestep(aggregate_3h))) # Yearly Storage Balance in each country
     
     #con_SB = m.ext[:constraints][:con_SB] = @constraint(m, sum(g_PtH[j,n] for j in 1:timestep(aggregate_3h), n in N) == sum(g_OCGT[j,n] for j in 1:timestep(aggregate_3h), n in N)) # Yearly Storage Balance over whole EU 

     # Add storage model with cyclic boundary conditions
     con_SB1 = m.ext[:constraints][:con_SB1] = @constraint(m, [i=S,j=J[2:end-1],n=N], e[i,j+1,n] - e[i,j,n] == eff_ch[i]*c[i,j,n] - (d[i,j,n]/eff_dis[i]))
     con_SB_init = m.ext[:constraints][:con_SB_init] = @constraint(m, [i=S,j=J[1],n=N], e[i,j,n] == (Emax_Storage[n][i]))
     con_SB_end = m.ext[:constraints][:con_SB_end] = @constraint(m, [i=["Batteries", "PHS_OL", "PHS_CL"],j=J[end],n=N], e[i,j,n] == (Emax_Storage[n][i]))
     #con_SB_end2 = m.ext[:constraints][:con_SB_end2] = @constraint(m, [i=["Hydro_Res", "Hydro_Pon"],j=J[end],n=N], e[i,j,n] == 0)
     #con_SB_begin = m.ext[:constraints][:con_SB_begin] = @constraint(m, [i=S,j=J[1],n=N], e[i,j+1,n] - (Emax_Storage[n][i]) == eff_ch[i]*c[i,j,n] - (d[i,j,n]/eff_dis[i]))
     #con_SB_end = m.ext[:constraints][:con_SB_end] = @constraint(m, [i=["Batteries", "PHS_OL", "PHS_CL"],j=J[end],n=N], (Emax_Storage[n][i]) - e[i,j,n] == eff_ch[i]*c[i,j,n] - (d[i,j,n]/eff_dis[i]))
     #con_SB_end2 = m.ext[:constraints][:con_SB_end2] = @constraint(m, [i=["Hydro_Res", "Hydro_Pon"],j=J[end],n=N], 0 - e[i,j,n] == eff_ch[i]*c[i,j,n] - (d[i,j,n]/eff_dis[i])) # Final Energy in the Hydro Res and Hydro Pon has to be equal to zero
 
 end

 # Build model
 println("Start BUILD MODEL")
 @time build_model!(m)
 println("BUILD MODEL Done")
 
 # Solve
 optimize!(m)

 # Stop the timer and calculate elapsed time
 elapsed_time = time() - t_start
 println("Elapsed time: $elapsed_time seconds")


 # check termination status
 print(
     """
 
     Termination status: $(termination_status(m))
 
     """
 )
 # check objective value
 @show value(m.ext[:objective])
 
 function calculate_country_capacity(data::Matrix{Float64}, countries::Vector{String}, weight::Float64)
    country_capacity = Dict{String, Float64}()
    for (j, country) in enumerate(countries)
        total_capacity = sum(data[:, j])
        country_capacity[country] = total_capacity * weight
    end
    return country_capacity
end
#countries_sum = ["ES", "FR", "DE", "BE", "UK", "NL", "DK", "NO", "SE", "FI", "IE", "CH", "AT", "IT", "PT"]
#calculate_country_capacity(Matrix(value.(m.ext[:variables][:cap])), countries_sum, 1.0)

 ## Step 5: Visualization
 
 N = m.ext[:sets][:Nodes]
 I = m.ext[:sets][:I]
 L_ac = m.ext[:sets][:AC_Lines] 
 L_dc = m.ext[:sets][:DC_Lines]
 S = m.ext[:sets][:S]
 IV = ["Solar", "WindOnshore", "WindOffshore"]
 IP = ["Solar", "WindOnshore", "WindOffshore", "Biomass"]
 ID = ["Nuclear", "SPP_lignite", "SPP_coal", "OCGT", "ICE", "Biofuel", "Hydro_RoR"]  
 ITOT = ["Solar", "WindOnshore", "WindOffshore", "Nuclear", "SPP_lignite", "SPP_coal", "OCGT", "ICE", "Biofuel", "Hydro_RoR"]

 plot_demand_2030 = plot_demand(dem) # plot demand national estimates 2030
 
 #gen_dict_abs_year = get_generators_capacity_storage_years(m, I, N, S, year)
 #gen_dict_abs = get_generators_capacity_storage(m, I, N, S)
 gen_dict_abs = get_generators_capacity_storage_fixed_network(m, ITOT, IV, ID, S, N)
 #gen_dict_rel = get_generators_capacity_storage_relative(m, I, N, S)
 data_matrix_abs = matrix_generators_data(gen_dict_abs)
 #data_matrix_rel = matrix_generators_data(gen_dict_rel)
 nodes = collect(keys(gen_dict_abs))
 generation_abs = plot_generator_capacities_new(gen_dict_abs, Generator_colors_dict, "Installed Generation & Storage Capacity (GW)")
 #generation_abs = plot_generator_capacities(data_matrix_abs, nodes, Generator_colors_storage, Generator_labels_storage, "Installed Generation & Storage Capacity (GW)")
 #generation_rel = plot_generator_capacities(data_matrix_rel, nodes, Generator_colors_storage, Generator_labels_storage, "Installed Generation & Storage Capacity (%)")
 
 power_flow_histogram = plot_power_flow_histogram(["UK", "FR"])  # ["ES", "FR"] ["UK", "FR"] ["NO", "NL"] ["PT", "ES"] ["BE", "DE"] ["SE", "NO"]
 trans_ac_dict = get_ac_transmission_capacity_fixed_network(m, L_ac)
 trans_dc_dict = get_dc_transmission_capacity_fixed_network(m, L_dc)
 trans_dict = merge(get_ac_transmission_capacity_fixed_network(m, L_ac), get_dc_transmission_capacity_fixed_network(m, L_dc))
 transmission_map_full = plot_transmission_network(trans_ac_dict, trans_dc_dict, countries_coord)
 transmission_map = plot_transmission_needed(trans_ac_dict, trans_dc_dict, countries_coord)
 
 plot(generation_abs, transmission_map, layout=(1,2), size=(1200,500))






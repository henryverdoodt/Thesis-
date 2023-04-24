## Master Thesis: Analyse the impact of climate change on a G&TEP of EU Power System
# Author: Henry Verdoodt
# Last update: April 20, 2023

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

using CSV, YAML, JuMP, DataFrames, Distributions, Gurobi, Images, Plots, PolyChaos, InteractiveUtils
using StatsPlots, Images

# Transmission Network West Europe
function plot_transmission_network()
	img 	= load("/Users/henryverdoodt/Documents/CODE/IMAGES/europe_map.jpeg");
    countries = Dict("Spain" => (380, 1400), "France" => (640, 1140), "Belgium" => (722, 960), "Germany" => (880, 960), 
                 "Netherlands" => (745, 885), "Denmark" => (866, 720), "Norway" => (890, 482), "United Kingdom" => (553, 835))

	plot(img, axis=([], false))
	plot!([countries["Spain"][1],countries["France"][1]],
		  [countries["Spain"][2],countries["France"][2]],
		  color="blue", linewidth=2, label="HVAC")
	plot!([countries["France"][1],countries["Belgium"][1]],
		  [countries["France"][2],countries["Belgium"][2]],
		  color="blue", linewidth=2, label=:none)
	plot!([countries["France"][1],countries["Germany"][1]],
		  [countries["France"][2],countries["Germany"][2]],
		  color="blue", linewidth=2, label=:none)
	plot!([countries["Germany"][1],countries["Belgium"][1]],
		  [countries["Germany"][2],countries["Belgium"][2]],
		  color="blue", linewidth=2, label=:none)
    plot!([countries["Belgium"][1],countries["Netherlands"][1]],
		  [countries["Belgium"][2],countries["Netherlands"][2]],
		  color="blue", linewidth=2, label=:none)
	plot!([countries["Netherlands"][1],countries["Germany"][1]],
		  [countries["Netherlands"][2],countries["Germany"][2]],
		  color="blue", linewidth=2, label=:none)
	plot!([countries["Germany"][1],countries["Denmark"][1]],
		  [countries["Germany"][2],countries["Denmark"][2]],
		  color="blue", linewidth=2, label=:none)

	plot!([countries["France"][1],countries["United Kingdom"][1]],
		  [countries["France"][2],countries["United Kingdom"][2]],
		  color="red", linewidth=2, label="HVDC")
	plot!([countries["United Kingdom"][1],countries["Belgium"][1]],
		  [countries["United Kingdom"][2],countries["Belgium"][2]],
		  color="red", linewidth=2, label=:none)
    plot!([countries["United Kingdom"][1],countries["Norway"][1]],
		  [countries["United Kingdom"][2],countries["Norway"][2]],
		  color="red", linewidth=2, label=:none)
	plot!([countries["Norway"][1],countries["Belgium"][1]],
		  [countries["Norway"][2],countries["Belgium"][2]],
		  color="red", linewidth=2, label=:none)
    plot!([countries["Norway"][1],countries["Netherlands"][1]],
		  [countries["Norway"][2],countries["Netherlands"][2]],
		  color="red", linewidth=2, label=:none)
	plot!([countries["Norway"][1],countries["Denmark"][1]],
		  [countries["Norway"][2],countries["Denmark"][2]],
		  color="red", linewidth=2, label=:none)
    plot!([countries["Netherlands"][1],countries["Denmark"][1]],
		  [countries["Netherlands"][2],countries["Denmark"][2]],
		  color="red", linewidth=2, label=:none)

	# end
	scatter!([nc[1] for nc in values(countries)],
			 [nc[2] for nc in values(countries)], markersize=3,
			 color="black", label=:none)
end

plot_transmission_network()

# Read the CSV file into a DataFrame
demand = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/DEMAND/PECD-country-demand_national_estimates-2025.csv", DataFrame)
solar = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/SOLAR/PECD-2021.3-country-LFSolarPV-2025.csv", DataFrame)
windon = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDON/PECD-2021.3-country-Onshore-2025.csv", DataFrame)
windoff = CSV.read("/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDOFF/PECD-2021.3-country-Offshore-2025.csv", DataFrame)
data = YAML.load_file(joinpath(@__DIR__, "overview_data.yaml"))
network = YAML.load_file(joinpath(@__DIR__, "network_west_europe.yaml"))
println("Done")

# Get the unique countries in the 'country' column
countries_demand = unique(demand.country) # ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "ME", "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UA", "UK"]
countries_solar= unique(solar.country)    # ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "ME", "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UA", "UK"]
countries_windon = unique(windon.country) # MT is in solar & demand but not in windon => # ["AL", "AT", "BA", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "ME", "MK", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "TR", "UA", "UK"]
countries_windoff= unique(windoff.country) # ["BE", "DE", "DK", "FI", "FR", "IE", "IT", "NL", "PT", "SE", "UK"]

############################# REFORMATING DATASETS (one column per node where each row is the hourly value (demand or cf)) #############################
# PARAMETERS
y = 2016.0
countries = ["ES", "FR", "BE", "DE", "NL", "UK", "DK", "NO"]

dem = DataFrame(col1 = Float64[])
for c in intersect(countries, countries_demand)
    demand1 = dropmissing(demand[demand.country .== c, [:year, :dem_MW]], disallowmissing=true)
    demand2 = dropmissing(demand1[demand1.year .== y, [:dem_MW]], disallowmissing=true)
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
    solar2 = dropmissing(solar1[solar1.year .== y, [:cf]], disallowmissing=true)
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
    windon2 = dropmissing(windon1[windon1.year .== y, [:cf]], disallowmissing=true)
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
    windoff2 = dropmissing(windoff1[windoff1.year .== y, [:cf]], disallowmissing=true)
    rename!(windoff2, :cf => c)
    if isempty(woff)
        woff = hcat(windoff2)
    else
        woff = hcat(woff, windoff2)
    end
end


## Step 2: create model & pass data to model
m = Model(optimizer_with_attributes(Gurobi.Optimizer))

# Step 2a: create sets
function define_sets!(m::Model, data::Dict, network::Dict)
    m.ext[:sets] = Dict()
    J = m.ext[:sets][:J] = 1:8760 # Timesteps
    I = m.ext[:sets][:I] = [id for id in keys(data["PowerSector"]) if id ∉ ["SPP_lignite", "SPP_coal", "CCGT_old"]] # Generators per type
    L_ac = m.ext[:sets][:AC_Lines] = [network["AC_Lines"]["AC_$(i)"]["Connection"] for i in 1:7] # AC Lines
    L_dc = m.ext[:sets][:DC_Lines] = [network["DC_Lines"]["DC_$(i)"]["Connection"] for i in 1:7] # DC Lines
    L = m.ext[:sets][:Lines] = union(m.ext[:sets][:AC_Lines],m.ext[:sets][:DC_Lines]) # Lines
    N = m.ext[:sets][:Nodes] =  [network["Nodes"]["Node_$(i)"] for i in 1:8]# Nodes
    return m # return model
end

# Step 2b: add time series
function process_time_series_data!(m::Model, dem::DataFrame, sol::DataFrame, won::DataFrame, woff::DataFrame)
    # extract the relevant sets
    I = m.ext[:sets][:I]
    J = m.ext[:sets][:J]
    N = m.ext[:sets][:Nodes]

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

    # Create parameter dictonary
    m.ext[:parameters] = Dict()

    # general parameters
    disc_r = m.ext[:parameters][:disc_r] = data["General"]["discount_rate"] # -, discount rate
    VOLL = m.ext[:parameters][:VOLL] = data["General"]["valueOfLostLoad"] # EUR/MWh, VOLL
    alphaCO2 = m.ext[:parameters][:alphaCO2] = data["ETS"]["P_calibration"] # €/ton, alphaCO2
    Bl_ac = m.ext[:parameters][:Bl_ac] = network["Bl_ac"]*10^(-6) # S/km, Susceptance of AC line 
    Bl_dc =  m.ext[:parameters][:Bl_dc] = network["Bl_dc"]*10^(-6) # S/km, Susceptance of DC line  
    Bl = m.ext[:parameters][:Bl] = cat(Bl_ac, Bl_dc, dims=1) # S/km, Susceptance of All line [USE cat() OR union() function??]
  
    # parameters of generators per unit
    d = data["PowerSector"]
    betaD = m.ext[:parameters][:betaD] = Dict(i => (d[SubString(i,1:length(i))]["fuelCosts"]/d[SubString(i,1:length(i))]["efficiency"]) for i in I) # EUR/MWh, Cost of generation of unit i
    deltaD = m.ext[:parameters][:deltaD] = Dict(i => (d[SubString(i,1:length(i))]["emissions"]/d[SubString(i,1:length(i))]["efficiency"]) for i in I) # ton/MWh, Emissions of generation of unit i
    GMAX = m.ext[:parameters][:GMAX] = Dict(i => d[SubString(i,1:length(i))]["gMAX"] for i in I) # MW, Maximum Power Output of one generation unit
    VC = m.ext[:parameters][:VC] = Dict(i => ((d[SubString(i,1:length(i))]["fuelCosts"] + d[SubString(i,1:length(i))]["emissions"]*alphaCO2)/d[SubString(i,1:length(i))]["efficiency"]) for i in I) # EUR/MWh, Variable Cost of unit i
    OC = m.ext[:parameters][:OC] = Dict(i => 10^3*d[SubString(i,1:length(i))]["OC"] for i in I) # EUR/MW, Overnight investment cost of a unit
    LifeT = m.ext[:parameters][:LifeT] = Dict(i => d[SubString(i,1:length(i))]["Lifetime"] for i in I) # years, Lifetime of Asset
    IC = m.ext[:parameters][:IC] = Dict(i => 10^3*(d[SubString(i,1:length(i))]["OC"]*disc_r)/(1-(1/(1+disc_r)^(d[SubString(i,1:length(i))]["Lifetime"]))) for i in I) # EUR/MWy, Investment Cost in generation unit i
    AF = m.ext[:parameters][:AF] = Dict(i => d[SubString(i,1:length(i))]["AF"] for i in I) # Availability facor of a generator

    
    # parameters of transmission lines AC
    Fl_MAX_AC = m.ext[:parameters][:Fl_MAX_AC] = [network["AC_Lines"]["AC_$(i)"]["Capacity"] for i in 1:7] # MW, AC Lines
    L_length_AC = m.ext[:parameters][:L_length_AC] = [network["AC_Lines"]["AC_$(i)"]["Length"] for i in 1:7] # km, Length of AC transmission line
    L_price_AC = m.ext[:parameters][:L_price_AC] = [network["AC_Lines"]["AC_$(i)"]["Price"] for i in 1:7] # EUR/km.MW, Price of AC line per unit length
    OC_var_AC = m.ext[:parameters][:OC_var_AC] = [network["AC_Lines"]["AC_$(i)"]["Price"] * network["AC_Lines"]["AC_$(i)"]["Length"] for i in 1:7] # EUR/MW, Investment Cost for AC Transmission line
    IC_var_AC = m.ext[:parameters][:IC_var_AC] = [(network["AC_Lines"]["AC_$(i)"]["Price"] * network["AC_Lines"]["AC_$(i)"]["Length"]*disc_r)/(1-(1/(1+disc_r)^(network["Lifetime_ac"]))) for i in 1:7] # EUR/MWy, Investment Cost of AC line

    # parameters of transmission lines DC 
    Fl_MAX_DC = m.ext[:parameters][:Fl_MAX_DC] = [network["DC_Lines"]["DC_$(i)"]["Capacity"] for i in 1:7] # MW, DC Lines
    L_length_DC = m.ext[:parameters][:L_length_DC] = [network["DC_Lines"]["DC_$(i)"]["Length"] for i in 1:7] # km, Length of DC transmission line
    L_price_DC = m.ext[:parameters][:L_price_DC] = [network["DC_Lines"]["DC_$(i)"]["Price"] for i in 1:7] # EUR/km.MW, Price of DC line per unit length
    OC_var_DC = m.ext[:parameters][:OC_var_DC] = [network["DC_Lines"]["DC_$(i)"]["Price"] * network["DC_Lines"]["DC_$(i)"]["Length"] for i in 1:7] # EUR/MW, Investment Cost for DC Transmission line
    IC_var_DC = m.ext[:parameters][:IC_var_DC] = [(network["DC_Lines"]["DC_$(i)"]["Price"] * network["DC_Lines"]["DC_$(i)"]["Length"]*disc_r)/(1-(1/(1+disc_r)^(network["Lifetime_dc"]))) for i in 1:7] # EUR/MWy, Investment Cost of DC line
    Fl_MAX = m.ext[:parameters][:Fl_MAX] = cat(Fl_MAX_AC, Fl_MAX_DC, dims=1) # MW, All Lines, Maximum capacity [USE cat() OR union() function??]
    # return model
    return m
end


# call functions
define_sets!(m, data, network)  
process_time_series_data!(m, dem, sol, won, woff)
process_parameters!(m, data, network)


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

#test = find_line_number(network["DC_Lines"], ["UK", "FR"])

## Step 3: construct your model
function build_model!(m::Model)
    # Create m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    J = m.ext[:sets][:J] # Timesteps
    I = m.ext[:sets][:I] # Generators per type
    ID = ["Nuclear", "CCGT_new", "OCGT", "ICE", "Biomass"]  # ID = ["Nuclear", "SPP_lignite", "SPP_coal", "CCGT_old", "CCGT_new", "OCGT", "ICE", "Biomass"]; As we assume Greenfield method we only take the new tech. and dont take coal into account
    IV = ["Solar", "WindOnshore", "WindOffshore"]
    L_ac = m.ext[:sets][:AC_Lines] # AC Lines
    L_dc = m.ext[:sets][:DC_Lines] # DC Lines
    L = m.ext[:sets][:Lines] # All Lines
    N = m.ext[:sets][:Nodes] # Nodes

    # Extract time series data
    D = m.ext[:timeseries][:D]
    AFWON = m.ext[:timeseries][:AFWON]
    AFWOFF = m.ext[:timeseries][:AFWOFF]
    AFS = m.ext[:timeseries][:AFS]

    # Extract parameters
        ### general parameters
    VOLL = m.ext[:parameters][:VOLL] 
    Bl_ac = m.ext[:parameters][:Bl_ac] 
    Bl_dc =  m.ext[:parameters][:Bl_dc] 
    Bl = m.ext[:parameters][:Bl]
    
        ### parameters of generators per unit
    VC = m.ext[:parameters][:VC]
    IC = m.ext[:parameters][:IC]
    AF = m.ext[:parameters][:AF]

        ### parameters of transmission lines
    IC_var_AC = m.ext[:parameters][:IC_var_AC]
    IC_var_DC = m.ext[:parameters][:IC_var_DC] 
   

    # create variables 
    g = m.ext[:variables][:g] = @variable(m, [i=I,j=J,n=N], lower_bound=0, base_name="generation") # Power produced by candidate generating unit i at time j in node n [MW]
    cap = m.ext[:variables][:cap] = @variable(m, [i=I,n=N], lower_bound=0, base_name="capacity of generating unit") # Capacity of candidate generating unit i in node n[MW]
    ens = m.ext[:variables][:ens] = @variable(m, [j=J,n=N], lower_bound=0, base_name="energy not served") # Energy Not Served at time j in node n [MWh] (OR Should I use Load shed of demand instead in MW?)
    pl = m.ext[:variables][:pl] = @variable(m, [j=J,l=L], lower_bound=0, base_name="power flow in transmission") # Power flow through Transmission Line l at time j [MW]
    θ = m.ext[:variables][:θ] = @variable(m, [j=J,n=N], lower_bound=0, base_name="voltage angle") # Voltage angle at node n and time j [rad]
    #xl = m.ext[:variables][:xl] = @variable(m, [l=L], Bin, lower_bound=0, base_name="decision built transmission") # Decision variable, built transmission line equal to 1 if not built equal to 0
    varlac = m.ext[:variables][:varlac] = @variable(m, [la=L_ac], lower_bound=0, base_name="capacity of ac line") # Capacity of AC line [MW]
    varldc = m.ext[:variables][:varldc] = @variable(m, [ld=L_dc], lower_bound=0, base_name="capacity of dc line") # Capacity of DC line [MW]


    # Objective
    #obj = m.ext[:objective] = @objective(m, Min, sum(IC[i]*cap[i,n] + (IC_l_AC[la]+IC_l_DC[ld])*xl[l] + VC[i]*g[i,j,n] + VOLL*ens[j,n] for i in I, j in J, n in N, la in 1:6, ld in 1:5, l in L))
    obj = m.ext[:objective] = @objective(m, Min, sum(IC[i]*cap[i,n] for i in I, n in N)  + sum(IC_var_AC[find_line_number(network["AC_Lines"], la)]*varlac[la] for la in L_ac) + sum(IC_var_DC[find_line_number(network["DC_Lines"], ld)]*varldc[ld] for ld in L_dc) + sum(VC[i]*g[i,j,n] for i in I, j in J, n in N) + sum(VOLL*ens[j,n] for j in J, n in N))

    # Matrice format to reduce calculation time
    #=
    # Define variables as matrices
    cap_mat = value.(reshape([cap[i,n] for i in I for n in N], (length(I), length(N))))
    g_mat = value.(reshape([g[i,j,n] for i in I for j in J for n in N], (length(I), length(J), length(N))))
    ens_mat = value.(reshape([ens[j,n] for j in J for n in N], (length(J), length(N))))
    varlac_vec = value.(varlac)
    varldc_vec = value.(varldc)

    # Calculate the objective using matrix operations
    obj = m.ext[:objective] = @objective(m, Min, dot(IC, cap_mat) + dot(IC_var_AC, varlac_vec) + dot(IC_var_DC, varldc_vec) + dot(VC, g_mat) + VOLL * sum(ens_mat))
     =#
    #=
    println(sum(IC[i] for i in I))
    println(sum(IC[i]*cap[i,n] for i in I, n in N))
    println(sum(IC_var_AC[find_line_number(network["AC_Lines"], la)]*varlac[la] for la in L_ac))
    find_line_number(network["AC_Lines"], ["Milaan", "Florence"])
    println(sum(IC_var_DC[ld] for ld in 1:length(L_dc)))
    println(la for la in L_ac)
    println(n for n in N)
    println(la for la in L_ac) #sum(IC_l_AC[la] for 
    println(IC_l_AC[2])
    println([network["AC_Lines"]["AC_$(la)"]["Connection"] for la in 1:6])
    =#
                                                                                            # sum(Pl[l] for l in L if l[1] == n; init=0)
    # Constraints
    con_DGl = m.ext[:constraints][:con_DGl] = @constraint(m, [i=ID, j=J, n=N], g[i,j,n] <= AF[i]*cap[i,n]) # Dispatchable generation limit
    con_SGl = m.ext[:constraints][:con_SGl] = @constraint(m, [i=["Solar"], j=J, n=N], g[i,j,n] <= AFS[Symbol(n)][j]*AF[i]*cap[i,n]) # Solar generation limit
    con_WONGl = m.ext[:constraints][:con_WONGl] = @constraint(m, [i=["WindOnshore"], j=J, n=N], g[i,j,n] <= AFWON[Symbol(n)][j]*AF[i]*cap[i,n]) # Wind generation limit
    con_WOFFGl = m.ext[:constraints][:con_WOFFGl] = @constraint(m, [i=["WindOffshore"], j=J, n=intersect(countries, countries_windoff)], g[i,j,n] <= AFWOFF[Symbol(n)][j]*AF[i]*cap[i,n]) # Wind generation limit
    con_CWO_WOFF = m.ext[:constraints][:con_CWO_WOFF] = @constraint(m, [i=["WindOffshore"], j=J, n=setdiff(countries, countries_windoff)], g[i,j,n] == 0.0) # No wind offshore in these countries
    #con_SOL_ES = m.ext[:constraints][:con_SOL_ES] = @constraint(m, [i=["Nuclear", "CCGT_new", "OCGT", "ICE", "Biomass", "WindOnshore"], j=J, n=["ES"]], g[i,j,n] == 0.0)

    con_MC = m.ext[:constraints][:con_MC] = @constraint(m, [j=J, n=N], sum(g[i,j,n] for i in I) + sum(pl[j,l] for l in L if l[1] == n; init=0) >= D[Symbol(n)][j] - ens[j,n] + sum(pl[j,l] for l in L if l[2] == n; init=0)) # Market Clearing constraint (if we assume curtailment of RES: replace == with >=) NOT OKAY SHOULD USE + P_RECEIVING - P_SENDING
    con_LoL = m.ext[:constraints][:con_LoL] = @constraint(m, [j=J,n=N], 0 <= ens[j,n] <= D[Symbol(n)][j]) # Loss of Load constraint
    #con_PFap_ac = m.ext[:constraints][:con_PFap_ac] = @constraint(m, [j=J,la=L_ac], pl[j,la] == Bl_ac*(θ[j,la[1]] - θ[j,la[2]])) # 'DC' Power Flow Approximation # θ should have a node as argument but l[x] is a node 
    #con_PFap_dc = m.ext[:constraints][:con_PFap_dc] = @constraint(m, [j=J,ld=L_dc], pl[j,ld] == Bl_dc*(θ[j,ld[1]] - θ[j,ld[2]])) # 'DC' Power Flow Approximation # θ should have a node as argument but l[x] is a node 
    con_varlac1 = m.ext[:constraints][:con_varlac1] = @constraint(m, [j=J,la=L_ac], pl[j,la] <= varlac[la]) # Upper Limit Power Flow AC
    con_varlac2 = m.ext[:constraints][:con_varlac2] = @constraint(m, [j=J,la=L_ac], -varlac[la] <= pl[j,la]) # Lower Limit Power Flow AC
    con_varldc1 = m.ext[:constraints][:con_varldc1] = @constraint(m, [j=J,ld=L_dc], pl[j,ld] <= varldc[ld]) # Upper Limit Power Flow DC
    con_varldc2 = m.ext[:constraints][:con_varldc2] = @constraint(m, [j=J,ld=L_dc], -varldc[ld] <= pl[j,ld]) # Lower Limit Power Flow DC
    #con_θb = m.ext[:constraints][:con_θb] = @constraint(m, [j=J,n=N], -pi <= θ[j,n] <= pi) # Bound Voltage angles
    #con_θref = m.ext[:constraints][:con_θref] = @constraint(m, [j=J], θ[j,"BE"] == 0.0) # Voltage angle at ref node
    #con_DGl = m.ext[:constraints][:con_DGl] = @constraint(m, [i=ID, j=J, n=N], g[i,j,n] <= AF[i]*cap[i,n]) # Dispatchable generation limit
    #con_SGl = m.ext[:constraints][:con_SGl] = @constraint(m, [i=["Solar"], j=J, n=N], g[i,j,n] <= AFS[Symbol(n)][j]*AF[i]*cap[i,n]) # Solar generation limit
    #con_WONGl = m.ext[:constraints][:con_WONGl] = @constraint(m, [i=["WindOnshore"], j=J, n=N], g[i,j,n] <= AFWON[Symbol(n)][j]*AF[i]*cap[i,n]) # Wind generation limit
    #con_WOFFGl = m.ext[:constraints][:con_WOFFGl] = @constraint(m, [i=["WindOffshore"], j=J, n=intersect(countries, countries_windoff)], g[i,j,n] <= AFWOFF[Symbol(n)][j]*AF[i]*cap[i,n]) # Wind generation limit
    #con_CWO_WOFF = m.ext[:constraints][:con_CWO_WOFF] = @constraint(m, [i=["WindOffshore"], j=J, n=setdiff(countries, countries_windoff)], g[i,j,n] == 0.0) # No wind offshore in these countries
    #con_SOL_ES = m.ext[:constraints][:con_SOL_ES] = @constraint(m, [i=["Nuclear", "CCGT_new", "OCGT", "ICE", "Biomass", "WindOnshore"], j=J, n=["ES"]], g[i,j,n] == 0.0)
end

# Build model
build_model!(m)

# Solve
optimize!(m)
# check termination status
print(
    """

    Termination status: $(termination_status(m))

    """
)
# check objective value
@show value(m.ext[:objective])

## Step 5: Visualization
 
N = m.ext[:sets][:Nodes]
I = m.ext[:sets][:I]
L_ac = m.ext[:sets][:AC_Lines] 
L_dc = m.ext[:sets][:DC_Lines]

# Dictionary to store the generators capacity for each node 
gen_dict = Dict{String, Dict{String, Float64}}()  
for node in N
    gen_dict[node] = Dict((gen => 0.0 for gen in I)...)
    for unit in I
        gen_dict[node][unit] = value.(m.ext[:variables][:cap][unit,node])
    end
end

# Create a dictionaris to store the AC and DC transmission capacity
trans_ac_dict = Dict(la => value.(m.ext[:variables][:varlac][la]) for la in L_ac)
trans_dc_dict = Dict(ld => value.(m.ext[:variables][:varldc][ld]) for ld in L_dc)

# Extract the node names and the generating unit types
nodes = collect(keys(gen_dict))
unit_types = collect(keys(gen_dict[nodes[1]]))

# Create the data matrix
data = zeros(length(nodes), length(unit_types))
for (i, node) in enumerate(nodes)
    for (j, unit_type) in enumerate(unit_types)
        data[i, j] = gen_dict[node][unit_type]
    end
end

# Color Palette
Generator_colors = [:green, :lightblue, :orange, :blue, :purple, :gold, :red, :gray]

# Plot generation stacked bar chart
generation = groupedbar(data,
bar_position = :stack,
bar_width=0.5,
bar_edges=true,
#xlabel="Node", 
ylabel="Capacity (MW)",
xticks=(1:length(nodes), nodes),
label= ["Biomass" "WindOffshore" "CCGT_new" "WindOnshore" "Nuclear" "Solar" "OCGT" "ICE"],
color_palette= Generator_colors)

# Transmission Network Italy
function plot_transmission_needed(dict1::Dict, dict2::Dict)
	img 	= load("/Users/henryverdoodt/Documents/CODE/IMAGES/europe_map.jpeg");
    countries = Dict("Spain" => (380, 1400), "France" => (640, 1140), "Belgium" => (722, 960), "Germany" => (880, 960), 
                 "Netherlands" => (745, 885), "Denmark" => (866, 720), "Norway" => (890, 482), "United Kingdom" => (553, 835))

	plot(img, axis=([], false))
    if dict1[["ES", "FR"]] != 0.0
        plot!([countries["Spain"][1],countries["France"][1]],
		    [countries["Spain"][2],countries["France"][2]],
		    color="blue", linewidth=2, label="HVAC")
    end
	if dict1[["FR", "BE"]] != 0.0
	    plot!([countries["France"][1],countries["Belgium"][1]],
		    [countries["France"][2],countries["Belgium"][2]],
		    color="blue", linewidth=2, label=:none)
    end
    if dict1[["FR", "DE"]] != 0.0
	    plot!([countries["France"][1],countries["Germany"][1]],
		    [countries["France"][2],countries["Germany"][2]],
		    color="blue", linewidth=2, label=:none)
    end
    if dict1[["BE", "DE"]] != 0.0
	    plot!([countries["Germany"][1],countries["Belgium"][1]],
		     [countries["Germany"][2],countries["Belgium"][2]],
		     color="blue", linewidth=2, label=:none)
    end
    if dict1[["BE", "NL"]] != 0.0
        plot!([countries["Belgium"][1],countries["Netherlands"][1]],
		    [countries["Belgium"][2],countries["Netherlands"][2]],
		    color="blue", linewidth=2, label=:none)
    end
    if dict1[["NL", "DE"]] != 0.0
	    plot!([countries["Netherlands"][1],countries["Germany"][1]],
		    [countries["Netherlands"][2],countries["Germany"][2]],
		    color="blue", linewidth=2, label=:none)
    end
    if dict1[["DK", "DE"]] != 0.0
	    plot!([countries["Germany"][1],countries["Denmark"][1]],
		    [countries["Germany"][2],countries["Denmark"][2]],
		    color="blue", linewidth=2, label=:none)
    end


    if dict2[["UK", "FR"]] != 0.0
	    plot!([countries["France"][1],countries["United Kingdom"][1]],
		    [countries["France"][2],countries["United Kingdom"][2]],
		    color="red", linewidth=2, label="HVDC")
    end
    if dict2[["UK", "BE"]] != 0.0
	    plot!([countries["United Kingdom"][1],countries["Belgium"][1]],
		    [countries["United Kingdom"][2],countries["Belgium"][2]],
		    color="red", linewidth=2, label=:none)
    end
    if dict2[["UK", "NO"]] != 0.0
        plot!([countries["United Kingdom"][1],countries["Norway"][1]],
		    [countries["United Kingdom"][2],countries["Norway"][2]],
		    color="red", linewidth=2, label=:none)
    end
    if dict2[["NO", "BE"]] != 0.0
	    plot!([countries["Norway"][1],countries["Belgium"][1]],
		    [countries["Norway"][2],countries["Belgium"][2]],
		    color="red", linewidth=2, label=:none)
    end
    if dict2[["NO", "NL"]] != 0.0
        plot!([countries["Norway"][1],countries["Netherlands"][1]],
		    [countries["Norway"][2],countries["Netherlands"][2]],
		    color="red", linewidth=2, label=:none)
    end
    if dict2[["NO", "DK"]] != 0.0
	    plot!([countries["Norway"][1],countries["Denmark"][1]],
		    [countries["Norway"][2],countries["Denmark"][2]],
		    color="red", linewidth=2, label=:none)
    end
    if dict2[["DK", "NL"]] != 0.0
        plot!([countries["Netherlands"][1],countries["Denmark"][1]],
		    [countries["Netherlands"][2],countries["Denmark"][2]],
		    color="red", linewidth=2, label=:none)
    end

	# Dots for Nodes
	scatter!([nc[1] for nc in values(countries)],
			 [nc[2] for nc in values(countries)], markersize=3,
			 color="black", label=:none)
end

function plot_transmission_capacities(dict1::Dict, dict2::Dict)
    countries = unique([k[1] for k in keys(dict1)] ∪ [k[2] for k in keys(dict1)] ∪ [k[1] for k in keys(dict2)] ∪ [k[2] for k in keys(dict2)])

    # Create a matrix of capacity values
    capacity_matrix = zeros(length(countries), length(countries))
    for (k, v) in dict1
        i = findfirst(countries .== k[1])
        j = findfirst(countries .== k[2])
        capacity_matrix[i, j] = v
    end
    for (k, v) in dict2
        i = findfirst(countries .== k[1])
        j = findfirst(countries .== k[2])
        capacity_matrix[i, j] += v
    end
    
    # Plot the heatmap
    heatmap(countries, countries, capacity_matrix, c=:blues, aspect_ratio=:equal, xlabel="To country", ylabel="From country", title="Capacity needed for different links (MW)")
end

# Plot transmission map and capacities
transmission = plot_transmission_network()
transmission_map = plot_transmission_needed(trans_ac_dict, trans_dc_dict)
transmission_cap = plot_transmission_capacities(trans_ac_dict, trans_dc_dict)

# Plot side by side
plot(generation, transmission_map, transmission_cap, layout=(1,3), size=(1800,600))



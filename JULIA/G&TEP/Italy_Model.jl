## Master Thesis: Analyse the impact of climate change on a G&TEP of EU Power System
# Author: Henry Verdoodt
# Last update: April 4, 2023

## Step 0: Activate environment - ensure consistency accross computers
#=
using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
Pkg.instantiate() # If a Manifest.toml file exist in the current project, download all the packages declared in that manifest. Else, resolve a set of feasible packages from the Project.toml files and install them.
Pkg.add("InteractiveUtils")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Distributions")
Pkg.add("Gurobi")
#Pkg.add("GLPK")
Pkg.add("Images")
Pkg.add("JuMP")
Pkg.add("Plots")
Pkg.add("PolyChaos")
Pkg.add("YAML") 
Pkg.add("StatsPlots")
=#

using CSV, YAML, JuMP, DataFrames, Distributions, Images, Plots, PolyChaos, InteractiveUtils #, GLPK
using StatsPlots, Gurobi


# Transmission Network Italy
function plot_transmission_network_italy()
	img 	= load("/Users/henryverdoodt/Documents/OneDrive - KU Leuven/CODE/JULIA/OPES/PLC_uncertainty/Italy.png");
	cities 	= Dict( "Milaan" 	=> (190,115),
					"Florence" 	=> (285,220),
					"Rome" 		=> (336,336),
					"Bari" 		=> (540,390),
					"Sardinia" 	=> (200,400),
					"Palermo" 	=> (380,565))
	plot(img, axis=([], false))
	plot!([cities["Milaan"][1],cities["Florence"][1]],
		  [cities["Milaan"][2],cities["Florence"][2]],
		  color="blue", linewidth=3, label="existing ac-lines")
	plot!([cities["Milaan"][1],cities["Rome"][1]],
		  [cities["Milaan"][2],cities["Rome"][2]],
		  color="blue", linewidth=3, label=:none)
	plot!([cities["Milaan"][1],cities["Bari"][1]],
		  [cities["Milaan"][2],cities["Bari"][2]],
		  color="blue", linewidth=3, label=:none)
	plot!([cities["Florence"][1],cities["Bari"][1]],
		  [cities["Florence"][2],cities["Bari"][2]],
		  color="blue", linewidth=3, label=:none)
	plot!([cities["Sardinia"][1],cities["Florence"][1]],
		  [cities["Sardinia"][2],cities["Florence"][2]],
		  color="yellow", linewidth=3, label="existing dc-lines")
	plot!([cities["Sardinia"][1],cities["Rome"][1]],
		  [cities["Sardinia"][2],cities["Rome"][2]],
		  color="yellow", linewidth=3, label=:none)
	plot!([cities["Florence"][1],cities["Rome"][1]],
		  [cities["Florence"][2],cities["Rome"][2]],
		  line=("blue", 3, :dash), label="proposed ac-lines")
	plot!([cities["Bari"][1],cities["Rome"][1]],
		  [cities["Bari"][2],cities["Rome"][2]],
		  line=("blue", 3, :dash), label=:none)
	plot!([cities["Palermo"][1],cities["Sardinia"][1]],
		  [cities["Palermo"][2],cities["Sardinia"][2]],
		  line=("yellow", 3, :dash), label="proposed dc-lines")
	plot!([cities["Palermo"][1],cities["Rome"][1]],
		  [cities["Palermo"][2],cities["Rome"][2]],
		  line=("yellow", 3, :dash), label=:none)
	# if x == 1
		plot!([cities["Bari"][1],cities["Palermo"][1]],
			  [cities["Bari"][2],cities["Palermo"][2]],
			  line=("yellow", 3, :dash), label=:none)
	# end
	scatter!([nc[1] for nc in values(cities)],
			 [nc[2] for nc in values(cities)],
			 color="orange", label=:none)
end

#plot_transmission_network_italy()

# Load input data
demand = CSV.read(joinpath(@__DIR__, "case_6_demand_10.csv"), DataFrame)
wind_cf = CSV.read(joinpath(@__DIR__, "case_6_wind_10.csv"), DataFrame)
pv_cf = CSV.read(joinpath(@__DIR__, "case_6_PV_10.csv"), DataFrame)
data = YAML.load_file(joinpath(@__DIR__, "overview_data.yaml"))
network = YAML.load_file(joinpath(@__DIR__, "network_italy.yaml"))
println("Done")

## Step 2: create model & pass data to model
#m = Model(optimizer_with_attributes(GLPK.Optimizer))
m = Model(Gurobi.Optimizer)

# Step 2a: create sets
function define_sets!(m::Model, data::Dict, network::Dict, demand::DataFrame, wind_cf::DataFrame, pv_cf::DataFrame)
    m.ext[:sets] = Dict()
    J = m.ext[:sets][:J] = 1:8760 # Timesteps
    I = m.ext[:sets][:I] = [id for id in keys(data["PowerSector"])] # Generators per type
    L_ac = m.ext[:sets][:AC_Lines] = [network["AC_Lines"]["AC_$(i)"]["Connection"] for i in 1:6] # AC Lines
    L_dc = m.ext[:sets][:DC_Lines] = [network["DC_Lines"]["DC_$(i)"]["Connection"] for i in 1:5] # DC Lines
    L = m.ext[:sets][:Lines] = union(m.ext[:sets][:AC_Lines],m.ext[:sets][:DC_Lines]) # Lines
    N = m.ext[:sets][:Nodes] =  [network["Nodes"]["Node_$(i)"] for i in 1:6]# Nodes
    return m # return model
end

# Step 2b: add time series
function process_time_series_data!(m::Model, demand::DataFrame, wind_cf::DataFrame, pv_cf::DataFrame)
    # extract the relevant sets
    I = m.ext[:sets][:I]
    J = m.ext[:sets][:J]
    N = m.ext[:sets][:Nodes]

    # create dictionary to store time series
    m.ext[:timeseries] = Dict()
    m.ext[:timeseries][:D] = Dict()
    m.ext[:timeseries][:AFW] = Dict()
    m.ext[:timeseries][:AFS] = Dict()

    # Demand TimeSeries; REF Demand := 2000 MWh
    m.ext[:timeseries][:D][:Florence] = demand.Florence[1:8760]*2000 
    m.ext[:timeseries][:D][:Milaan] = demand.Milaan[1:8760]*2000 
    m.ext[:timeseries][:D][:Rome] = demand.Rome[1:8760]*2000 
    m.ext[:timeseries][:D][:Bari] = demand.Bari[1:8760]*2000 
    m.ext[:timeseries][:D][:Sardinia] = demand.Sardinia[1:8760]*2000 
    m.ext[:timeseries][:D][:Palermo] = demand.Palermo[1:8760]*2000 

    # Wind_AF TimeSeries
    m.ext[:timeseries][:AFW][:Florence] = wind_cf.Florence[1:8760]
    m.ext[:timeseries][:AFW][:Milaan] = wind_cf.Milaan[1:8760]
    m.ext[:timeseries][:AFW][:Rome] = wind_cf.Rome[1:8760]
    m.ext[:timeseries][:AFW][:Bari] = wind_cf.Bari[1:8760]
    m.ext[:timeseries][:AFW][:Sardinia] = wind_cf.Sardinia[1:8760]
    m.ext[:timeseries][:AFW][:Palermo] = wind_cf.Palermo[1:8760]

    # PV_AF TimeSeries
    m.ext[:timeseries][:AFS][:Florence] = pv_cf.Florence[1:8760]
    m.ext[:timeseries][:AFS][:Milaan] = pv_cf.Milaan[1:8760]
    m.ext[:timeseries][:AFS][:Rome] = pv_cf.Rome[1:8760]
    m.ext[:timeseries][:AFS][:Bari] = pv_cf.Bari[1:8760]
    m.ext[:timeseries][:AFS][:Sardinia] = pv_cf.Sardinia[1:8760]
    m.ext[:timeseries][:AFS][:Palermo] = pv_cf.Palermo[1:8760]

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
    Fl_MAX_AC = m.ext[:parameters][:Fl_MAX_AC] = [network["AC_Lines"]["AC_$(i)"]["Capacity"] for i in 1:6] # MW, AC Lines
    L_length_AC = m.ext[:parameters][:L_length_AC] = [network["AC_Lines"]["AC_$(i)"]["Length"] for i in 1:6] # km, Length of AC transmission line
    L_price_AC = m.ext[:parameters][:L_price_AC] = [network["AC_Lines"]["AC_$(i)"]["Price"] for i in 1:6] # EUR/km.MW, Price of AC line per unit length
    OC_var_AC = m.ext[:parameters][:OC_var_AC] = [network["AC_Lines"]["AC_$(i)"]["Price"] * network["AC_Lines"]["AC_$(i)"]["Length"] for i in 1:6] # EUR/MW, Investment Cost for AC Transmission line
    IC_var_AC = m.ext[:parameters][:IC_var_AC] = [(network["AC_Lines"]["AC_$(i)"]["Price"] * network["AC_Lines"]["AC_$(i)"]["Length"]*disc_r)/(1-(1/(1+disc_r)^(network["Lifetime_ac"]))) for i in 1:6] # EUR/MWy, Investment Cost of AC line

    # parameters of transmission lines DC 
    Fl_MAX_DC = m.ext[:parameters][:Fl_MAX_DC] = [network["DC_Lines"]["DC_$(i)"]["Capacity"] for i in 1:5] # MW, DC Lines
    L_length_DC = m.ext[:parameters][:L_length_DC] = [network["DC_Lines"]["DC_$(i)"]["Length"] for i in 1:5] # km, Length of DC transmission line
    L_price_DC = m.ext[:parameters][:L_price_DC] = [network["DC_Lines"]["DC_$(i)"]["Price"] for i in 1:5] # EUR/km.MW, Price of DC line per unit length
    OC_var_DC = m.ext[:parameters][:OC_var_DC] = [network["DC_Lines"]["DC_$(i)"]["Price"] * network["DC_Lines"]["DC_$(i)"]["Length"] for i in 1:5] # EUR/MW, Investment Cost for DC Transmission line
    IC_var_DC = m.ext[:parameters][:IC_var_DC] = [(network["DC_Lines"]["DC_$(i)"]["Price"] * network["DC_Lines"]["DC_$(i)"]["Length"]*disc_r)/(1-(1/(1+disc_r)^(network["Lifetime_dc"]))) for i in 1:5] # EUR/MWy, Investment Cost of DC line
    Fl_MAX = m.ext[:parameters][:Fl_MAX] = cat(Fl_MAX_AC, Fl_MAX_DC, dims=1) # MW, All Lines, Maximum capacity [USE cat() OR union() function??]
    # return model
    return m
end


# call functions
define_sets!(m, data, network, demand, wind_cf, pv_cf)
process_time_series_data!(m, demand, wind_cf, pv_cf)
process_parameters!(m, data, network)


function find_line_number(network::Dict, cities::Vector)
    for (name, line) in network
        if line["Connection"] == cities
            # Return the line number if the cities match
            return parse(Int, split(name, "_")[2])
        end
    end
    # Return Nothing if no match is found
    return nothing
end

#test = find_line_number(network["DC_Lines"], ["Rome", "Sardinia"])

## Step 3: construct your model
function build_model!(m::Model)
    # Create m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Get the list of generator types to remove
    types_to_remove = ["SPP_lignite", "SPP_coal", "CCGT_old"]
    # Remove the specified generator types from the set I
    m.ext[:sets][:I] = filter(x -> !(x in types_to_remove), m.ext[:sets][:I])

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
    AFW = m.ext[:timeseries][:AFW]
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

    # Constraints
    con_MC = m.ext[:constraints][:con_MC] = @constraint(m, [j=J, n=N], sum(g[i,j,n] for i in I) + sum(pl[j,l] for l in L) == D[Symbol(n)][j] - ens[j,n]) # Market Clearing constraint (if we assume curtailment of RES: replace == with >=) NOT OKAY SHOULD USE + P_RECEIVING - P_SENDING
    con_LoL = m.ext[:constraints][:con_LoL] = @constraint(m, [j=J,n=N], 0 <= ens[j,n] <= D[Symbol(n)][j]) # Loss of Load constraint
    con_PFap_ac = m.ext[:constraints][:con_PFap_ac] = @constraint(m, [j=J,la=L_ac], pl[j,la] == Bl_ac*(θ[j,la[1]] - θ[j,la[2]])) # 'DC' Power Flow Approximation # θ should have a node as argument but l[x] is a node 
    con_PFap_dc = m.ext[:constraints][:con_PFap_dc] = @constraint(m, [j=J,ld=L_dc], pl[j,ld] == Bl_dc*(θ[j,ld[1]] - θ[j,ld[2]])) # 'DC' Power Flow Approximation # θ should have a node as argument but l[x] is a node 
    con_varlac1 = m.ext[:constraints][:con_varlac1] = @constraint(m, [j=J,la=L_ac], pl[j,la] <= varlac[la]) # Upper Limit Power Flow AC
    con_varlac2 = m.ext[:constraints][:con_varlac2] = @constraint(m, [j=J,la=L_ac], -varlac[la] <= pl[j,la]) # Lower Limit Power Flow AC
    con_varldc1 = m.ext[:constraints][:con_varldc1] = @constraint(m, [j=J,ld=L_dc], pl[j,ld] <= varldc[ld]) # Upper Limit Power Flow DC
    con_varldc2 = m.ext[:constraints][:con_varldc2] = @constraint(m, [j=J,ld=L_dc], -varldc[ld] <= pl[j,ld]) # Lower Limit Power Flow DC
    con_θb = m.ext[:constraints][:con_θb] = @constraint(m, [j=J,n=N], -pi <= θ[j,n] <= pi) # Bound Voltage angles
    con_θref = m.ext[:constraints][:con_θref] = @constraint(m, [j=J], θ[j,"Rome"] == 0.0) # Voltage angle at ref node
    con_DGl = m.ext[:constraints][:con_DGl] = @constraint(m, [i=ID, j=J, n=N], g[i,j,n] <= AF[i]*cap[i,n]) # Dispatchable generation limit
    con_SGl = m.ext[:constraints][:con_SGl] = @constraint(m, [i=["Solar"], j=J, n=N], g[i,j,n] <= AFS[Symbol(n)][j]*AF[i]*cap[i,n]) # Solar generation limit
    con_WGl = m.ext[:constraints][:con_WGl] = @constraint(m, [i=["WindOnshore", "WindOffshore"], j=J, n=N], g[i,j,n] <= AFW[Symbol(n)][j]*AF[i]*cap[i,n]) # Wind generation limit
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

# Plot transmission map
transmission = plot_transmission_network_italy()

# Plot side by side
plot(generation, transmission, layout=(1,2), size=(800,400))
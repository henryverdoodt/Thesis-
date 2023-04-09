## Practical Session II: ED & UC problems
# Author: Kenneth Bruninx, Sebastian Gonzato
# Last update: November 2, 2020

## Step 0: Activate environment - ensure consistency accross computers
using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
Pkg.instantiate() # If a Manifest.toml file exist in the current project, download all the packages declared in that manifest. Else, resolve a set of feasible packages from the Project.toml files and install them.

##  Step 1: input data
using CSV
using DataFrames
using YAML

data = YAML.load_file(joinpath(@__DIR__, "data.yaml"))
ts = CSV.read(joinpath(@__DIR__, "timeseries.csv"), DataFrame)
println("Done")

## Step 2: create model & pass data to model
using JuMP
using Gurobi
m = Model(optimizer_with_attributes(Gurobi.Optimizer))

# Step 2a: create sets
function define_sets!(m::Model, data::Dict, ts::DataFrame)
    m.ext[:sets] = Dict()

    # Time steps
    J = m.ext[:sets][:J] = 1:24

    # Dispatchable generators per type
    IDtype = m.ext[:sets][:IDtype] = [id for id in keys(data["dispatchableGenerators"])]
    # Dispatchable generators per unit
    ID = Array{Union{Nothing,String}}(nothing,0)
    for idtype in IDtype, i in 1:data["dispatchableGenerators"][idtype]["numberOfUnits"]
       ID = m.ext[:sets][:ID] = push!(ID,string(idtype,"_$(i)"))
    end

    # Variable generators
    IV = m.ext[:sets][:IV] = [iv for iv in keys(data["variableGenerators"])]

    # All variable and dispatchable generators, per unit
    I = m.ext[:sets][:I] = union(m.ext[:sets][:IV],m.ext[:sets][:ID])

    # return model
    return m
end

# Step 2b: add time series
function process_time_series_data!(m::Model, ts::DataFrame)
    # extract the relevant sets
    IV = m.ext[:sets][:IV]
    J = m.ext[:sets][:J]

    # create dictionary to store time series
    m.ext[:timeseries] = Dict()
    m.ext[:timeseries][:AF] = Dict()

    # add time series to dictionary
    m.ext[:timeseries][:D] = ts.Load[1:24]
    m.ext[:timeseries][:AF][IV[1]] = ts.LFW[1:24]
    m.ext[:timeseries][:AF][IV[2]] = ts.LFS[1:24]

    # return model
    return m
end

# step 2c: process input parameters
function process_parameters!(m::Model, data::Dict)
    # extract sets
    ID = m.ext[:sets][:ID]
    IV = m.ext[:sets][:IV]

    # Create parameter dictonary
    m.ext[:parameters] = Dict()

    # general parameters
    m.ext[:parameters][:VOLL] = data["valueOfLostLoad"] # VOLL
    m.ext[:parameters][:alphaCO2] = data["CO2Price"] # alphaCO2
    m.ext[:parameters][:timestep] = data["timeStepLength"] # timestep

    # parameters of dispatchable generators per unit
    d = data["dispatchableGenerators"]
    m.ext[:parameters][:MUT] = Dict(i => d[SubString(i,1:length(i)-2)]["minUpTime"] for i in ID)
    m.ext[:parameters][:GmaxD] = Dict(i => d[SubString(i,1:length(i)-2)]["maxPowerOutput"] for i in ID)
    m.ext[:parameters][:GminD] = Dict(i => d[SubString(i,1:length(i)-2)]["minStableOperatingPoint"] for i in ID) # Capacity of unit i
    m.ext[:parameters][:fuelcost] = Dict(i => d[SubString(i,1:length(i)-2)]["fuelcost"] for i in ID)
        
    #for problem EC
    m.ext[:parameters][:betaD] = Dict(i => (d[SubString(i,1:length(i)-2)]["fuelcost"]/d[SubString(i,1:length(i)-2)]["effmax"]) for i in ID) # Cost of generation of unit i
    m.ext[:parameters][:deltaD] = Dict(i => (d[SubString(i,1:length(i)-2)]["carbonintensity"]/d[SubString(i,1:length(i)-2)]["effmax"]) for i in ID) # CO2 intensity (adapted with eff max)

    # parameters of variable generators
    d = data["variableGenerators"]
    m.ext[:parameters][:GmaxV] = Dict(i => d[i]["installedCapacity"] for i in IV)
    m.ext[:parameters][:betaV] = Dict(i => 0 for i in IV) # cost of gen
    m.ext[:parameters][:deltaV] = Dict(i => 0 for i in IV) # CO2 intensity

    #for problem UC
    d = data["dispatchableGenerators"]
    m.ext[:parameters][:STE] = Dict(i => d[SubString(i,1:length(i)-2)]["startupenergy"] for i in ID)   

    m.ext[:parameters][:alphaD] = Dict(i => (d[SubString(i,1:length(i)-2)]["fuelcost"]/d[SubString(i,1:length(i)-2)]["effmin"]) for i in ID) # Cost of generation of unit i at min operating point
    m.ext[:parameters][:betaD] = Dict(i => (d[SubString(i,1:length(i)-2)]["fuelcost"]/d[SubString(i,1:length(i)-2)]["effmax"]) for i in ID) # Marginal Cost of generation of unit i

    m.ext[:parameters][:gammaD] = Dict(i => (d[SubString(i,1:length(i)-2)]["carbonintensity"]/d[SubString(i,1:length(i)-2)]["effmin"]) for i in ID) # CO2 intensity at min operating point
    m.ext[:parameters][:deltaD] = Dict(i => (d[SubString(i,1:length(i)-2)]["carbonintensity"]/d[SubString(i,1:length(i)-2)]["effmax"]) for i in ID) # Marginal CO2 intensity 

    #for problem UC with technical constraints
    d = data["dispatchableGenerators"]
    m.ext[:parameters][:SU] = Dict(i => d[SubString(i,1:length(i)-2)]["startupRate"]  for i in ID) # start-up rate of unit i, MW
    m.ext[:parameters][:SD] = Dict(i => d[SubString(i,1:length(i)-2)]["shutdownRate"] for i in ID) # shut-down rate of unit i, MW
    m.ext[:parameters][:RU] = Dict(i => (d[SubString(i,1:length(i)-2)]["rampUpRate"] * d[SubString(i,1:length(i)-2)]["maxPowerOutput"] * 60) for i in ID) # upward ramping rate of unit i, MW
    m.ext[:parameters][:RD] = Dict(i => (d[SubString(i,1:length(i)-2)]["rampDownRate"] * d[SubString(i,1:length(i)-2)]["maxPowerOutput"] * 60) for i in ID) # downward ramping rate of unit i, MW
    m.ext[:parameters][:MUT] = Dict(i => d[SubString(i,1:length(i)-2)]["minUpTime"] for i in ID)   # min uptime of unit i, h
    m.ext[:parameters][:MDT] = Dict(i => d[SubString(i,1:length(i)-2)]["minDownTime"] for i in ID) # min downtime of unit i, h

    # parameters storage
    d = data["Storage"]
    m.ext[:parameters][:Pmax] = d["Pmax"] # MW -> ̅D = ̅C = Pmax
    m.ext[:parameters][:Emax] = d["Emax"] # MWh -> ̅E
    m.ext[:parameters][:eff] = d["eff"] # Round-trip efficiency -> ηc = ηd = sqrt(eff)
    m.ext[:parameters][:Einit] = d["Einit"] # MWh, starting value state of charge -> e_0
    m.ext[:parameters][:Efinal] = d["Efinal"] # MWh, ending value state of charge -> e_J

    # return model
    return m
end

# call functions
define_sets!(m, data, ts)
process_time_series_data!(m, ts)
process_parameters!(m, data)

## Step 3: construct your model
function build_basic_ED_model!(m::Model)
    # Create m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    I = m.ext[:sets][:I] #Generators
    J = m.ext[:sets][:J] #TimeSteps
    ID = m.ext[:sets][:ID]
    IV = m.ext[:sets][:IV]

    # Extract time series data
    D = m.ext[:timeseries][:D]
    AF = m.ext[:timeseries][:AF]

    # Extract parameters
    VOLL = m.ext[:parameters][:VOLL]
    alphaCO2 = m.ext[:parameters][:alphaCO2]
    MUT = m.ext[:parameters][:MUT] 
    GmaxD = m.ext[:parameters][:GmaxD] 
    betaD = m.ext[:parameters][:betaD]
    deltaD = m.ext[:parameters][:deltaD]
    GmaxV = m.ext[:parameters][:GmaxV] 
    betaV = m.ext[:parameters][:betaV] 
    deltaV = m.ext[:parameters][:deltaV] 

    # create variables (g is already created for you, don't change this)
    g = m.ext[:variables][:g] = @variable(m, [i=I,j=J], lower_bound=0, base_name="generation")
    ens = m.ext[:variables][:ens] = @variable(m, [j=J], lower_bound=0, base_name="energy not served")
    fcd = m.ext[:variables][:fcd] = @variable(m, [i=ID,j=J], lower_bound=0, base_name="fuel cost")
    ccd = m.ext[:variables][:ccd] = @variable(m, [i=ID,j=J], lower_bound=0, base_name="emission cost")

    # Create affine expressions (= linear combinations of variables)
    #m.ext[:expressions][:emissioncost] = @expression(m, [i=ID,j=J], ccd[i,j] == alphaCO2*deltaD[i]*g[i,j])
    #m.ext[:expressions][:fuelcost] = @expression(m, [i=ID,j=J], fcd[i,j] == betaD[i]*g[i,j])
    # Issue with the expression Format
    # Calling it an Expression or an Constraint is only a way to define it doesnt change anything

    # Objective
    obj = m.ext[:objective] = @objective(m, Min, sum(fcd[i,j]+ccd[i,j]+VOLL*ens[j] for i in ID, j in J))
    
    
    # Constraints
    # Formulate constraints (market clearing constraint "con2a" is given - don't change this name, but add your energy not served variable and demand parameter)
    con2a = m.ext[:constraints][:con2a] = @constraint(m, [j=J], sum(g[i,j] for i in I) == D[j] - ens[j])
    con2c = m.ext[:constraints][:con2c] = @constraint(m, [j=J], 0 <= ens[j] <= D[j])
    con3a = m.ext[:constraints][:con3a] = @constraint(m, [i=IV, j=J], 0 <= g[i,j] <= AF[i][j]*GmaxV[i])
    con3b = m.ext[:constraints][:con3b] = @constraint(m, [i=ID, j=J], 0 <= g[i,j] <= GmaxD[i])
    # Create affine Expressions (= linear combinations of variables)
    exp1b = m.ext[:constraints][:exp1b] = @constraint(m, [i=ID, j=J], fcd[i,j] == betaD[i]*g[i,j])
    exp1c = m.ext[:constraints][:exp1c] = @constraint(m, [i=ID, j=J], ccd[i,j] == alphaCO2*deltaD[i]*g[i,j])

end

function build_basic_ED_storage_model!(m::Model)
    # Build the basic model
    build_basic_ED_model!(m)

    # Extract sets
    I = m.ext[:sets][:I] #Generators
    J = m.ext[:sets][:J] #TimeSteps
    ID = m.ext[:sets][:ID]
    IV = m.ext[:sets][:IV]

    # Extract time series data
    D = m.ext[:timeseries][:D]
    AF = m.ext[:timeseries][:AF]

    # Extract parameters
    Pmax = m.ext[:parameters][:Pmax] 
    Emax = m.ext[:parameters][:Emax] 
    eff = m.ext[:parameters][:eff] 
    Einit = m.ext[:parameters][:Einit] 
    Efinal = m.ext[:parameters][:Efinal]
    T = m.ext[:parameters][:timestep]


    # Extract variables
    g = m.ext[:variables][:g]
    ens = m.ext[:variables][:ens]

    # create variables (state-of-charge "e", charging "c" and discharging "d" are given, don't change these names)
    c =  m.ext[:variables][:c] = @variable(m, [j=J], lower_bound=0, upper_bound= Pmax, base_name="charging")
    d =  m.ext[:variables][:d] = @variable(m, [j=J], lower_bound=0, upper_bound= Pmax, base_name="discharging")
    e =  m.ext[:variables][:e] = @variable(m, [j=J], lower_bound=0, upper_bound= Emax, base_name="SOC")
    # Einit = m.ext[:variables][:Einit] = @variable(m, [j=J], lower_bound=0, upper_bound= Emax, base_name="Initial Energy")
    # Efinal = m.ext[:variables][:Efinal] = @variable(m, [j=J], lower_bound=0, upper_bound= Emax, base_name="Final Energy")

    # Remove market clearing constraint using "delete!"
    for x in J
        delete(m,m.ext[:constraints][:con2a][x])
    end

    # Add new market clearing constraint
    con2a = m.ext[:constraints][:con2a] = @constraint(m, [j=J], sum(g[i,j] for i in I) + d[j] - c[j] == D[j] - ens[j])

    # Add storage model with cyclic boundary conditions
    m.ext[:constraints][:con4a] = @constraint(m, [j=J[2:end-1]], e[j+1] - e[j] == T*sqrt(eff)*c[j] - T*(d[j]/sqrt(eff)))
    m.ext[:constraints][:con4b] = @constraint(m, [j=J[1]], e[j+1] - Einit == T*sqrt(eff)*c[j] - T*(d[j]/sqrt(eff)))
    m.ext[:constraints][:con4c] = @constraint(m, [j=J[end]], Efinal - e[j] == T*sqrt(eff)*c[j] - T*(d[j]/sqrt(eff)))

    return m
end

function build_basic_UC_model!(m::Model)
    # Build the basic model
    build_basic_ED_model!(m)

    # Extract sets
    I = m.ext[:sets][:I]
    ID = m.ext[:sets][:ID]
    IV = m.ext[:sets][:IV]

    # Extract time series
    J = m.ext[:sets][:J] # Time steps
    D = m.ext[:timeseries][:D]


    # Extract parameters
    alphaD = m.ext[:parameters][:alphaD]
    GminD = m.ext[:parameters][:GminD]
    GmaxD = m.ext[:parameters][:GmaxD]
    gammaD = m.ext[:parameters][:gammaD]
    alphaC02 = m.ext[:parameters][:alphaCO2]
    deltaD = m.ext[:parameters][:deltaD]
    STE = m.ext[:parameters][:STE]
    betaD = m.ext[:parameters][:betaD]
    VOLL = m.ext[:parameters][:VOLL]
    fuelcost = m.ext[:parameters][:fuelcost]

    #extract variables
    g = m.ext[:variables][:g]
    ens = m.ext[:variables][:ens]
    fcd = m.ext[:variables][:fcd]
    ccd = m.ext[:variables][:ccd]

    # Create variables (commitment status "z" and start-up "v" are given, don't change this)
    z = m.ext[:variables][:z] = @variable(m, [i=ID,j=J], binary=true, base_name="commitment")
    v = m.ext[:variables][:v] = @variable(m, [i=ID,j=J], binary=true, base_name="start_up")
    w = m.ext[:variables][:w] = @variable(m, [i=ID,j=J], binary=true, base_name="shut_down")
    scu = m.ext[:variables][:scu] = @variable(m, [i=ID,j=J],lower_bound = 0 ,base_name="start_up_cost")

    # Remove constraints & expressions that need to be adapted
    for x in J
        for y in ID
            delete(m,m.ext[:constraints][:exp1b][y,x])
            delete(m,m.ext[:constraints][:exp1c][y,x])
            delete(m,m.ext[:constraints][:con3b][y,x])
        end
    end
    # Create affine expressions (= linear combinations of variables)
    cond5a = m.ext[:constraints][:con5a] = @constraint(m, [i=ID,j=J], fcd[i,j] == alphaD[i]*z[i,j]+betaD[i]*(g[i,j]-z[i,j]*GminD[i]))
    cond5b = m.ext[:constraints][:con5b] = @constraint(m, [i=ID,j=J], ccd[i,j] == alphaC02*(z[i,j]*gammaD[i]+deltaD[i]*(g[i,j]-z[i,j]*GminD[i])))
    cond5c = m.ext[:constraints][:con5b] = @constraint(m, [i=ID,j=J], scu[i,j] == STE[i]*fuelcost[i]*v[i,j])

    # Formulate objective
    obj = m.ext[:objective] = @objective(m, Min, sum(fcd[i,j]+ccd[i,j]+scu[i,j]+VOLL*ens[j] for i in ID, j in J))

    # Formulate constraints
    cond5d = m.ext[:constraints][:cond5d] = @constraint(m, [i=ID,j=J], GminD[i]*z[i,j]<=g[i,j])
    cond5g = m.ext[:constraints][:cond5g] = @constraint(m, [i=ID,j=J], g[i,j]<=GmaxD[i]*z[i,j])
    cond5e = m.ext[:constraints][:cond5e] = @constraint(m, [i=ID,j=J[2:end]], z[i,j-1]-z[i,j]+v[i,j]-w[i,j]==0)
    cond5f = m.ext[:constraints][:cond5f] = @constraint(m, [i=ID,j=J[1]], 1-z[i,j]+v[i,j]-w[i,j]==0)

    return m
end

function build_basic_UC_technical_constraints_model!(m::Model)
    # Build the basic model
    build_basic_UC_model!(m)

    # Extract sets
    I = m.ext[:sets][:I]
    ID = m.ext[:sets][:ID]
    IV = m.ext[:sets][:IV]

    # Extract time series
    J = m.ext[:sets][:J] # Time steps
    D = m.ext[:timeseries][:D]

    # Extract parameters
    SU = m.ext[:parameters][:SU] 
    SD = m.ext[:parameters][:SD]
    RU = m.ext[:parameters][:RU] 
    RD = m.ext[:parameters][:RD] 
    MUT = m.ext[:parameters][:MUT]
    MDT = m.ext[:parameters][:MDT]

    # Extract variables
    g = m.ext[:variables][:g]
    z = m.ext[:variables][:z]
    v = m.ext[:variables][:v]
    w = m.ext[:variables][:w]

    # Remove constraints & expressions that need to be adapted
    for x in J
        for y in ID
            delete(m,m.ext[:constraints][:cond5d][y,x])
            delete(m,m.ext[:constraints][:cond5g][y,x])
        end
    end

    # Additional technical constraints (Minimum up-time constraint "con3d" is given, don't change this name)
    con3b2 = m.ext[:constraints][:con3b2] = @constraint(m, [i=ID,j=J[2:end]], g[i,j] - g[i,j-1] <= RU[i]*z[i,j] + ((SU[i]-RU[i])*v[i,j]))
    con3c2 = m.ext[:constraints][:con3c2] = @constraint(m, [i=ID,j=J[2:end]], g[i,j-1] - g[i,j] <= RD[i]*z[i,j] + (SD[i]*w[i,j]))
    con3d = m.ext[:constraints][:con3d] = Dict()
    for i in ID
        for j in J[MUT[i]:end]
            con3d[i,j] = @constraint(m,
                z[i,j] >= sum(v[i,j-jj] for jj in 0:MUT[i]-1)
            )
        end
    end
    con3e = m.ext[:constraints][:con3e] = Dict()
    for i in ID
        for j in J[MDT[i]:end]
            con3e[i,j] = @constraint(m,
                1 - z[i,j] >= sum(w[i,j-jj] for jj in 0:MDT[i]-1)
            )
        end
    end

    return m
end

# Build model: uncomment the relevant function-call, depending on which model you're running
# build_basic_ED_model!(m)
build_basic_ED_storage_model!(m)
# build_basic_UC_model!(m)
# build_basic_UC_technical_constraints_model!(m)

## Step 4: solve
optimize!(m)
# check termination status
print(
    """

    Termination status: $(termination_status(m))

    """
)
# check objective value
@show value(m.ext[:objective])
@show value(m.ext[:objective])/sum(m.ext[:timeseries][:D])

## Step 5: Visualization
using Plots
using StatsPlots

# sets
J = m.ext[:sets][:J]
I = m.ext[:sets][:I]
ID = m.ext[:sets][:ID]

# parameters
D = m.ext[:timeseries][:D]

# variables/expressions
Ivec = [i for  i in I]
g = value.(m.ext[:variables][:g])
gvec = [g[i,j] for  i in I, j in J]

if haskey(m.ext[:variables],:z)
    z = value.(m.ext[:variables][:z])
    uc = [z[i,j] for i in ID, j in J]
    uc = [zeros(2,24); uc] # plotting purposes -- added zero rows for wind and solar to make legend match
else
    λ = dual.(m.ext[:constraints][:con2a])
    λvec = [λ[j] for j in J]
end

# Check if storage is present in the model
if haskey(m.ext[:variables],:e)
    d = value.(m.ext[:variables][:d])
    c = value.(m.ext[:variables][:c])
    dcvec = [d[j]-c[j] for j in J]
    gvec = [transpose(dcvec); gvec]
    Ivec =["Storage"; Ivec]
end

if haskey(m.ext[:variables],:z)
# UC schedule
p1 = groupedbar(transpose(uc[:,:]),bar_position = :stack,xlabel="Timesteps [-]",ylabel="Commitment [# units]",label="",legend=:outertopright);
else
# electricity price price
p1 = plot(J,λvec, xlabel="Timesteps [-]",ylabel="λ [EUR/MWh]",label="",legend=:outertopright);
end
# dispatch
p2 = groupedbar(transpose(gvec[:,:]),bar_position = :stack,label=permutedims(Ivec),legend=:outertopright);
plot!(J,D, xlabel="Timesteps [-]",ylabel="Generation [MWh]", linewidth = 2,label="");
plot(p1, p2, layout = (1,2))
plot!(size=(900,400))

#=
if haskey(m.ext[:variables],:z) # unit commitment
    if haskey(m.ext[:constraints],:con3d) # technical constraints in the model
    savefig("BasicUCTechCon")
    else
    savefig("BasicUC")
    end
else
    if haskey(m.ext[:variables],:e) # storage in the model
    savefig("BasicEDStorage")
    else
    savefig("BasicED")
    end
end
=#

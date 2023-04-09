## Step 0: Activate environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.add("Ipopt")
Pkg.add("PowerModels")
Pkg.add("JuMP")
using PowerModels, Ipopt, JuMP

# Define solver
ipopt = optimizer_with_attributes(Ipopt.Optimizer)

## Step 1: input data
case_5 = "/Users/henryverdoodt/Documents/CODE/JULIA/OPES/PLC3_OPF/case5.m"
#case_5 = "case5.m"
## Here we are using the parser of Powermodels for convenience 
data = PowerModels.parse_file(case_5)

## Step 2: create the JuMP model & pass data to model
m = Model(ipopt)

# Step 2a: create sets
function define_sets!(m::Model, data::Dict)
    m.ext[:sets] = Dict()
    # ac nodes
    N = m.ext[:sets][:N] = 1:length(data["bus"])
    # ac branches
    B = m.ext[:sets][:B] = 1:length(data["branch"]) 
    # ac topology from side (i->j) and to side (j->i)
    B_ac_fr = m.ext[:sets][:B_ac_fr] = [(b, data["branch"][string(b)]["f_bus"], data["branch"][string(b)]["t_bus"]) for b in B] 
    B_ac_to = m.ext[:sets][:B_ac_to] = [(b, data["branch"][string(b)]["t_bus"], data["branch"][string(b)]["f_bus"]) for b in B]
    # Build union set of both sets above
    B_ac = m.ext[:sets][:B_ac] = [B_ac_fr; B_ac_to]
    # Branch connectivity to buses, e.g. which branches are connected to a certain node, used in nodal power balance equations
    bus_arcs = Dict((i, Tuple{Int,Int,Int}[]) for i in N)
    for (l,i,j) in B_ac
        push!(bus_arcs[i], (l,i,j))
    end
    B_arcs = m.ext[:sets][:B_arcs] = bus_arcs
    # Set of generators
    G = m.ext[:sets][:G] = 1:length(data["gen"])
    # Generator connectivity, e.g. which generators are connected to a certain node, used in nodal power balance equations
    bus_gens = Dict((i, Int[]) for i in N)
    for (g, gen) in data["gen"]
        push!(bus_gens[gen["gen_bus"]], parse(Int, g))
    end
    G_ac = m.ext[:sets][:G_ac] = bus_gens
    # Set of loads
    L = m.ext[:sets][:L] = 1:length(data["load"])
    # Load connectivity, e.g. which loads are connected to a certain node, used in nodal power balance equations
    bus_loads = Dict((i, Int[]) for i in N)
    for (l, load) in data["load"]
        push!(bus_loads[load["load_bus"]], parse(Int, l))
    end
    L_ac = m.ext[:sets][:L_ac] = bus_loads 
    return m
end

# step 2b: process input parameters
function process_parameters!(m::Model, data::Dict)
    # extract sets
    N = m.ext[:sets][:N]
    B = m.ext[:sets][:B]
    G = m.ext[:sets][:G]
    L = m.ext[:sets][:L]

    # Create parameter dictionary
    m.ext[:parameters] = Dict()

    # Bus parameters
    vmmin = m.ext[:parameters][:vmmin] = Dict(i => data["bus"][string(i)]["vmin"] for i in N) # minimum voltage magnitude
    vmmax = m.ext[:parameters][:vmmax] = Dict(i => data["bus"][string(i)]["vmax"] for i in N) # maximum voltage magnitude
    vamin = m.ext[:parameters][:vamin] = Dict(i => -pi for i in N) # Arbitrary limit of -pi for minimum bus voltage angle
    vamax = m.ext[:parameters][:vamax] = Dict(i =>  pi for i in N) # Arbitrary limit of  pi for maximum bus voltage angle

    # Branch parameters
    r = m.ext[:parameters][:br] = Dict(b => data["branch"][string(b)]["br_r"] for b in B) # branch resistance
    x = m.ext[:parameters][:bx] = Dict(b => data["branch"][string(b)]["br_x"] for b in B) # branch reactance
    gb =  m.ext[:parameters][:gb] = Dict(b => real(1 / (data["branch"][string(b)]["br_r"] + data["branch"][string(b)]["br_x"]im)) for b in B) # branch series conductance
    bb =  m.ext[:parameters][:bb] = Dict(b => imag(1 / (data["branch"][string(b)]["br_r"] + data["branch"][string(b)]["br_x"]im)) for b in B) # branch series admittance
    gfr = m.ext[:parameters][:gb_sh_fr] = Dict(b => data["branch"][string(b)]["g_fr"] for b in B) # branch shunt conductance from side i -> j
    bfr = m.ext[:parameters][:bb_sh_fr] = Dict(b => data["branch"][string(b)]["b_fr"] for b in B) # branch shunt susceptance from side i -> j
    gto = m.ext[:parameters][:gb_sh_to] = Dict(b => data["branch"][string(b)]["g_to"] for b in B) # branch shunt conductance to side j -> i
    bto = m.ext[:parameters][:bb_sh_to] = Dict(b => data["branch"][string(b)]["b_to"] for b in B) # branch shunt susceptance to side  j -> j
    smax = m.ext[:parameters][:smax] = Dict(b => data["branch"][string(b)]["rate_a"] for b in B)  # branch rated power in pu
    angmin = m.ext[:parameters][:angmin] = Dict(b => data["branch"][string(b)]["angmin"] for b in B) # minimum voltage angle difference over branch
    angmax = m.ext[:parameters][:angmax] = Dict(b => data["branch"][string(b)]["angmax"] for b in B) # maximum voltage angle difference over branch
    
    # load parameters: Assuming a fixed demand!
    pd = m.ext[:parameters][:pd] = Dict(l => data["load"][string(l)]["pd"] for l in L)  # active power demand in pu
    qd = m.ext[:parameters][:qd] = Dict(l => data["load"][string(l)]["qd"] for l in L)  # reactive power demand in pu
    
    # generator parameters:
    pmax = m.ext[:parameters][:pmax] = Dict(g => data["gen"][string(g)]["pmax"] for g in G)  # maximum active power in pu
    pmin = m.ext[:parameters][:pmin] = Dict(g => data["gen"][string(g)]["pmin"] for g in G)  # minimum active power in pu
    qmax = m.ext[:parameters][:qmax] = Dict(g => data["gen"][string(g)]["qmax"] for g in G)  # maximum reactive power in pu
    qmin = m.ext[:parameters][:qmin] = Dict(g => data["gen"][string(g)]["qmin"] for g in G)  # minimum reactive power in pu
    ag = m.ext[:parameters][:ag] = Dict(g => data["gen"][string(g)]["cost"][2] for g in G)   # auxiliary generation cost
    bg = m.ext[:parameters][:bg] = Dict(g => data["gen"][string(g)]["cost"][1] for g in G)   # linear generation cost
    # return model
    return m
end

# call functions
define_sets!(m, data)
process_parameters!(m, data)

## Step 3: construct your model
function build_model!(m::Model)
    # Create m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    N = m.ext[:sets][:N]
    B = m.ext[:sets][:B]
    B_ac_fr = m.ext[:sets][:B_ac_fr]
    B_ac_to = m.ext[:sets][:B_ac_to]
    G = m.ext[:sets][:G]
    G_ac = m.ext[:sets][:G_ac]
    L = m.ext[:sets][:L]
    L_ac = m.ext[:sets][:L_ac]
    B_ac = m.ext[:sets][:B_ac]
    B_arcs = m.ext[:sets][:B_arcs]

    # Extract parameters
    vmmin = m.ext[:parameters][:vmmin]
    vmmax = m.ext[:parameters][:vmmax]
    vamin = m.ext[:parameters][:vamin]
    vamax = m.ext[:parameters][:vamax]
    gb =  m.ext[:parameters][:gb]
    bb =  m.ext[:parameters][:bb] 
    gfr = m.ext[:parameters][:gb_sh_fr]
    bfr = m.ext[:parameters][:bb_sh_fr]
    gto = m.ext[:parameters][:gb_sh_to]
    bto = m.ext[:parameters][:bb_sh_to]
    smax = m.ext[:parameters][:smax]
    angmin = m.ext[:parameters][:angmin]
    angmax = m.ext[:parameters][:angmax]
    pd = m.ext[:parameters][:pd]
    qd = m.ext[:parameters][:qd]
    pmax = m.ext[:parameters][:pmax]
    pmin = m.ext[:parameters][:pmin]
    qmax = m.ext[:parameters][:qmax]
    qmin = m.ext[:parameters][:qmin]
    ag = m.ext[:parameters][:ag]
    bg = m.ext[:parameters][:bg]

    ## create variables 
    # bus variables
    vm = m.ext[:variables][:vm] = @variable(m, [i=N], lower_bound = vmmin[i], upper_bound = vmmax[i], base_name = "vm") # voltage magnitude
    va = m.ext[:variables][:va] = @variable(m, [i=N], lower_bound = vamin[i], upper_bound = vamax[i], base_name = "va") # voltage angle

    # generator variables
    pg = m.ext[:variables][:pg] = @variable(m, [g=G], lower_bound = pmin[g], upper_bound = pmax[g], base_name = "pg") # active power
    qg = m.ext[:variables][:qg] = @variable(m, [g=G], lower_bound = qmin[g], upper_bound = qmax[g], base_name = "qg") # reactive power

    # branch variables
    pb = m.ext[:variables][:pb] = @variable(m, [(b,i,j) in B_ac], lower_bound = -smax[b], upper_bound = smax[b], base_name = "pb") # from side active power flow (i->j)
    qb = m.ext[:variables][:qb] = @variable(m, [(b,i,j) in B_ac], lower_bound = -smax[b], upper_bound = smax[b], base_name = "qb") # from side reactive power flow (i->j)
     
    ## Objective
    m.ext[:objective] = @objective(m, Min, sum([ag[g] + bg[g]*pg[g] for g in G]))

    # Power flow constraints in from and to direction
    m.ext[:constraints][:pbij] = @NLconstraint(m, [(b,i,j) = B_ac_fr], pb[(b, i, j)] ==  (gb[b] + gfr[b])*vm[i]^2 - (gb[b] * vm[i] * vm[j] * cos(va[i] - va[j])) - (bb[b] * vm[i] * vm[j] * sin(va[i] - va[j]))) # active power i to j
    m.ext[:constraints][:qbij] = @NLconstraint(m, [(b,i,j) = B_ac_fr], qb[(b, i, j)] == -(bb[b] + bfr[b])*vm[i]^2 + (bb[b] * vm[i] * vm[j] * cos(va[i] - va[j])) - (gb[b] * vm[i] * vm[j] * sin(va[i] - va[j]))) # reactive power i to j
    m.ext[:constraints][:pbji] = @NLconstraint(m, [(b,j,i) = B_ac_to], pb[(b, j, i)] ==  (gb[b] + gto[b])*vm[j]^2 - (gb[b] * vm[j] * vm[i] * cos(va[j] - va[i])) - (bb[b] * vm[j] * vm[i] * sin(va[j] - va[i]))) # active power j to i
    m.ext[:constraints][:qbji] = @NLconstraint(m, [(b,j,i) = B_ac_to], qb[(b, j, i)] == -(bb[b] + bto[b])*vm[j]^2 + (bb[b] * vm[j] * vm[i] * cos(va[j] - va[i])) - (gb[b] * vm[j] * vm[i] * sin(va[j] - va[i]))) # reactive power j to i

    # Thermal limits for the branches
    m.ext[:constraints][:sij] = @NLconstraint(m, [(b,i,j) = B_ac_fr], pb[(b, i, j)]^2 + qb[(b, i, j)]^2 <= smax[b]^2)
    m.ext[:constraints][:sji] = @NLconstraint(m, [(b,j,i) = B_ac_to], pb[(b, j, i)]^2 + qb[(b, j, i)]^2 <= smax[b]^2)

    # Branch angle limits
    m.ext[:constraints][:thetaij] = @constraint(m, [(b,i,j) = B_ac_fr], va[i] - va[j] <= angmax[b])
    m.ext[:constraints][:thetaji] = @constraint(m, [(b,i,j) = B_ac_fr], va[i] - va[j] >= angmin[b])
    m.ext[:constraints][:thetaij] = @constraint(m, [(b,j,i) = B_ac_to], va[j] - va[i] <= angmax[b])
    m.ext[:constraints][:thetaji] = @constraint(m, [(b,j,i) = B_ac_to], va[j] - va[i] >= angmin[b])

    # Kirchhoff's current law, e.g., nodal power balance
    m.ext[:constraints][:p_balance] = @constraint(m, [i in N], sum(pg[g] for g in G_ac[i]) - sum(pd[l] for l in L_ac[i]) == sum(pb[(b,i,j)] for (b,i,j) in B_arcs[i]) )
    m.ext[:constraints][:q_balance] = @constraint(m, [i in N], sum(qg[g] for g in G_ac[i]) - sum(qd[l] for l in L_ac[i]) == sum(qb[(b,i,j)] for (b,i,j) in B_arcs[i]) )

    # Voltage angle on reference bus = 0, reference bus is bus 4 in this case
    m.ext[:constraints][:varef] = @constraint(m, va[4] == 0)
    
    return m 
end

# Build model
build_model!(m)

# Optimize model
optimize!(m)

# Print solution for inspection
solution_summary(m, verbose=true)


######### Validation is key! Validate solution with PowerModels
result = PowerModels.solve_opf(case_5, ACPPowerModel, ipopt)
# Print PowerModels solution for inspection
print_summary(result["solution"])
###################################
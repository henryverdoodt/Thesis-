## Master Thesis: Analyse the impact of climate change on a G&TEP of EU Power System
# Author: GPT
# Last update: March 24, 2023
using Pkg
Pkg.activate(@__DIR__)
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
Pkg.add("GLPK")

using CSV, YAML, JuMP, GLPK, DataFrames, Distributions, Gurobi, Images, Plots, PolyChaos, InteractiveUtils

# Transmission Network Italy
begin
	img 	= load("/Users/henryverdoodt/Documents/CODE/JULIA/OPES/PLC_uncertainty/Italy.png");
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

# Load input data
demand = CSV.read(joinpath(@__DIR__, "case_6_demand_10.csv"), DataFrame)
wind_cf = CSV.read(joinpath(@__DIR__, "case_6_wind_10.csv"), DataFrame)
pv_cf = CSV.read(joinpath(@__DIR__, "case_6_PV_10.csv"), DataFrame)
data = YAML.load_file(joinpath(@__DIR__, "overview_data.yaml"))

# Define the sets
N = Set(["Milaan", "Florence", "Rome", "Bari", "Palermo", "Sardinia"]); # Nodes
L_ac = Set([("Florence","Rome"), ("Rome","Bari"), ("Milaan","Florence"), ("Milaan","Rome"), ("Milaan","Bari"), ("Florence","Bari")]); # Candidate AC Lines
L_dc = Set([("Rome","Palermo"), ("Sardinia","Palermo"), ("Bari", "Palermo"), ("Florence","Sardinia"), ("Rome","Sardinia")]); # Candidate DC Lines
L = union(L_ac, L_dc); # Union of Lines
T = Set(collect(1:8760)) # Timesteps
I = 

## Define the parameters
begin
	M = 100.0
	θmax = 0.4
end;

# Deterministic Model
m = Model(Gurobi.Optimizer)

# Define decision variables
@variable(m, g[node in 1:6, tech in keys(data["generation"])], lowerbound=0)
@variable(m, t[(i,j) in [("Florence","Rome"), ("Rome","Bari"), ("Milaan","Florence"), ("Milaan","Rome"),
                         ("Milaan","Bari"), ("Florence","Bari"), ("Rome","Palermo"), ("Sardinia","Palermo"), 
                         ("Bari", "Palermo"), ("Florence","Sardinia"), ("Rome","Sardinia")]])

# Define objective function
@objective(m, Min, sum(g[node, tech] * data["generation"][tech]["cost"] for node in 1:6, tech in keys(data["generation"])) 
                        + sum(t[(i,j)] * data["transmission"]["cost"] for (i,j) in [("Florence","Rome"), ("Rome","Bari"), 
                        ("Milaan","Florence"), ("Milaan","Rome"), ("Milaan","Bari"), ("Florence","Bari"), ("Rome","Palermo"), 
                        ("Sardinia","Palermo"), ("Bari", "Palermo"), ("Florence","Sardinia"), ("Rome","Sardinia")]))

# Define constraints
for node in 1:6
    @constraint(m, sum(g[node,tech] for tech in keys(data["generation"])) == demand[node,:][2:end] * 1e6)
end
for tech in keys(data["generation"])
    for node in 1:6
        @constraint(m, g[node,tech] <= data["generation"][tech]["capacity"])
    end
end
for (i,j) in [("Florence","Rome"), ("Rome","Bari"), ("Milaan","Florence"), ("Milaan","Rome"), ("Milaan","Bari"), ("Florence","Bari"), ("Rome","Palermo"), ("Sardinia","Palermo"), ("Bari", "Palermo"), ("Florence","Sardinia"), ("Rome","Sardinia")]
    @constraint(m, t[(i,j)] <= data["transmission"]["capacity"])
end

# Solve optimization problem
optimize!(m)

# Print optimal solution
println("Optimal generation capacity:")
for node in 1:6
    for tech in keys(data["generation"])
        println("Node $node, $tech: ", value(g[node,tech]))
    end
end
println("Optimal transmission capacity:")
for (i,j) in [("Florence","Rome"), ("Rome","Bari"), ("Milaan","Florence"), ("Milaan","Rome"), ("Milaan","Bari"), ("Florence","Bari"), ("Rome","Palermo"), ("Sardinia","Palermo"), ("Bari", "Palermo"), ("Florence","Sardinia"), ("Rome","Sardinia")]
    println("Line ($i,$j): ", value(t[(i,j)]))
end
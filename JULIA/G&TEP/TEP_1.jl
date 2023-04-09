### A Pluto.jl notebook ###

# v0.17.3
using Pkg
Pkg.add("Markdown")
Pkg.add("InteractiveUtils")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Distributions")
Pkg.add("Gurobi")
Pkg.add("Images")
Pkg.add("JuMP")
Pkg.add("Plots")
Pkg.add("PolyChaos")

using Markdown
using InteractiveUtils
using CSV, DataFrames, Distributions, Gurobi, Images, JuMP, Plots, PolyChaos

# Lecture 11 - Optimisation under uncertainty - Transmission expansion planning

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

## Define the sets

### Nodes
N = Set(["Milaan", "Florence", "Rome", "Bari", "Palermo", "Sardinia"]);

### Existing Lines
#Bᵃᶜ = Set([("Milaan","Florence"), ("Milaan","Rome"), ("Milaan","Bari"), ("Florence","Bari")]); # Existing AC Lines
#Bᵈᶜ = Set([("Florence","Sardinia"), ("Rome","Sardinia")]); # Existing DC Lines

### Candidate Lines
Cᵃᶜ = Set([("Florence","Rome"), ("Rome","Bari"), ("Milaan","Florence"), ("Milaan","Rome"), ("Milaan","Bari"), ("Florence","Bari")]);
Cᵈᶜ = Set([("Rome","Palermo"), ("Sardinia","Palermo"), ("Bari", "Palermo"), ("Florence","Sardinia"), ("Rome","Sardinia")]);

### Union subset for lines
#B = union(Bᵃᶜ, Bᵈᶜ);
C = union(Cᵃᶜ, Cᵈᶜ);
#L = union(C, C);

## Define the parameters

begin
	M = 100.0
	θmax = 0.4
end;

### Conventional Generators
Pgmax = Dict("Florence" => 4.0, "Rome" => 4.0);
Cg = Dict("Florence" => 8760 * 10 * 100 * 50,   # Cost: hours; years; MVA; €/MVA (gas)
		  "Rome" => 8760 * 10 * 100 * 20);	    # Cost: hours; years; MVA; €/MVA (nuclear)

### Load
Pd_dm = Dict("Milaan" => 0.5, "Florence" => 0.55, "Rome" => 0.6, "Bari" => 0.65, "Sardinia" => 0.7, "Palermo" => 0.2); # equivalent model: expectancy of future load

### Lines

Plmax = Dict(# existing ac-lines                                         # max power flow and min power flow in the other direction
			 ("Milaan","Florence") => 1.0, ("Milaan","Rome") => 0.8, 
			 ("Milaan","Bari") => 0.5, ("Florence","Bari") => 1.0, 
			 # existing dc-lines
			 ("Florence","Sardinia") => 1.0, ("Rome","Sardinia") => 1.0,
			 # candidate ac-lines
			 ("Florence","Rome") => 1.0, ("Rome","Bari") => 1.0,
			 # candidate dc-lines
			 ("Rome","Palermo") => 2.3, ("Sardinia","Palermo") => 2.3, 
			 ("Bari", "Palermo") => 2.3);

Bl = Dict(  # existing ac-lines  										# inverse of X of line, Susceptance
			("Milaan","Florence") => 2.5, ("Milaan","Rome") => 2.0, 
			("Milaan","Bari") => 5.0, ("Florence","Bari") => 2.5,
			# candidate ac-lines
			("Florence","Rome") => 5.0, ("Rome","Bari") => 5.0);

Cc = Dict(  # candidate ac-lines										# Cost associate to build the line: (lenght &/OR MW) * (€/MW)
			("Florence","Rome") => 126 * 4.0e6, ("Rome","Bari") => 211 * 4.0e6,
			# candidate dc-lines
			("Rome","Palermo") => 233 * 1.7e6, 							
			("Sardinia","Palermo") => 244 * 1.7e6, 
			("Bari", "Palermo") => 237 * 1.7e6);


## Deterministic Model, i.e., certainty equivalent problem where Pl = 𝔼[ξ] 
dm = Model(Gurobi.Optimizer)

### Variables
@variable(dm, -θmax <= θ[n ∈ N] <= θmax); # nodal voltage angle
@variable(dm, Pl[l ∈ L]); # active power flow over branches
@variable(dm, Pg[n ∈ N; n ∈ keys(Pgmax)]); # active power from conventional generators
@variable(dm, x[c ∈ C], Bin); # investment decision for candidate lines

### Constraints

# KCL
@constraint(dm, [n ∈ N], sum(Pg[m] for m in [n] if m ∈ keys(Pgmax); init=0) 
						+ sum(Pl[c] for c in C if c[1] == n; init=0) 
						== 
						Pd_dm[n]
						+ sum(Pl[c] for c in C if c[2] == n; init=0));

# Reference voltage angle
@constraint(dm, θ["Rome"] == 0.0);

# Generator limit
begin
	@constraint(dm, [n ∈ N; n ∈ keys(Pgmax)], 0.0 <= Pg[n])
	@constraint(dm, [n ∈ N; n ∈ keys(Pgmax)], Pg[n] <= Pgmax[n])
end;

# Angle difference over all lines
begin
	@constraint(dm, [c ∈ C], -θmax <= θ[c[1]] - θ[c[2]])
	@constraint(dm, [c ∈ C], θ[c[1]] - θ[c[2]] <= θmax)
end;

# Active power flow limit over all lines
begin
	#@constraint(dm, [b ∈ B], -Plmax[b] <= Pl[b])
	#@constraint(dm, [b ∈ B], Pl[b] <= Plmax[b])
	@constraint(dm, [c ∈ C], -x[c] * Plmax[c] <= Pl[c])
	@constraint(dm, [c ∈ C], Pl[c] <= x[c] * Plmax[c])
end;

# Active power flow over all ac lines
begin
	#@constraint(dm, [b ∈ Bᵃᶜ], Pl[b] == Bl[b] * (θ[b[1]] - θ[b[2]]))
	@constraint(dm, [c ∈ Cᵃᶜ], Pl[c] >= Bl[c] * (θ[c[1]] - θ[c[2]]) - (1-x[c]) * M)
	@constraint(dm, [c ∈ Cᵃᶜ], Pl[c] <= Bl[c] * (θ[c[1]] - θ[c[2]]) + (1-x[c]) * M)
end;

### Objective
@objective(dm, Min,  sum(Cc[c] * x[c] for c ∈ C; init=0) 
					+ sum(Cg[n] * Pg[n] for n ∈ N if n ∈ keys(Pgmax); init=0));

optimize!(dm)

### Result for the deterministic model

begin
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
	plot!([0.0,0.0],[0.0,0.0], line=("blue", 3, :dash), label="proposed ac-lines")
	plot!([0.0,0.0],[0.0,0.0], line=("yellow", 3, :dash), label="proposed dc-lines")
	if value.(x)[("Florence","Rome")] >= 0.5
	plot!([cities["Florence"][1],cities["Rome"][1]],
		  [cities["Florence"][2],cities["Rome"][2]],
		  line=("blue", 3, :dash), label="proposed ac-lines")
	end
	if value.(x)[("Rome","Bari")] >= 0.5
	plot!([cities["Bari"][1],cities["Rome"][1]],
		  [cities["Bari"][2],cities["Rome"][2]],
		  line=("blue", 3, :dash), label=:none)
	end
	if value.(x)[("Sardinia","Palermo")] >= 0.5
	plot!([cities["Palermo"][1],cities["Sardinia"][1]],
		  [cities["Palermo"][2],cities["Sardinia"][2]],
		  line=("yellow", 3, :dash), label="proposed dc-lines")
	end
	if value.(x)[("Rome","Palermo")] >= 0.5
		plot!([cities["Palermo"][1],cities["Rome"][1]],
			  [cities["Palermo"][2],cities["Rome"][2]],
			  line=("yellow", 3, :dash), label=:none)
	end
	if value.(x)[("Bari","Palermo")] >= 0.5
		plot!([cities["Bari"][1],cities["Palermo"][1]],
			  [cities["Bari"][2],cities["Palermo"][2]],
			  line=("yellow", 3, :dash), label=:none)
	end
	scatter!([nc[1] for nc in values(cities)],
			 [nc[2] for nc in values(cities)],
			 color="orange", label=:none)
end


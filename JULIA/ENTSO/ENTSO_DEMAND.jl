## Step 0: Activate environment - ensure consistency accross computers
using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
Pkg.instantiate() # If a Manifest.toml file exist in the current project, download all the packages declared in that manifest. Else, resolve a set of feasible packages from the Project.toml files and install them.
#Pkg.add("CSV")
#Pkg.add("DataFrames")
#Pkg.add("Plots")

##  Step 1: input data
using CSV
using DataFrames
using Plots

data = CSV.read(joinpath(@__DIR__, "NT2040_Demand_CY1984.csv"), DataFrame; header=false)

# Extract year, month, day, and hour data
year = data[!, 1]
month = data[!, 2]
day = data[!, 3]
hour = data[!, 4]


println(size(data, 1)) 
println(size(data, 2)) 
println(typeof(parse(Float64,data[4,6])))


#println(sum(parse(Float64, data[3, 5:8])))

# Calculate total demand for each hour
#total_demand = [sum(parse(Float64, data[i, 5:(size(data, 2))])) for i in 2:size(data, 1)]
#total_demand = [sum(data[i, 5:(size(data, 2))]) for i in 1:size(data, 1)]
row1 = (sum(data[1,i]) for i in 5:60)

total_demand = [data[i, sum(data[i,j])] for i in 1:size(data, 1), j in 5:size(data, 2)]

# Calculate total demand for each hour
#total_demand = [sum(data[i, 5:60]) for i in 2:size(data, 1)]

# Create timestamp vector
timestamp = [DateTime(year[i], month[i], day[i], hour[i], 0, 0) for i in 1:size(data, 1)]

# Plot demand vs time
plot(timestamp, total_demand, xlabel="Time", ylabel="Electricity Demand", legend=false)





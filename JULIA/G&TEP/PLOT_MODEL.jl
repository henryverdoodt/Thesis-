## Master Thesis: Analyse the impact of climate change on a G&TEP of EU Power System
# Author: Henry Verdoodt
# Last update: April 4, 2023

#############################################################################################################################################################
###############################                                      VISUALIZATION                                            ###############################
#############################################################################################################################################################

using CSV, YAML, JuMP, DataFrames, Distributions, Gurobi, Images, Plots, PolyChaos, InteractiveUtils, StatsPlots


# Dictionary to store the generators capacity for each node in MW
function get_generators_capacity(m::Model, I::Vector{String}, N::Vector{String})
    gen_dict = Dict{String, Dict{String, Float64}}()  
    for node in N
        gen_dict[node] = Dict((gen => 0.0 for gen in I)...)
        for unit in I
            gen_dict[node][unit] = value.(m.ext[:variables][:cap][unit, node])
        end
    end
    return gen_dict
end

# Create a dictionaries to store the AC transmission capacity
function get_ac_transmission_capacity(m::Model, L_ac::Vector{Vector{String}})
    trans_ac_dict = Dict(la => value.(m.ext[:variables][:varlac][la]) for la in L_ac)
    return trans_ac_dict
end

# Create a dictionaries to store the AC and DC transmission capacity
function get_dc_transmission_capacity(m::Model, L_dc::Vector{Vector{String}})
    trans_dc_dict = Dict(ld => value.(m.ext[:variables][:varldc][ld]) for ld in L_dc)
    return trans_dc_dict
end

# Generation per type given in GW
function matrix_generators_data(gen_dict::Dict{String, Dict{String, Float64}})
    # Extract the node names and the generating unit types
    nodes = collect(keys(gen_dict))
    unit_types = collect(keys(gen_dict[nodes[1]]))

    # Create the data matrix
    data = zeros(length(nodes), length(unit_types))
    for (i, node) in enumerate(nodes)
        for (j, unit_type) in enumerate(unit_types)
            data[i, j] = gen_dict[node][unit_type]/1000
        end
    end
    #cols_to_remove = findall(all.(x -> x == 0.0, eachcol(data)))
    #data = data[:, setdiff(1:size(data, 2), cols_to_remove)]
    return data
end

Generator_colors = Dict("Biomass" => :green, "WindOffshore" => :blue, "CCGT_new" => :orange, "WindOnshore" => :lightblue, "Nuclear" => :purple, "Solar" => :gold, "OCGT" => :red, "ICE" => :gray)
Generator_colors = [:green, :lightblue, :orange, :blue, :purple, :gold, :red, :gray]    # [:green, :blue, :lightblue, :gold]  
Generator_labels = ["Biomass" "WindOnshore" "CCGT_new" "WindOffshore" "Nuclear" "Solar" "OCGT" "ICE"]   # ["Biomass" "WindOffshore" "WindOnshore" "Solar"]  

function plot_generator_capacities(data_matrix::Matrix{Float64}, nodes::Vector{String}, Generator_colors::Vector{Symbol}, Generator_labels::Matrix{String})
    # Plot generation stacked bar chart
    groupedbar(data_matrix,
        bar_position = :stack,
        bar_width=0.5,
        bar_edges=true,
        ylabel= "Capacity (MW)",
        xticks=(1:length(nodes), nodes),
        label= Generator_labels,
        color_palette= Generator_colors,
        title="Installed Generation Capacity (GW)")
end


# Transmission Network West EU 
countries_coord = Dict("ES" => (380, 1400), "FR" => (640, 1140), "BE" => (722, 960), "DE" => (880, 960), 
                 "NL" => (745, 885), "DK" => (866, 720), "NO" => (890, 482), "UK" => (553, 835), "PT" => (210, 1380), 
                 "IT" => (930, 1330), "CH" => (800,1140), "AT" => (1000,1120), "IE" => (380,780), "SE" => (1000, 550), "FI" => (1230, 400))

#countries_coord = Dict("ES" => (380, 1400), "FR" => (640, 1140), "BE" => (722, 960), "DE" => (880, 960), "NL" => (745, 885), "DK" => (866, 720), "NO" => (890, 482), "UK" => (553, 835))

function plot_transmission_network(dict1::Dict{Vector{String}, Float64}, dict2::Dict{Vector{String}, Float64}, countries::Dict{String, Tuple{Int64, Int64}})
    img 	= load("/Users/henryverdoodt/Documents/CODE/IMAGES/Other/europe_map.jpeg")

    plot(img, axis=([], false), title="Transmission Lines")

    for (from, to) in keys(dict1)
        plot!([countries[from][1], countries[to][1]],
                [countries[from][2], countries[to][2]],
                color="blue", linewidth= 2, label=:none) 
    end

    for (from, to) in keys(dict2)
        plot!([countries[from][1], countries[to][1]],
                [countries[from][2], countries[to][2]],
                color="red", linewidth= 2, label=:none) 
    end
 
    scatter!([nc[1] for nc in values(countries)],
             [nc[2] for nc in values(countries)],
             markersize=3, color="black", label=:none)
end

function plot_transmission_needed(dict1::Dict{Vector{String}, Float64}, dict2::Dict{Vector{String}, Float64}, countries::Dict{String, Tuple{Int64, Int64}})
    img 	= load("/Users/henryverdoodt/Documents/CODE/IMAGES/Other/europe_map.jpeg")

    plot(img, axis=([], false), title="Transmission Lines")

    for (from, to) in keys(dict1)
        if dict1[[from, to]] != 0.0
            GW = round(dict1[[from, to]]/1000, digits=2)
            plot!([countries[from][1], countries[to][1]],
                  [countries[from][2], countries[to][2]],
                  color="blue", linewidth= ((dict1[[from, to]]/maximum(values(dict1)))*4 + 2), label="$from - $to : $(GW) GW", legend=:topright) 
        end
    end

    for (from, to) in keys(dict2)
        if dict2[[from, to]] != 0.0
            GW = round(dict2[[from, to]]/1000, digits=2)
            plot!([countries[from][1], countries[to][1]],
                  [countries[from][2], countries[to][2]],
                  color="red", linewidth= ((dict2[[from, to]]/maximum(values(dict1)))*4 + 2), label="$from - $to : $(GW) GW", legend=:topright) 
        end
    end
 
    scatter!([nc[1] for nc in values(countries)],
             [nc[2] for nc in values(countries)],
             markersize=3, color="black", label=:none)
end

# Heatmap with transmission capacities
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



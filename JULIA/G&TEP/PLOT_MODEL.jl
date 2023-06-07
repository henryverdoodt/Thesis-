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

function get_generators_capacity_storage(m::Model, I::Vector{String}, N::Vector{String}, S::Vector{String})
    gen_dict = Dict{String, Dict{String, Float64}}()  
    for node in N
        gen_dict[node] = Dict((gen => 0.0 for gen in cat(I, S, dims=1))...)
        for unit in I
            gen_dict[node][unit] = value.(m.ext[:variables][:cap][unit, node])
        end
        gen_dict[node]["PtH"] = value.(m.ext[:variables][:cap_PtH][node])
        gen_dict[node]["OCGT_H"] = value.(m.ext[:variables][:cap_OCGT][node])
    end
    return gen_dict
end

function get_generators_capacity_storage_fixed_network(m::Model, I::Vector{String}, IV::Vector{String}, ID::Vector{String}, N::Vector{String})
    gen_dict = Dict{String, Dict{String, Float64}}()  
    for node in N
        gen_dict[node] = Dict((gen => 0.0 for gen in I)...)
        for unit in IV
            gen_dict[node][unit] = value.(m.ext[:variables][:cap][unit, node])
        end
        for unit in ID
            gen_dict[node][unit] = value.(m.ext[:parameters][:cap_exist][node][unit])
        end
    end
    return gen_dict
end

function get_generators_capacity_storage_years(m::Model, I::Vector{String}, N::Vector{String}, S::Vector{String}, year::Int64)
    gen_dict = Dict{String, Dict{String, Dict{Int64, Float64}}}()  
    for node in N
        gen_dict[node] = Dict((gen => Dict{Int64, Float64}() for gen in cat(I, S, dims=1))...)
        for unit in I
            gen_dict[node][unit][year] = value.(m.ext[:variables][:cap][unit, node])
        end
        gen_dict[node]["PtH"][year] = value.(m.ext[:variables][:cap_PtH][node])
        gen_dict[node]["OCGT_H"][year] = value.(m.ext[:variables][:cap_OCGT][node])
    end
    return gen_dict
end

function get_generators_capacity_storage_relative(m::Model, I::Vector{String}, N::Vector{String}, S::Vector{String})
    gen_dict = Dict{String, Dict{String, Float64}}()  
    for node in N
        total_cap = sum(value.(m.ext[:variables][:cap][unit, node]) for unit in I) + value.(m.ext[:variables][:cap_PtH][node]) + value.(m.ext[:variables][:cap_OCGT][node])
        gen_dict[node] = Dict((gen => 0.0 for gen in cat(I, S, dims=1))...)
        for unit in I
            gen_dict[node][unit] = (value.(m.ext[:variables][:cap][unit, node]) / total_cap) * 100000
        end
        gen_dict[node]["PtH"] = (value.(m.ext[:variables][:cap_PtH][node]) / total_cap) * 100000
        gen_dict[node]["OCGT_H"] = (value.(m.ext[:variables][:cap_OCGT][node]) / total_cap) * 100000
    end
    return gen_dict
end

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

# Create a dictionaries to store the AC transmission capacity
function get_ac_transmission_capacity(m::Model, L_ac::Vector{Vector{String}})
    trans_ac_dict = Dict(la => value.(m.ext[:variables][:varlac][la]) for la in L_ac)
    return trans_ac_dict
end

function get_ac_transmission_capacity_fixed_network(m::Model, L_ac::Vector{Vector{String}})
    trans_ac_dict = Dict(la => value.(float(m.ext[:parameters][:Fl_MAX_AC][find_line_number(network["AC_Lines"], la)])) for la in L_ac)
    return trans_ac_dict
end


# Create a dictionaries to store the AC and DC transmission capacity
function get_dc_transmission_capacity(m::Model, L_dc::Vector{Vector{String}})
    trans_dc_dict = Dict(ld => value.(m.ext[:variables][:varldc][ld]) for ld in L_dc)
    return trans_dc_dict
end

function get_dc_transmission_capacity_fixed_network(m::Model, L_dc::Vector{Vector{String}})
    trans_dc_dict = Dict(ld => value.(float(m.ext[:parameters][:Fl_MAX_DC][find_line_number(network["DC_Lines"], ld)])) for ld in L_dc)
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
    return data
end

#Generator_colors = Dict("Biomass" => :green, "WindOffshore" => :blue, "CCGT_new" => :orange, "WindOnshore" => :lightblue, "Nuclear" => :purple, "Solar" => :gold, "OCGT" => :red, "ICE" => :gray)
Generator_colors = [:green, :blue, :orange, :lightblue, :purple, :gold, :red, :gray]    # [:green, :blue, :lightblue, :gold]  
Generator_labels = ["Biomass" "WindOffshore" "CCGT_new" "WindOnshore" "Nuclear" "Solar" "OCGT" "ICE"]   # ["Biomass" "WindOffshore" "WindOnshore" "Solar"]  

Generator_colors_storage = [:green, :blue, :orange, :magenta, :lightblue, :purple, :gold, :red, :gray, :cyan, :brown, :black, :lime, :teal]    # [:green, :blue, :lightblue, :gold]  
Generator_labels_storage = ["Biomass" "WindOffshore" "CCGT_new" "OCGT_H" "WindOnshore" "Nuclear" "Solar" "OCGT" "ICE" "PtH" "SPP_lignite" "SPP_coal" "Biofuel" "Hydro_RoR"]   # ["Biomass" "WindOffshore" "WindOnshore" "Solar"]
Generator_colors_dict = Dict(zip(Generator_labels_storage, Generator_colors_storage))

function plot_generator_capacities(data_matrix::Matrix{Float64}, nodes::Vector{String}, Generator_colors::Vector{Symbol}, Generator_labels::Matrix{String}, title::String)
    # Plot generation stacked bar chart
    groupedbar(data_matrix,
        bar_position = :stack,
        bar_width=0.5,
        bar_edges=true,
        ylabel= "Capacity (MW)",
        xticks=(1:length(nodes), nodes),
        label= Generator_labels,
        color_palette= Generator_colors,
        title= title,
        legend= :outertopright)
end


function plot_generator_capacities_new(gen_dict::Dict{String, Dict{String, Float64}}, Generator_colors::Dict{String, Symbol}, title::String)
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

    # Create the color palette and labels
    colors = [Generator_colors[unit_type] for unit_type in unit_types]
    labels = [unit_type for unit_type in unit_types]

    # Plot generation stacked bar chart
    groupedbar(data,
        bar_position = :stack,
        bar_width = 0.5,
        bar_edges = true,
        ylabel = "Capacity (MW)",
        xticks = (1:length(nodes), nodes),
        label = hcat(labels...),
        color_palette = colors,
        title = title,
        legend = :outertopright)
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



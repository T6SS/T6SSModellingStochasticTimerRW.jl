"""
Takes type `SolutionVarParDom` and returns the so-called theoretical distance.
"""
function get_t6ss_walk_theoretical_dist(sol::SolutionVarParDom)
    distance = 2*sqrt(
        (sol.parameters.λ₁*sol.parameters.h/sol.variables.Δx^2*sol.variables.T)/
        (sol.parameters.λ₀ + sol.parameters.λ₁))
        return distance        
end

"""
Biologically Inspired
Takes type `SolutionVarParDom` and returns the so-called theoretical distance.
"""
function get_t6ss_walk_theoretical_dist(::T6ssBiologicalInspired,sol::SolutionVarParDom)
    distance = 2*sqrt(
        log(2)*(sol.parameters.h/sol.parameters.λ₀))
        return distance        
end




        


# I can get the theoretical distance per group. I need to add the corresponding value 
# for each parameter and variable and then pair it with the outputted theoretical distance
# Once I have this, I need to make sure the appropriate label is available, but this could be 
# done through the plotting and saving functions.

function get_theoretical_distance(data,meta_data)
    return @chain data begin
        read_solution_from_memory(_,SolutionVarParDom)
        (get_value(_,meta_data[:data_type], meta_data[:data_value]),get_t6ss_walk_theoretical_dist(_))
    end
end

function get_theoretical_tuple(grouped_file_paths,iters_get_values)
    theoretical = []
    num_params = size(grouped_file_paths,1)
    for i in 1:num_params
        push!(theoretical,map(x -> get_theoretical_distance(x,iters_get_values[i]),grouped_file_paths[i]))
    end
    return theoretical
end







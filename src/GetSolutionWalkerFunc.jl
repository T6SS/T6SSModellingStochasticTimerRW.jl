###################################################
#### Solve the T6SS model and compute distance ####
###################################################
//
    # loads a .jl file of parameters
    # solves the T6SS equation
    # returns distances
    function solve(parameter_file::String)
        include(parameter_file)
        machine_positions = get_walker_position(λ₀Δt,λ₁Δt,hΔt,t⃗,N,walker,xwidth,ylength,boundary)
        return get_distance_count_sample_experimental(machine_positions,xwidth,ylength,N,boundary)
    end

    function solve(parameter_file::String;walker_ic::String = "random")
        include(parameter_file)
        if walker_ic == "random"
            machine_positions = get_walker_position(λ₀Δt,λ₁Δt,hΔt,t⃗,N,xwidth,ylength,boundary)
        else 
            machine_positions = get_walker_position(λ₀Δt,λ₁Δt,hΔt,t⃗,N,walker,xwidth,ylength,boundary)
        end
        return get_distance_count_sample_experimental(machine_positions,xwidth,ylength,N,boundary)
    end
 //
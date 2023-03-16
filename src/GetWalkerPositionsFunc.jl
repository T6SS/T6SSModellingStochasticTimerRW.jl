##############################
#### T6SS Model functions ####
##############################
//
    """
    Get the state time series
    """        
    function get_state_time_series(
        var::Variables,
        par::Parameters)
        state = s₀()
        states = []
        
        for i in t⃗(var)
            state = Telegraph.update(state,par.λ₀,par.λ₁)
            push!(states,state)
        end
        return states
    end

    

    """
    Run the T6SS for one cycle of State 0 and State 1
    Inspired for model that matches biological data
    """        
    function get_state_time_series(
        ::T6ssBiologicalInspired,
        var::Variables,
        par::Parameters)

        
        state = s₀()
        states = []
        
        cycle = 1
        while cycle < 3
            state = Telegraph.update(state,par.λ₀,par.λ₁)
        
            push!(states,state)
            sz_states = size(states,1)
            if cycle == 1 && (states[sz_states] == states[1])
                continue
            elseif cycle == 1 && (states[sz_states] != states[1])
                cycle = 2
            elseif cycle == 2 && (states[sz_states] == states[1])
                cycle  = 3
                states = states[1:sz_states-1]
            else
                # do nothing
            end
        end
        return states
    end
    





    """
    Get the complete random walk
    """        
    function get_all_walker_position(
        var::Variables,
        par::Parameters,
        dom::Domain)
        
        walker = rand_walker(dom,var)
        #w₀ = walker
        state = s₀()
        walkers = []
        for i in t⃗(var)
            push!(walkers,walker)
            state = Telegraph.update(state,par.λ₀,par.λ₁)
            walker = walker_update(walker,state,var,par,dom)
        end
        
        return RandomWalker.position.(walkers)
    end


    """
    Get the complete random walk and all of the states
    """        
    function get_all_walker_position(
        var::Variables,
        par::Parameters,
        dom::Domain)
        
        walker = rand_walker(dom,var)
        #w₀ = walker
        state = s₀()
        walkers = []
        states = []
        for i in t⃗(var)
            push!(walkers,walker)
            push!(states,state)
            state = Telegraph.update(state,par.λ₀,par.λ₁)
            walker = walker_update(walker,state,var,par,dom)
        end
        
        return T6ssOutput(states,walkers)
    end



    """
    Run the T6SS for one cycle of State 0 and State 1
    Inspired for model that matches biological data
    """        
    function get_all_walker_position(
        ::T6ssBiologicalInspired,
        var::Variables,
        par::Parameters,
        dom::Domain)

        walker = rand_walker(dom,var)
        
        state = s₀()
        states = []
        walkers = []
        cycle = 1
        while cycle < 3
            state = Telegraph.update(state,par.λ₀,par.λ₁)
            walker = walker_update(walker,state,var,par,dom)
            push!(walkers,walker)
            push!(states,state)
            sz_states = size(states,1)
            if cycle == 1 && (states[sz_states] == states[1])
                continue
            elseif cycle == 1 && (states[sz_states] != states[1])
                cycle = 2
            elseif cycle == 2 && (states[sz_states] == states[1])
                cycle  = 3
                walkers = walkers[1:sz_states-1]
                states = states[1:sz_states-1]
            else
                # do nothing
            end
        end
            
            return (states = states, walk = RandomWalker.position.(walkers))
        end

    """
    Get the first and last position of the random walk. Telegraph series is coupled to the random walk.
    This is a single iteration. The step is generated from the Multinomila distribution.
    """        
    function get_walker_position(
        var::Variables,
        par::Parameters,
        dom::Domain)
        
        walker = rand_walker(dom,var)
        w₀ = walker
        state = s₀()
        
        for i in t⃗(var)
            #state = Telegraph.update(state,par)
            #walker = walk_update(Random,walker,state,var,par,dom)
            state = Telegraph.update(state,par.λ₀,par.λ₁)
            walker = walker_update(walker,state,var,par,dom)
        end
        
        return WalkerFirstLast(w₀,walker)
    end


    """
        Generate the first and last position of the random walk. The step and telegraph series is all drawn from the same the random variable. 
    """
    function get_walker_position(
        var::Variables,
        par::Parameters,
        dom::Domain)
        
        # assume state = 0
        #  
        Step2D(0,1)


        walker = rand_walker(dom,var)
        w₀ = walker
        state = s₀()
        
        for i in t⃗(var)
            #state = Telegraph.update(state,par)
            #walker = walk_update(Random,walker,state,var,par,dom)
            state = Telegraph.update(state,par.λ₀,par.λ₁)
            walker = walker_update(walker,state,var,par,dom)
        end
        
        return WalkerFirstLast(w₀,walker)
    end

    """
    Run the T6SS for one cycle of State 0 and State 1
    Inspired for model that matches biological data
    """        
    function get_walker_position(
        ::T6ssBiologicalInspired,
        var::Variables,
        par::Parameters,
        dom::Domain)

        walker = rand_walker(dom,var)
        
        state = s₀()
        states = []
        walkers = []
        cycle = 1
        while cycle < 3
            state = Telegraph.update(state,par.λ₀,par.λ₁)
            walker = walker_update(walker,state,var,par,dom)
            push!(walkers,walker)
            push!(states,state)
            sz_states = size(states,1)
            if cycle == 1 && (states[sz_states] == states[1])
                continue
            elseif cycle == 1 && (states[sz_states] != states[1])
                cycle = 2
            elseif cycle == 2 && (states[sz_states] == states[1])
                cycle  = 3
                walkers = walkers[1:sz_states-1]
                states = states[1:sz_states-1]
            else
                # do nothing
            end
        end
            
            return WalkerFirstLast(walkers[1],walkers[end])
        end
        






    # Single iteration walker is predefined
    function get_walker_position(λ₀Δt,λ₁Δt,hΔt,t⃗,walker,xwidth,ylength,boundary)
        state = s₀()
        positions = []
        states = []
        T6SS = []
        push!(states = state)
        push!(positions,RandomWalker.position(walker))
        (machine_states = states, machine_positions = positions)
        if boundary == "no_boundary"
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = RandomWalker.update(walker, state,hΔt, xwidth, ylength)
                push!(states,state)
                push!(positions,RandomWalker.position(walker))
            end
        elseif boundary == "double_periodic"
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = RandomWalker.updateperiod(walker, state,hΔt, xwidth, ylength)
                push!(states,state)
                push!(positions,RandomWalker.position(walker))
            end
        elseif boundary == "no_flux"
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = RandomWalker.updatenoflux(walker, state,hΔt, xwidth, ylength)
                push!(states,state)
                push!(positions,RandomWalker.position(walker))
            end
        elseif boundary == "cell_boundary"
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = RandomWalker.updatecell(walker, state,hΔt, xwidth, ylength)
                push!(states,state)
                push!(positions,RandomWalker.position(walker))
            end
        else
            #nothing
        end
        return (machine_states = states, machine_positions = positions)
    end

    # single iteration, walker is chosen at random
    function get_walker_position(λ₀Δt,λ₁Δt,hΔt,t⃗,xwidth,ylength,boundary)
        walker = get_random_position(xwidth,ylength;step=1)
        state = s₀()
        positions = []
        states = []
        push!(positions,RandomWalker.position(walker))
        push!(states,state)
        (machine_states = states, machine_positions = positions)
        if boundary == "no_boundary"
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = RandomWalker.update(walker, state,hΔt, xwidth, ylength)
                push!(states,state)
                push!(positions,RandomWalker.position(walker))
            end
        elseif boundary == "double_periodic"
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = RandomWalker.updateperiod(walker, state,hΔt, xwidth, ylength)
                push!(states,state)
                push!(positions,RandomWalker.position(walker))
            end
        elseif boundary == "no_flux"
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = RandomWalker.updatenoflux(walker, state,hΔt, xwidth, ylength)
                push!(states,state)
                push!(positions,RandomWalker.position(walker))
            end
        elseif boundary == "cell_boundary"
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = RandomWalker.updatecell(walker, state,hΔt, xwidth, ylength)
                push!(states,state)
                push!(positions,RandomWalker.position(walker))
            end
        else
            #nothing
        end
        return (machine_states = states, machine_positions = positions)
    end

    # N iterations, walker predefined all single parameters
    function get_walker_position(
        λ₀Δt::Float64,
        λ₁Δt::Float64,
        hΔt::Float64,
        t⃗,
        N::Int64,
        walker,xwidth,ylength,boundary)
        state = s₀()
        last_position = []
        if boundary == "no_boundary"
            for j in 1:10^N
                for i in t⃗
                    state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                    walker = RandomWalker.update(walker, state,hΔt, xwidth, ylength)
                end
                push!(last_position,walker)
                walker = Walker2D(0,0)
            end
        elseif boundary == "double_periodic"
            for j in 1:10^N
                for i in t⃗
                    state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                    walker = RandomWalker.updateperiod(walker, state,hΔt, xwidth, ylength)
                end
                push!(last_position,walker)
                walker = Walker2D(0,0)
            end
        elseif boundary == "no_flux"
            for j in 1:10^N
                for i in t⃗
                    state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                    walker = RandomWalker.updatenoflux(walker, state,hΔt, xwidth, ylength)
                end
                push!(last_position,walker)
                walker = Walker2D(0,0)
            end
        elseif boundary == "cell_boundary"
            for j in 1:10^N
                for i in t⃗
                    state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                    walker = RandomWalker.updatecell(walker, state,hΔt, xwidth, ylength)
                end
                push!(last_position,walker)
                walker = Walker2D(0,0)
            end
        else
            #nothing
        end
        return (machine_positions = last_position)
    end

    # N iterations, walker sampled at random
    # single scalar value parameters
    function get_walker_position(
        λ₀Δt::Float64,
        λ₁Δt::Float64,
        hΔt::Float64,
        t⃗::StepRangeLen,
        N::Int64,
        xwidth::Int64,
        ylength::Int64,
        boundary::String)
        state = s₀()
        last_position = []
        if boundary == "no_boundary"
            for j in 1:10^N
                walker = get_random_position(xwidth,ylength;step=1)
                w₀ = walker
                for i in t⃗
                    state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                    walker = RandomWalker.update(walker, state,hΔt, xwidth, ylength)
                end
                push!(last_position,(first_walker = w₀,last_walker = walker))
            end
        elseif boundary == "double_periodic"
            for j in 1:10^N
                walker = get_random_position(xwidth,ylength;step=1)
                w₀ = walker
                for i in t⃗
                    state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                    walker = RandomWalker.updateperiod(walker, state,hΔt, xwidth, ylength)
                end
                push!(last_position,(first_walker = w₀,last_walker = walker))
            end
        elseif boundary == "no_flux"
            for j in 1:10^N
                walker = get_random_position(xwidth,ylength;step=1)
                w₀ = walker
                for i in t⃗
                    state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                    walker = RandomWalker.updatenoflux(walker, state,hΔt, xwidth, ylength)
                end
                push!(last_position,(first_walker = w₀,last_walker = walker))
            end
        elseif boundary == "cell_boundary"
            for j in 1:10^N
                walker = get_random_position(xwidth,ylength;step=1)
                w₀ = walker
                for i in t⃗
                    state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                    walker = RandomWalker.updatecell(walker, state,hΔt, xwidth, ylength)
                end
                push!(last_position,(first_walker = w₀,last_walker = walker))
            end
        else
            #nothing
        end
        return last_position
    end


        # N iterations, walker sampled at random
    # single scalar value parameters



    # Simulated results
    function simulate_random_walk(xwidth,ylength,N)
        w() = (first_walker = get_random_position(xwidth,ylength;step=1),
        last_walker = get_random_position(xwidth,ylength;step=1))
        sols = [w() for i in 1:10^N]
        return (machine_positions = sols)
    end

    # Updated function
    # Single parameters
    function get_walker_position(
        λ₀Δt::Float64,
        λ₁Δt::Float64,
        hΔt::Float64,
        t⃗::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
        N::Int64,
        xwidth::Int64,
        ylength::Int64,
        boundary::String,
        boundary_conditional::Function)
        
        next_position = boundary_conditional(
        boundary,(
            RandomWalker.update,
            RandomWalker.updateperiod,
            RandomWalker.updatenoflux,
            RandomWalker.updatecell))
    
        last_position = []
        for j in 1:10^N
            walker = get_random_position(xwidth,ylength;step=1)
            w₀ = walker
            state = s₀()
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = next_position(walker, state,hΔt, xwidth, ylength)
            end
            push!(last_position,(first_walker = w₀,last_walker = walker))
        end
        
        return last_position
    end

    # Updated function
    # Iterate over hΔt
    function get_walker_position(
        λ₀Δt::Float64,
        λ₁Δt::Float64,
        hΔt::StepRangeLen,
        t⃗::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
        N::Int64,
        xwidth::Int64,
        ylength::Int64,
        boundary::String,
        boundary_conditional::Function)
        
        next_position = boundary_conditional(
        boundary,(
            RandomWalker.update,
            RandomWalker.updateperiod,
            RandomWalker.updatenoflux,
            RandomWalker.updatecell))
    
        last_position = []
        for j in 1:10^N
            walker = get_random_position(xwidth,ylength;step=1)
            w₀ = walker
            state = s₀()
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = next_position(walker, state,hΔt, xwidth, ylength)
            end
            push!(last_position,(first_walker = w₀,last_walker = walker))
        end
        
        return last_position
    end

    # Updated function
    # Iterate over λ₀Δt
    function get_walker_position(
        λ₀Δt::StepRangeLen,
        λ₁Δt::Float64,
        hΔt::Float64,
        t⃗::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
        N::Int64,
        xwidth::Int64,
        ylength::Int64,
        boundary::String,
        boundary_conditional::Function)
        
        next_position = boundary_conditional(
        boundary,(
            RandomWalker.update,
            RandomWalker.updateperiod,
            RandomWalker.updatenoflux,
            RandomWalker.updatecell))
    
        last_position = []
        for j in 1:10^N
            walker = get_random_position(xwidth,ylength;step=1)
            w₀ = walker
            state = s₀()
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = next_position(walker, state,hΔt, xwidth, ylength)
            end
            push!(last_position,(first_walker = w₀,last_walker = walker))
        end
        
        return last_position
    end

    # Updated function
    # Single parameters
    function get_walker_position(
        λ₀Δt::Float64,
        λ₁Δt::Float64,
        hΔt::Float64,
        t⃗::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
        N::Int64,
        xwidth::Int64,
        ylength::Vector,
        boundary::String,
        boundary_conditional::Function)
        
        next_position = boundary_conditional(
        boundary,(
            RandomWalker.update,
            RandomWalker.updateperiod,
            RandomWalker.updatenoflux,
            RandomWalker.updatecell))
    
        last_position = []
        for j in 1:10^N
            walker = get_random_position(xwidth,ylength;step=1)
            w₀ = walker
            state = s₀()
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = next_position(walker, state,hΔt, xwidth, ylength)
            end
            push!(last_position,(first_walker = w₀,last_walker = walker))
        end
        
        return last_position
    end


    function get_walker_position(
        λ₀Δt::Float64,
        λ₁Δt::Float64,
        hΔt::Float64,
        t⃗::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
        N::Int64,
        xwidth::Int64,
        ylength::Vector,
        boundary::String,
        boundary_conditional::Function)
        
        next_position = boundary_conditional(
        boundary,(
            RandomWalker.update,
            RandomWalker.updateperiod,
            RandomWalker.updatenoflux,
            RandomWalker.updatecell))
    
        last_position = []
        for j in 1:10^N
            walker = get_random_position(xwidth,ylength;step=1)
            w₀ = walker
            state = s₀()
            for i in t⃗
                state = Telegraph.update(state,λ₀Δt,λ₁Δt)
                walker = next_position(walker, state,hΔt, xwidth, ylength)
            end
            push!(last_position,(first_walker = w₀,last_walker = walker))
        end
        
        return last_position
    end

//
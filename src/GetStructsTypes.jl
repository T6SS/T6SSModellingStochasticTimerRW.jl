#################
#### Structs ####
#################
//
    struct SolutionStructure
        machine_positions::Vector{Tuple{Int64, Int64}}
        indx::Vector{Int64}
        distances_experimental::Vector{Int64}
        count_experimental::Vector{Int64}
        distance_sample::Vector{Int64}
        count_sample::Vector{Int64}
        count_norm::Vector{Float64}
    end

    struct DistanceExperSampleStruct
        experimental::Vector{Union{Float64,Int64}}
        sample::Vector{Union{Float64,Int64}}
    end

    struct DistanceStructure
        parameter_range::Union{StepRangeLen,Vector}
        mean_distances::Vector{Float64}
    end
    
    

    abstract type T6ssFunctionInputs end
    abstract type T6ssFunctionOutput end

    struct T6ssOutput <: T6ssFunctionOutput
        states
        walk
    end


    struct T6ssBiologicalInspired <: T6ssFunctionInputs end

    mutable struct Variables <: T6ssFunctionInputs
        Δt::Float64
        Δx::Union{Int64,Float64}
        T::Int64
        X::Union{Int64,Float64}
        aspect::Union{Int64,Float64}
        boundary::String
    end

    mutable struct  Iterations    
        simulation_iterations::Int64
        parameter_iterations::Int64
    end
    
    struct Parameters <: T6ssFunctionInputs
        λ₀::Union{Float64,Int64}
        λ₁::Union{Float64,Int64}
        h::Union{Float64,Int64}
    end

    struct Domain
        width::Union{Float64,Int64}
        length::Union{Float64,Int64}
    end

    struct WalkerFirstLast
        w₁::Walker2D
        w₂::Walker2D
    end

    struct T6ssDistance
        distance::Union{Int64,Float64}
    end
    
    
    struct Distance
        distance_x::Union{Int64,Float64}
        distance_y::Union{Int64,Float64}
    end

    struct SolutionVarParDom 
        variables::Variables
        parameters::Parameters
        domain::Domain
        iterations::Iterations
        experimental::Vector{Union{Float64,Int64}}
        sample::Vector{Union{Float64,Int64}}
    end
    
    struct T6ssDataGenerateInput
        Δt::Float64
        Δx::Union{Int64,Float64}
        T::Int64
        X::Int64
        aspect::Union{Int64,Float64}
        boundary::String
        simulation_iterations::Int64
        parameter_iterations::Int64
        λ₀::Union{Float64,Int64}
        λ₁::Union{Float64,Int64}
        h::Union{Float64,Int64}
        dif_min::Float64
        dif_max::Float64
        λmax::Int64
        hmax::Int64
        amax::Int64
        data_root::String
    end
//

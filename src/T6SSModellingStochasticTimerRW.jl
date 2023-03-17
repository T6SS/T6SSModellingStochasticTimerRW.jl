module T6SSModellingStochasticTimerRW

using Reexport

@reexport using Revise
@reexport using Profile
@reexport using RandomWalker
@reexport using Telegraph
@reexport using AlgebraOfGraphics
@reexport using JSON3
@reexport using Distances
@reexport using Colors
@reexport using StructTypes
@reexport using Chain
@reexport using CairoMakie
@reexport using Distributed
@reexport using Statistics
@reexport using StatsBase
@reexport using ThreadsX
@reexport using UnPack
@reexport using Setfield
@reexport using NaturalSort


export
    SolutionStructure,
    DistanceExperSampleStruct,
    DistanceStructure,
    T6ssFunctionInputs,
    T6ssFunctionOutput,
    T6ssOutput,
    T6ssBiologicalInspired,
    Variables,
    Iterations,
    Domain,
    WalkerFirstLast,
    T6ssDistance,
    Distance,
    SolutionVarParDom,
    T6ssDataGenerateInput,
    generate_t6ss_data,
    get_t6ss_walk_theoretical_dist,
    get_theoretical_distance,
    get_theoretical_tuple,
    get_sample_distribution,
    get_periodic_boundary_distance,
    get_no_periodic_boundary_distance,
    get_distance,
    get_walker_distance,
    get_experimental_sample_dist_vecs,
    set_range_iter,
    get_count,
    find,
    write_solution_to_file,
    read_solution_from_memory,
    join_figure_number_letter,
    join_name_letters,
    append_forward_slash_figure,
    join_figure_number_tuple_letters,
    data_path,
    even,
    t⃗,
    x⃗,
    dimensions,
    Δ,
    s₀,
    LitRate,
    param_depend_distplacement,
    get_value,
    get_mean_distance,
    get_mean_vector,
    get_solutions_vector,
    solve,
    rand_walker,
    boundary_conditional,
    get_boundary_distance_functions,
    boundary_conditions,
    Telegraph.update,
    T6ss_next_choice,
    get_state_time_series,
    get_all_walker_position,
    get_walker_position,
    simulate_random_walk,
    run_simulation,
    solve_iterations,
    visualise_save_solution,
    view_distance_and_mean,
    view_mean_heatmap,
    view_scatter_mean,
    save_figure


include("GetFileAndNamingFuncjl")
include("GenerateT6ssData.jl")
include("GenerateTheoreticalData.jl")
include("GetDistanceFunc.jl")
include("GetHelperFuncjl")
include("GetSolutionWalkerFunc.jl")
include("GetStructsTypes.jl")
include("GetWalkerHelperFunc.jl")
include("GetWalkerPositionsFunc.jl")
include("ViewMultipleParamsIters.jl")
include("ViewWalkerDistSolutionFunc.jl")
include("GetHistogramDistancesT6ss.jl")




end

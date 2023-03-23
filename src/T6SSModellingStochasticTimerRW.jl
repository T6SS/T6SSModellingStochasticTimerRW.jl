module T6SSModellingStochasticTimerRW

using Reexport

@reexport using Revise
@reexport using Base
@reexport using Profile
@reexport using Telegraph
@reexport using RandomWalker
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
    generate_figure_3,
    generate_figure_4,
    generate_figure_5

include("GetStructsTypes.jl")
include("GetFileAndNamingFunc.jl")
include("GenerateT6ssData.jl")
include("GenerateTheoreticalData.jl")
include("GetDistanceFunc.jl")
include("GetHelperFunc.jl")
include("GetSolutionWalkerFunc.jl")
include("GetWalkerHelperFunc.jl")
include("GetWalkerPositionsFunc.jl")
include("ViewMultipleParamsIters.jl")
include("ViewWalkerDistSolutionFunc.jl")
include("GetHistogramDistancesT6ss.jl")
include("GenerateFigures.jl")



end

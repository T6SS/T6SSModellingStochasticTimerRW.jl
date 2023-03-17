# T6SSModellingStochasticTimerRW

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fieldofnodes.github.io/ModellingStochasticTimerRW.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fieldofnodes.github.io/ModellingStochasticTimerRW.jl/dev/)
[![Build Status](https://github.com/fieldofnodes/ModellingStochasticTimerRW.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fieldofnodes/ModellingStochasticTimerRW.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/fieldofnodes/ModellingStochasticTimerRW.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/fieldofnodes/ModellingStochasticTimerRW.jl)


## Installation
This package depends on two other unregistered packages which I wrote, `Telegraph.jl` and `RandomWalker.jl`. To install this package, `T6SSModellingStochasticTimerRW.jl` open julia in your preferred way.

### `REPL`
```julia
] add https://github.com/fieldofnodes/Telegraph.jl
] add https://github.com/fieldofnodes/RandomWalker.jl
] add https://github.com/fieldofnodes/T6SSModellingStochasticTimerRW.jl
```
### Script
```julia
using Pkg
Pkg.add(url = "https://github.com/fieldofnodes/Telegraph.jl")
Pkg.add(url = "https://github.com/fieldofnodes/RandomWalker.jl")
Pkg.add(url = "https://github.com/fieldofnodes/T6SSModellingStochasticTimerRW.jl")
```

Once complete run

```julia
using Telegraph
using RandomWalker
using T6SSModellingStochasticTimerRW
```


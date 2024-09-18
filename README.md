# GeneralizedBM.jl

GeneralizedBM uses the tight-binding model and the [continnum model](
https://doi.org/10.48550/arXiv.2406.15712
) to simulate the band structure of two models of twisted bilayer graphen: a [Wannierized physics model](https://doi.org/10.1103/PhysRevB.98.075106), and a [simplified TBG model](https://doi.org/10.1063/5.0115771).

## Installation
GeneralizedBM.jl is an unregistered package and therefore needs to be downloaded or cloned to the user's local computer first, and then installed by running

```julia
julia> cd("your-local-path/GeneralizedBM.jl")
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```

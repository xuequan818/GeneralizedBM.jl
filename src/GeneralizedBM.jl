module GeneralizedBM

using Arpack, LinearAlgebra, KrylovKit
using Printf, Plots, Plots.PlotMeasures, LaTeXStrings
using StaticArrays, SparseArrays
using ForwardDiff

include("TBL.jl")
include("basis.jl")
include("hopping.jl")
include("hamiltonian.jl")

end # module GeneralizedBM

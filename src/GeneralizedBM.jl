module GeneralizedBM

using Arpack, LinearAlgebra, KrylovKit
using Printf, Plots, Plots.PlotMeasures, LaTeXStrings
using StaticArrays, SparseArrays
using Roots, TaylorDiff, ForwardDiff
using FastGaussQuadrature, Dierckx

include("TBL.jl")
include("basis.jl")
include("hopping_struct.jl")
include("hopping_toy.jl")
include("hopping_tbg.jl")
include("hopping_BM.jl")
include("Hamiltonian.jl")
include("band.jl")

end # module GeneralizedBM

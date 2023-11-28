module GeneralizedBM

using Arpack, LinearAlgebra, KrylovKit
using Printf, Plots, Plots.PlotMeasures, LaTeXStrings
using StaticArrays, SparseArrays
using Roots, TaylorDiff
using FastGaussQuadrature, Dierckx

include("TBL.jl")
include("basis.jl")
include("hopping.jl")
include("hopping_tbg.jl")
include("HamiltonianGBM.jl")
include("HamiltonianMS.jl")
include("band.jl")

end # module GeneralizedBM

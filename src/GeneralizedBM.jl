module GeneralizedBM

using Arpack, LinearAlgebra, KrylovKit
using Printf, Plots, Plots.PlotMeasures, LaTeXStrings
using StaticArrays, SparseArrays
using ForwardDiff, Roots

include("TBL.jl")
include("basis.jl")
include("hopping.jl")
include("HamiltonianBM.jl")
include("HamiltonianMS.jl")

end # module GeneralizedBM

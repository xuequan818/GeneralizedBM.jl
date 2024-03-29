using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
rcut = 0.3 # cutoff of the basis

# define the TBL model
Lat = TBLG(θ; a=2.46)
tau = 3
hop = hopToy(Lat; τ = tau)
basis = Basis(rcut, Lat);
pind, P = band_plot(Lat, basis, hop, 0.005; num=15)
plot!(P, title=L"\theta = %$θ^\circ, \tau = %$tau")

using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
# define the TBL model
Lat = TBLG(θ; a=2.46)
tau = 1
hop = hopTBG(Lat;τinter=tau)

rcut = 30.0 # cutoff of the basis
basis = Basis(rcut, Lat);
P = band_plot(Lat, basis, hop, 0.02; num=10)
plot!(P, title=L"\theta = %$θ^\circ, \tau = %$tau")

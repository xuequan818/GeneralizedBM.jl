using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 

# define the TBG model
Lat = TBLG(θ; a=2.46)
@time hop = hopTBG(Lat);

rcut = 0.2 # cutoff of the basis
basis = Basis(rcut, Lat);
pind, P = band_plot(Lat, basis, hop, 0.02; nE = 8, num=15)
plot!(P, ylims=(-0.26, 0.3), title="Wannierized")

savefig("pics/ms_tbg.pdf")
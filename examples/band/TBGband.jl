using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 

# define the TBG model
Lat = TBLG(θ; a=2.46)
@time hop = hopTBG(Lat);

rcut = 20.0 # cutoff of the basis
basis = Basis(rcut, Lat);
pind, P = band_plot(Lat, basis, hop, 0.02; num=15)
plot!(P, title=L"%$θ^\circ")

savefig("pics/ms_tbg.pdf")
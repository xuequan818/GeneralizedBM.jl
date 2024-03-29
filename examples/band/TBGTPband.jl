using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
rcut = 0.2 # cutoff of the basis
p1 = 2 # intralayer expansion order
p2 = 1 # interlayer expansion order
tau = 2 # interlayer hopping truncation

# define the TBL model
Lat = TBLG(θ;a=2.46);
@time hop = hopTBG(Lat; Pintra=p1, Pinter=p2, τinter = tau)
basis = Basis(rcut, Lat);

pind, P = band_plot(Lat, basis, hop, 0.01;nE=8, num=15)
plot!(P, ylims=(-0.26, 0.3), title=L"(%$p1,%$p2,%$tau)")
savefig("pics/212_tbg.pdf")



using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
rcut = 20. # cutoff of the basis
p1 = 2 # intralayer expansion order
p2 = 1 # interlayer expansion order
tau = 2 # interlayer hopping truncation

# define the TBL model
Lat = TBLG(θ;a=2.46);
@time hop = hopTBG(Lat; Pintra=p1, Pinter=p2, τinter = tau)
basis = Basis(rcut, Lat);

pind, P = band_plot(Lat, basis, hop, 0.01; num=15)
plot!(P, title=L"%$θ^\circ,\,\, (%$p1,%$p2,%$tau)")
lens!([pind[2] - 0.15, pind[2] + 0.15], [-0.015, 0.025], inset=(1, bbox(0.3, 0.15, 0.2, 0.2)), subplot=2, ticks=nothing, box=:on)



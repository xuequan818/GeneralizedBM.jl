using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
nE = 3
rcut = 30.

# define the TBL model
Lat = TBLG(θ);
basis = Basis(rcut, Lat);

# generate exact band at Gamma point
q = [Lat.KM[1][1] + norm(Lat.KM[1] - Lat.KM[2]) * sqrt(3) / 2, 0.0]
hop = hopGBM(Lat)
@time H0 = ham_MS(Lat, basis, hop, q)
E0 = band(H0, nE)

e = []
tau = collect(1:6)
for t in tau
    hop = hopGBM(Lat; τ = t)
    @time H = ham_MST(Lat, basis, hop, q)
    E = band(H, nE)
    push!(e, norm(E - E0, Inf))
end
P = plot(tau, e, yscale=:log10, ylabel="Error", xlabel=L"\tau",guidefontsize=22, color=:black, title=L"\theta = %$θ^\circ", label="", tickfontsize=20, legendfontsize=20, legend=:topright, grid=:off, box=:on, size=(740, 620), titlefontsize=30, right_margin=3mm, top_margin=3mm, lw=2, marker=:circle, markersize=8, markercolor=:white,markerstrokecolor=:black)

using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
nE = 20

# define the TBG model
Lat = TBLG(θ; a=2.46)
@time hop = hopTBG(Lat)

# generate exact band at Gamma point
q = [Lat.KM[1][1] + norm(Lat.KM[1] - Lat.KM[2]) * sqrt(3) / 2, 0.0]
basis = Basis(norm(Lat.KM[1]) *0.8, Lat);
@time H0 = hamiltonian(Lat, basis, hop, q)
fv = 0.02
E0 = band(H0, nE; fv=fv)
E0 = E0[1][E0[2]]

e = []
rcut = collect(0.12:0.06:0.42)
for r in rcut
    @time H = hamiltonian(Lat, Basis(r, Lat), hop, q)
    E = band(H, nE; fv=fv)
    E = E[1][E[2]]
    push!(e, norm(E - E0, 2))
end
P = plot(rcut, e, yscale=:log10, ylabel="Error", xlabel="r",guidefontsize=22, color=:black, title=L"\theta = %$θ^\circ", label="", tickfontsize=20, legendfontsize=20, legend=:topright, grid=:off, box=:on, size=(740, 620), titlefontsize=30, right_margin=3mm, top_margin=3mm, lw=2, marker=:circle, markersize=8, markercolor=:white,markerstrokecolor=:black)

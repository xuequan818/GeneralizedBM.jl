using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
nE = 3

# define the TBL model
Lat = TBLG(θ;a=2.46);
hop = hopToy(Lat)

# generate exact band at Gamma point
q = [Lat.KM[1][1] + norm(Lat.KM[1] - Lat.KM[2]) * sqrt(3) / 2, 0.0]
basis = Basis(norm(Lat.KM[1])*2/3, Lat);
@time H0 = hamiltonian(Lat, basis, hop, q)
fv = 0.0053
E0 = band(H0, nE; fv=fv)

e = []
rcut = collect(0.1:0.05:0.4)
for r in rcut
    @time H = hamiltonian(Lat, Basis(r, Lat), hop, q)
    E = band(H, nE; fv = fv)
    push!(e, norm(E - E0, Inf))
end
P = plot(rcut, e, yscale=:log10, ylabel="Error", xlabel="r",guidefontsize=22, color=:black, title=L"\theta = %$θ^\circ", label="", tickfontsize=20, legendfontsize=20, legend=:topright, grid=:off, box=:on, size=(740, 600), titlefontsize=30, left_margin=2mm,right_margin=2mm, top_margin=3mm, lw=2, marker=:circle, markersize=8, markercolor=:white,markerstrokecolor=:black)

using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
nE = 3

# define the TBG model
Lat = TBLG(θ; a=2.46)
basis = Basis(20., Lat);

# generate exact band at Gamma point
q = [Lat.KM[1][1] + norm(Lat.KM[1] - Lat.KM[2]) * sqrt(3) / 2, 0.0]
@time hop = hopTBG(Lat; model="MS");
@time H0 = ham_MS(Lat, basis, hop, q)
fv = 0.02
E0 = band(H0, nE; fv=fv)

e = []
tau = collect(1:6)
for t in tau
    @time hop = hopTBG(Lat; model="MST", τinter=t)
    @time H = ham_MST(Lat, basis, hop, q)
    E = band(H, nE; fv=fv)
    push!(e, norm(E - E0, Inf))
end
P = plot(tau, e, yscale=:log10, ylabel="Error", xlabel=L"\tau",guidefontsize=22, color=:black, title=L"\theta = %$θ^\circ", label="", tickfontsize=20, legendfontsize=20, legend=:topright, grid=:off, box=:on, size=(740, 620), titlefontsize=30, right_margin=3mm, top_margin=3mm, lw=2, marker=:circle, markersize=8, markercolor=:white,markerstrokecolor=:black)

savefig("hoptrunc_tbg.pdf")

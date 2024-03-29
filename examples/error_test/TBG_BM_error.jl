using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
rcut = 0.2 # cutoff of the basis
tau = 2
nE = 3

# define the TBG model
Lat = TBLG(θ;a=2.46);
basis = Basis(rcut, Lat);
@time hop = hopTBG(Lat;τinter=tau)

fv = 0.02
@time band_info = band_path(Lat, basis, hop, fv; num=10, nE=nE, tol=0.1)

M = collect(0:3)
e1 = []
for m in M
    println("intra order = $(m)")
    fv = m > 0 ? 0.02 : 0.001
	@time hop = hopTBG(Lat; Pintra = m, τinter=tau)
    band_info_test = band_path(Lat, basis, hop, fv; num=10, nE=nE, tol=0.1)
    push!(e1, norm(band_info[1] - band_info_test[1], Inf))
end
P = plot(M, e1, yscale=:log10, ylabel="Error", guidefontsize=22, color=:black, title="", label=L"\mathfrak{m}", tickfontsize=20, legendfontsize=20, legend=:topright, grid=:off, box=:on, size=(740, 600), titlefontsize=30, left_margin=2mm, right_margin=2mm, top_margin=3mm, lw=2, marker=:circle, markersize=8, markercolor=:white, markerstrokecolor=:black)

e2 = []
for m in M
    println("inter order = $(m)")
	fv = m > 0 ? 0.02 : 0.002
	@time hop = hopTBG(Lat; Pinter=m, τinter=tau)
    band_info_test = band_path(Lat, basis, hop, fv; num=10, nE=nE, tol=0.1)
    push!(e2, norm(band_info[1] - band_info_test[1], Inf))
end
plot!(M, e2, xlabel="Order", ylabel="Error", label=L"\mathfrak{n}", color=:red, lw=2, marker=:circle, markersize=8, markercolor=:white, markerstrokecolor=:red)

savefig("pics/expansion_tbg_2t_path.pdf")

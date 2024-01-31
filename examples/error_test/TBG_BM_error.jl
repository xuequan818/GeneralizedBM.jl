using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
rcut = 20. # cutoff of the basis
tau = 2
nE = 3

# define the TBG model
Lat = TBLG(θ;a=2.46);
basis = Basis(rcut, Lat);
@time hop = hopTBG(Lat;τinter=tau)

# generate band at momentum q
q = [Lat.KM[1][1] + norm(Lat.KM[1] - Lat.KM[2]) * sqrt(3) / 2, 0.0]
@time H0 = hamiltonian(Lat, basis, hop, q)
fv = 0.02
E0 = band(H0, nE;fv=fv)

M = collect(0:3)
e1 = []
for m in M
    println("intra order = $(m)")
    fv = m > 0 ? 0.02 : 0.002
	@time hop = hopTBG(Lat; Pinter=5, Pintra = m, τinter=tau)
    @time H = hamiltonian(Lat, basis, hop, q)
    E = band(H, nE;fv=fv)
	push!(e1, norm(E - E0, Inf))
end
P = plot(M, e1, yscale=:log10, xticks=[0, 3, 6], ylabel="Error", guidefontsize=22, color=:black, title=L"\theta = %$θ^\circ", label=L"\mathfrak{m}", tickfontsize=20, legendfontsize=20, legend=:topright, grid=:off, box=:on, size=(740, 620), titlefontsize=30, right_margin=3mm, top_margin=3mm, lw=2, marker=:circle, markersize=8, markercolor=:white,markerstrokecolor=:black)

e2 = []
for m in M
    println("inter order = $(m)")
	fv = m > 0 ? 0.02 : 0.002
	@time hop = hopTBG(Lat; Pinter=m, Pintra = 5, τinter=tau)
    @time H = hamiltonian(Lat, basis, hop, q)
    E = band(H, nE;fv=fv)
	push!(e2, norm(E - E0, Inf))
end
plot!(M, e2, xlabel="Order", ylabel="Error", label=L"\mathfrak{n}", color=:red, lw=2, marker=:circle, markersize=8, markercolor=:white, markerstrokecolor=:red)

savefig("expansion_tbg_2t.pdf")

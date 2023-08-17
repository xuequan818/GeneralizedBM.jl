using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
rcut = 30. # cutoff of the basis
tau = 4
nE = 3

# define the TBL model
Lat = TBLG(θ);
basis = Basis(rcut, Lat);
hop = hopGBM(Lat;τ=tau)

# generate band at momentum q
q = [Lat.KM[1][1] + norm(Lat.KM[1] - Lat.KM[2]) * sqrt(3) / 2, 0.0]
@time H0 = ham_MST(Lat, basis, hop, q)
E0 = band(H0, nE)

M = collect(0:8)
e1 = []
for m in M
	hop = hopGBM(Lat; Pinter=10, Pintra = m, τ=tau)
    @time H = ham_GBM(Lat, basis, hop, q)
    E = band(H, nE)
	push!(e1, norm(E - E0, Inf))
end
P = plot(M, e1, yscale=:log10, xticks=[0, 4, 8], ylabel="Error", guidefontsize=22, color=:black, title=L"\theta = %$θ^\circ", label=L"P_{intra}", tickfontsize=20, legendfontsize=20, legend=:topright, grid=:off, box=:on, size=(740, 620), titlefontsize=30, right_margin=3mm, top_margin=3mm, lw=2, marker=:circle, markersize=8, markercolor=:white,markerstrokecolor=:black)

e2 = []
for m in M
	fv = m > 0 ? 0.01 : 0.
	hop = hopGBM(Lat; Pinter=m, Pintra = 10, τ=tau)
    @time H = ham_GBM(Lat, basis, hop, q)
    E = band(H, nE;fv=fv)
	push!(e2, norm(E - E0, Inf))
end
plot!(M, e2, ylabel="Error", label = L"P_{inter}", color=:red, lw=2, marker=:circle, markersize=8, markercolor=:white,markerstrokecolor=:red)

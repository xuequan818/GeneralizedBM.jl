using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
nE = 3
rcut = 0.3

# define the TBL model
Lat = TBLG(θ;a=2.46);
basis = Basis(rcut, Lat);

hop = hopToy(Lat)
fv = 0.0053
@time band_info = band_path(Lat, basis, hop, fv; num=10, nE=nE, tol=0.1)
#P = band_plot(band_info)[2]

e = []
tau = collect(1:6)
for t in tau
    hop = hopToy(Lat; τ = t)
    band_info_test = band_path(Lat, basis, hop, fv; num=10, nE=nE, tol=0.1)
    push!(e, norm(band_info[1] - band_info_test[1], Inf))
end
P = plot(tau, e, yscale=:log10, ylabel="Error", xlabel=L"\tau",guidefontsize=22, color=:black, title="", label="", tickfontsize=20, legendfontsize=20, legend=:topright, grid=:off, box=:on, size=(740, 600), titlefontsize=30, left_margin=2mm,right_margin=2mm, top_margin=3mm, lw=2, marker=:circle, markersize=8, markercolor=:white,markerstrokecolor=:black)

savefig("pics/hoptrunc_path.pdf")

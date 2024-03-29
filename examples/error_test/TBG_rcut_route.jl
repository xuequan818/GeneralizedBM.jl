using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra
using DelimitedFiles


θ = 1.1 # twist angle 
# define the TBG model
Lat = TBLG(θ; a=2.46)
@time hop = hopTBG(Lat)

# generate band structure for symmetric path
basis = Basis(norm(Lat.KM[1]) * 0.8, Lat);
fv = 0.02
nE = 32
#band_info = band_path(Lat, basis, hop, fv; num=10, nE=nE, tol = 0.1)
#band_info = band_info[1]
#P = band_plot(band_info)[2]
band_info = readdlm("bs.txt", '\t');

function band_cut(Eq, Emax)
    lb = findfirst(x -> x > -Emax, Eq)
    up = findfirst(x -> x > Emax, Eq)

    lb, up-1
end

rcut = collect(0.12:0.001:0.42)
band_info_test = [@time band_path(Lat, Basis(r, Lat), hop, fv; num=10, nE=nE) for r in rcut]
#band_plot(band_info_test[4])[2]

Erange = collect(0.03:0.001:0.65)
e = zeros(length(Erange), length(rcut))
for i = 1 : length(Erange), j = 1 : length(rcut)
    bs_cut = map(x -> band_cut(band_info[:, x], Erange[i]), 1:size(band_info, 2))

    #bs_cut = map(x -> band_cut(band_info[:, x], 0.42), 1:size(band_info, 2))

    band_test = band_info_test[j][1]
    er_vec = [norm(band_info[:, l][bs_cut[l][1]:bs_cut[l][2]] - band_test[:, l][bs_cut[l][1]:bs_cut[l][2]],Inf) for l = 1:length(bs_cut)]

    e[i,j] = norm(er_vec,Inf)
end
#P = plot(rcut, e[1,:], yscale=:log10, ylabel="Error", xlabel="r", guidefontsize=22, color=:black, title=L"\theta = %$θ^\circ", label="", tickfontsize=20, legendfontsize=20, legend=:topright, grid=:off, box=:on, size=(740, 620), titlefontsize=30, right_margin=3mm, top_margin=3mm, lw=2, marker=:circle, markersize=8, markercolor=:white, markerstrokecolor=:black)
#plot!(rcut, e[20 ,:])
#heatmap(rcut, Erange, log10.(e), grid=:on, size=(750, 600), tickfontsize=13, legendfontsize=18, guidefontsize=22, titlefontsize=32, zticks=nothing, border=:none, right_margin=7mm, left_margin=3mm, cbar=:best)
#plot!(camera=(0, 90))
using Interpolations
rcut2 = collect(rcut[1]:0.0002:rcut[end])
e2 = zeros(size(Erange,1),size(rcut2,1))
for i = 1:size(e,1)
    spl = linear_interpolation(rcut, log10.(e[i,:]))
    e2[i,:] = spl(rcut2)
end
P1 = heatmap(rcut2, Erange, e2, grid = :off, size=(1000, 700), left_margin=3mm,right_margin=6mm, bottom_margin=2mm,colorbar=:true, fill=:true, box=:on, tickfontsize=20, legendfontsize=18, guidefontsize=26, titlefontsize=30, colorbar_tickfontsize=16, xlabel=L"\Lambda", ylabel=L"\Sigma")
for yi in range(extrema(Erange)..., 20)
    annotate!(maximum(rcut2) + 0.035, yi, text("█", :white, 45))
end
for (i, yi) in enumerate(-0.029 .+0.0818 .* collect(8: -1:1))
    annotate!(maximum(rcut2) + 0.039 , yi, text(L"10^{-%$i}", :black, 22))
end
P1

savefig("pics/rtest_tbg.pdf")

writedlm("error.txt", e, '\t');
using GeneralizedBM
using Printf, Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
Lat = TBLG(θ; a=2.46)

hftrl, hftim = hopTBG_inter_strength(Lat)

L = 5.
xx = collect(-L:0.05:L)
yy = xx
r = [[x, y] for x in xx for y in yy]
vtot = [reshape(sqrt.(hftrl[i].(r).^2 .+ hftim[i].(r).^2), length(xx), length(yy)) for i = 1:4]
P = heatmap(xx, yy, vtot[1], size=(700, 600), c=cgrad([:white, :red, :orange, :yellow]), grid=:off, aspect_ratio=:equal, fill=:true, xlims=(-L, L), clims=(0, 0.04), box=:on, title=L"\hat{h}_{AA}^{12}(q)", tickfontsize=18, legendfontsize=18, guidefontsize=18, titlefontsize=20,bottom_margin=2mm, ticks=[-5, 0, 5], right_margin=8mm)


hftToy = hopToy_inter_strength(Lat)
vtoy = reshape(hftToy.(r), length(xx), length(yy))
P = heatmap(xx, yy, vtoy, size=(700, 600), c=cgrad([:white, :red, :orange, :yellow]), grid=:off, aspect_ratio=:equal, fill=:true, xlims=(-L, L), clims=(0, 0.04), box=:on, title=L"\hat{h}_{12}(q)", tickfontsize=18, legendfontsize=18, guidefontsize=18, titlefontsize=20, ticks=[-5, 0, 5], right_margin=8mm)

τ =6
K1 = Lat.KM[1]
Bτ, indτ, numτ = BtauGen(τ, K1, Lat.latR[1])
G1τ = Bτ * Lat.latR[1]'
Btx = K1[1] .+ G1τ[1:3,1]
Bty = K1[2] .+ G1τ[1:3, 2]
Bx = K1[1] .+ G1τ[:,1]
By = K1[2] .+ G1τ[:, 2]
plot!(P, Bx, By, markersize=6, st=:scatter, color=:white, xlims = (-L,L), ylims = (-L,L),label="", aspect_ratio=:equal)
plot!(P, Btx, Bty, markersize=6, st=:scatter, color=:black, label="", aspect_ratio=:equal)
savefig("inter_coup_toy.pdf")

using KrylovKit
using GeneralizedBM
using Plots
using LinearAlgebra

θ = 1. # twist angle 
rcut = 50.0 # cutoff of the basis

# define the TBL model
Lat = TBLG(θ)
basis = Basis(rcut, Lat);
q = (Lat.KM[1] + Lat.KM[2])/2

M = collect(1:10)
e = []
for m in M
	hop = hopGBM(Lat; Pinter=m)
	@time Hms = hamInter_MST(Lat, basis, hop, q)
	@time Hbm = hamInter_GBM(Lat, basis, hop, q)
	push!(e, norm(Hms-Hbm, 2))
end

plot(M,log.(e), label="inter error")
plot!(M, log(norm(Lat.latM, 2)) * M,label="")
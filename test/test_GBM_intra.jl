using KrylovKit
using GeneralizedBM
using Plots
using LinearAlgebra

θ = 1. # twist angle 
rcut = 30.0 # cutoff of the basis

# define the TBL model
Lat = TBLG(θ)
basis = Basis(rcut, Lat);
q = (Lat.KM[1] + Lat.KM[2])/2

M = collect(1:10)
e = []
for m in M
	hop = hopGBM(Lat; Pintra=m)
	@time Hms = hamIntra_MS(Lat, basis, hop, q)
	@time Hbm = hamIntra_GBM(Lat, basis, hop, q)
	push!(e, norm(Hms-Hbm, Inf))
end

plot(M,log.(e), label="intra error")
plot!(M, log(norm(Lat.latM, Inf)) * M,label="")

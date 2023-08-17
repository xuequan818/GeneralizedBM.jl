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
hop = hopGBM(Lat; τ = 4)
@time Hms = hamInter_MST(Lat, basis, hop, q)
Ems = band(Hms, 3)


M = collect(1:10)
e = []
for m in M
    hop = hopGBM(Lat; Pinter=m, τ=4)
	@time Hbm = hamInter_GBM(Lat, basis, hop, q)
    Ebm = band(Hbm, 3)
    push!(e, norm(Ems - Ebm, Inf))
end

plot!(M,log10.(e), label="inter error")
#plot!(M, log(norm(Lat.latM, Inf)) * M, label="")
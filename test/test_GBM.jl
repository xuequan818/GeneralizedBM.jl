using KrylovKit
using GeneralizedBM
using Plots
using LinearAlgebra


θ = 1. # twist angle 
rcut = 60.0 # cutoff of the basis

# define the TBL model
Lat = TBLG(θ)
p1 = 1
p2 = 2
hop = hopGBM(Lat; Pintra=p1, Pinter=p2)
basis = Basis(rcut, Lat);

# generate hamiltonian at momentum q
q = Lat.KM[1]
@time Hms = ham_MST(Lat, basis, hop, q)
@time Hbm = ham_GBM(Lat, basis, hop, q)

# solve the eigen problem
n_eigs = 20
n_E = 4
@time Ems, Ums = eigsolve(Hms, n_eigs, EigSorter(norm; rev=false); krylovdim=n_eigs + 50);
Ems = Ems[1:2nE]
sort!(Ems)
@time Ebm, Ubm = eigsolve(Hbm, n_eigs, EigSorter(norm; rev=false); krylovdim=n_eigs + 50);
Ebm = Ebm[1:2nE]
sort!(Ebm)
e = norm(Ebm - Ems, Inf)
println(" Taylor order = ($(p1), $(p2))  Error = $(e)")
P1 = plot(Ems, label = "MS")
plot!(P1, Ebm, label = "GBM($(p1),$(p2))")


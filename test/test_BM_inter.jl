using KrylovKit
#using GeneralizedBM
using Plots
using LinearAlgebra


θ = 1. # twist angle 
rcut = 100.0 # cutoff of the basis

# define the TBL model
Lat = TBLG(θ)
h = hopBM(Lat)
basis = Basis(rcut, Lat);

# generate hamiltonian at momentum q
q = Lat.KM[1]
@time Hms = hamInter_MS(basis, h, Lat, q)#ham_MS(Lat, basis, h, q);
@time Hbm = hamInter_BM(basis, h, Lat)#ham_MS(Lat, basis, h, q);

# solve the eigen problem
n_eigs = 10
@time Ems, Ums = eigsolve(Hms, n_eigs, EigSorter(norm; rev=false); krylovdim=n_eigs + 50);
Ems = Ems[1:n_eigs]
@show sort!(Ems);
@time Ebm, Ubm = eigsolve(Hbm, n_eigs, EigSorter(norm; rev=false); krylovdim=n_eigs + 50);
Ebm = Ebm[1:n_eigs]
@show sort!(Ebm);
p1 = plot(Ems, label = "MS")
plot!(p1, Ebm, label = "BM")


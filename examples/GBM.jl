using KrylovKit
using GeneralizedBM
using Plots
using LinearAlgebra


θ = 1.1 # twist angle 
rcut = 60. # cutoff of the basis

# define the TBL model
Lat = TBLG(θ);
p1 = 2
p2 = 1
tau = 4
hop = hopGBM(Lat; Pintra=p1, Pinter=p2, τ=tau)
basis = Basis(rcut, Lat);

# generate hamiltonian at momentum q
A = Lat.KM[1]
B = Lat.KM[2]
C = [B[1], B[2] + norm(A - B)]
@time Hbm = ham_GBM(Lat, basis, hop, C);

# solve the eigen problem
n_eigs = 10
@time Ebm, Ubm = eigsolve(Hbm, n_eigs, EigSorter(norm; rev=false); krylovdim=n_eigs + 50);
@show sort!(Ebm)
p1 = plot(Ebm)

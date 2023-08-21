using KrylovKit
using GeneralizedBM
using Plots
using LinearAlgebra


θ = 1.1 # twist angle 
rcut = 30. # cutoff of the basis

# define the TBL model
Lat = TBLG(θ);
p1 = 1
p2 = 0
tau = 1
hop = hopGBM(Lat; Pintra=p1, Pinter=p2, τ=tau)
basis = Basis(rcut, Lat);

# generate hamiltonian at momentum q
A = Lat.KM[1]
B = Lat.KM[2]
C = [B[1], B[2] + norm(A - B)]
@time Hbm = ham_GBM(Lat, basis, hop, A); 

# solve the eigen problem
n_eigs = 5
@time Ebm, Ubm = eigsolve(Hbm, n_eigs, EigSorter(norm; rev=false); krylovdim=n_eigs + 100);
@show sort!(Ebm)
p1 = plot(Ebm)
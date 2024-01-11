using KrylovKit
using GeneralizedBM
using Plots
using LinearAlgebra


θ = 1. # twist angle 
rcut = 20.0 # cutoff of the basis

# define the TBL model
Lat = TBLG(θ;a = 2.46)
hop = hopGBM(Lat)
basis = Basis(rcut, Lat);

# generate hamiltonian at momentum q
A = Lat.KM[1]
B = Lat.KM[2]
C = [B[1], B[2] + norm(A - B)]
@time Hms = ham_MS(Lat, basis, hop, A);

# solve the eigen problem
n_eigs = 10
@time Ems, Ums = eigsolve(Hms, n_eigs, EigSorter(norm; rev=false); krylovdim=n_eigs + 50);
@show sort!(Ems);
p1 = plot(Ems)

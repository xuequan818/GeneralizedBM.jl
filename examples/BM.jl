using KrylovKit
#using GeneralizedBM
using Plots
using LinearAlgebra


θ = 1.1 # twist angle 
rcut = 140. # cutoff of the basis

# define the TBL model
Lat = TBLG(θ; Lz = 0.45);
h = hopBM(Lat; t = 2.6)
basis = Basis(rcut,Lat);

# generate hamiltonian at momentum q
q = Lat.KM[1]
@time Hbm = ham_BM(Lat, basis, h, q);

# solve the eigen problem
n_eigs = 10
@time E, U = eigsolve(Hbm, n_eigs, EigSorter(norm; rev=false); krylovdim=n_eigs + 50);
@show sort!(E)
p1 = plot(E)

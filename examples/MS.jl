using KrylovKit
using GeneralizedBM
using Plots
using LinearAlgebra


θ = 1.1 # twist angle 
rcut = 70.0 # cutoff of the basis

# define the TBL model
Lat = TBLG(θ)
h = hopBM(Lat)
basis = Basis(rcut, Lat);

# generate hamiltonian at momentum q
q = Lat.KM[1]
@time Hms = ham_MS(Lat, basis, h, q);

# solve the eigen problem
n_eigs = 10
@time E, U = eigsolve(Hms, n_eigs, EigSorter(norm; rev=false); krylovdim=n_eigs + 50);
@show sort!(E);
p1 = plot(E)

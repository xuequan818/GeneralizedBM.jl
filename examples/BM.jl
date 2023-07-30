using KrylovKit
using GeneralizedBM
using Plots
using LinearAlgebra


θ = 1.05 # twist angle 
rcut = 200. # cutoff of the basis

# define the TBL model
Lat = TBLG(θ);
h = hopBM(Lat)
basis = Basis(rcut,Lat);

# generate hamiltonian at momentum q
q = Lat.KM[2] 
@time H = hamiltonian(Lat, basis, h, q);

# solve the eigen problem
n_eigs = 100
@time E, U = eigsolve(H, n_eigs, EigSorter(norm; rev=false); krylovdim=n_eigs + 50);
sort!(E)
plot(E)

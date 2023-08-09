using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.1 # twist angle 
rcut = 60.0 # cutoff of the basis

# define the TBL model
Lat = TBLG(θ)
hop = hopGBM(Lat)
#h = hopBM(Lat; Kt=map(x -> Lat.KM[2], 1:3));
basis = Basis(rcut, Lat);

# build the symmetric path (K->Γ->M->K)
A = Lat.KM[1]
B = [A[1] + norm(Lat.KM[1] - Lat.KM[2])*sqrt(3)/2, 0.0]
C = [A[1],0.]

num = 4
qAB = path(A, B, Int(round(sqrt(3) * num)));
qBC = path(B, C, 2num);
qCA = path(C,A,num);
qx = vcat(qAB[1], qBC[1][2:end], qCA[1][2:end-1])
qy = vcat(qAB[2], qBC[2][2:end], qCA[2][2:end-1])

cols = collect(palette(:tab10));
P1 = plot(qAB[1], qAB[2], st=:scatter, aspect_ratio=:equal, xlims=[minimum(qx) - 0.1, maximum(qx) + 0.1], ylims=[minimum(qy) - 0.1, maximum(qy) + 0.1], color=cols[1], label=L"K\to \Gamma")
for (X,i,pt) in zip([qBC,qCA],2:3,[L"\Gamma\to M", L"M\to K"] )
plot!(P1, [X[1]], [X[2]], st=:scatter, color = cols[i], label = pt)
end
P1
P2 = plot(qx, qy, st=:scatter, aspect_ratio=:equal, xlims=[minimum(qx) - 0.1, maximum(qx) + 0.1], ylims=[minimum(qy) - 0.1, maximum(qy) + 0.1])


# generate the band structure
Eq = []
n_eigs = 14
nE = 4
for (q1, q2, i) in zip(qx, qy, 1:length(qx))
    println(" $(i)-th q of $(length(qx)) q-points")
    H = ham_MST(Lat, basis, hop, [q1, q2])
    @time E, U = eigsolve(H, n_eigs, EigSorter(norm; rev=false); krylovdim=n_eigs + 50)
    s = sortperm(E)
    s1 = findfirst(x -> x == 1, s)
    s2 = findfirst(x -> x == 2, s)
    sort!(E)
    a1 = E[s1] 
    a2 = E[s2]
    l1 = a1*a2 < 0 ? findfirst(x->x>0.,E) : (a1 > 0 ? s2 : s1) 
    l1 = findfirst(x -> x > 0.0, E)
    E1 = maximum([minimum([abs(E[l1-j] - E[l1+j-1]), abs(E[l1-j] + E[l1+j-1])]) for j = 1:nE])
    E2 = maximum([minimum([abs(E[l1+1-j] - E[l1+j]), abs(E[l1-j+1] + E[l1+j])]) for j = 1:nE])
    E3 = maximum([minimum([abs(E[l1-j-1] - E[l1+j-2]), abs(E[l1-j-1] + E[l1+j-2])]) for j = 1:nE])
    l2 = findmin([E1, E2, E3])[2]
    l = [l1, l1 + 1, l1 - 1][l2]
    append!(Eq, E[l-nE:l+nE-1])
end
Eq = reshape(Eq, 2nE, length(qx))
Eq = hcat(Eq, Eq[:,1])
pind = [1, length(qAB[1]), length(qAB[1]) + length(qBC[1]) - 1, length(qx)+1]
pname = [L"K", L"\Gamma", L"M", L"K"]
P3 = plot(Eq[1, :], ylims=[-1.1*maximum(abs.(Eq)), 1.1*maximum(abs.(Eq))], ylabel="Energy", guidefontsize=22, color=cols[1], title=L"%$θ^\circ", label="", tickfontsize=20, legendfontsize=20, xticks=(pind, pname),
legend=:topright, grid=:off, box=:on, size=(740, 620), titlefontsize=30, right_margin=3mm, top_margin=3mm, lw = 3)
for i = 2:2nE
	plot!(P3, Eq[i,:],label="", lw = 3)
end
P3 
using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.05 # twist angle 
rcut = 30. # cutoff of the basis
p1 = 1 # intralayer expansion order
p2 = 0 # interlayer expansion order
tau = 1 # interlayer hopping truncation

# define the TBL model
Lat = TBLG(θ);
hop = hopToy(Lat; Pintra=p1, Pinter=p2, τ=tau)
basis = Basis(rcut, Lat);

# build the symmetric path (K->K'->Γ->Γ->K)
A = Lat.KM[1]
B = Lat.KM[2]
C = [B[1],B[2] + norm(A-B)]
D = [A[1] + sqrt(3)/2 * norm(A-B),0.]

num = 10
qAB = path(A,B,num);
qBC = path(B,C,num);
qCD = path(C,D,Int(round(sqrt(3)*num)));
qDA = path(D,A,num);
qx = vcat(qAB[1], qBC[1][2:end], qCD[1][2:end], qDA[1][2:end-1])
qy = vcat(qAB[2], qBC[2][2:end], qCD[2][2:end], qDA[2][2:end-1])

cols = collect(palette(:tab10));
P1 = plot(qAB[1], qAB[2], st=:scatter, aspect_ratio=:equal, xlims=[minimum(qx) - 0.1, maximum(qx) + 0.1], ylims=[minimum(qy) - 0.1, maximum(qy) + 0.1], color=cols[1], label=L"K\to K'")
for (X,i,pt) in zip([qBC,qCD,qDA],2:4,[L"K'\to \Gamma", L"\Gamma\to \Gamma", L"\Gamma\to K"] )
plot!(P1, [X[1]], [X[2]], st=:scatter, color = cols[i], label = pt)
end
P1
P2 = plot(qx, qy, st=:scatter, aspect_ratio=:equal, xlims=[minimum(qx) - 0.1, maximum(qx) + 0.1], ylims=[minimum(qy) - 0.1, maximum(qy) + 0.1])


# generate the band structure
Eq = []
nE = 3
fv = p2 > 0 ? 0.01 : 0.0
for (q1, q2, i) in zip(qx, qy, 1:length(qx))
    println(" $(i)-th q of $(length(qx)) q-points")
    @time H = hamiltonian(Lat, basis, hop, [q1, q2])
    @time E = band(H, nE; fv=fv)
    append!(Eq, E)
end
Eq = reshape(Eq, 2nE, length(qx))
Eq = hcat(Eq, Eq[:,1])
pind = [1, length(qAB[1]), length(qAB[1]) + length(qBC[1]) - 1, length(qAB[1]) + length(qBC[1]) + length(qCD[1])- 2, length(qx)+1]
pname = [L"K", L"K'", L"\Gamma", L"\Gamma", L"K"]
P3 = plot(Eq[1, :], ylims=[-1.1 * maximum(abs.(Eq)), 1.1 * maximum(abs.(Eq))], ylabel="Energy", guidefontsize=22, color=cols[1], title=L"\theta = %$θ^\circ,\,\, P = (%$p1,%$p2,%$tau)", label="", tickfontsize=20, legendfontsize=20, xticks=(pind, pname), legend=:topright, grid=:off, box=:on, size=(740, 620), titlefontsize=30, right_margin=3mm, top_margin=3mm, lw=1)
for i = 2:2nE
	plot!(P3, Eq[i,:],label="", lw = 1)	
end
P3 

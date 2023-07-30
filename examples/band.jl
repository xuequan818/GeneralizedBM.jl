using KrylovKit
using GeneralizedBM
using Plots, Plots.PlotMeasures, LaTeXStrings
using LinearAlgebra

θ = 1.05 # twist angle 
rcut = 200. # cutoff of the basis

# define the TBL model
Lat = TBLG(θ);
h = hopBM(Lat);
#h = hopBM(Lat; Kt=map(x -> Lat.KM[2], 1:3));
basis = Basis(rcut, Lat);
Hinter = hamInter(basis, h, Lat)

# build the symmetric path (K->K'->Γ->Γ->K)
A = Lat.KM[1]
B = Lat.KM[2]
C = [B[1],B[2] + norm(A-B)]
D = [A[1] + sqrt(3)/2 * norm(A-B),0.]

function path(A::Vector{Float64}, B::Vector{Float64}, factor::Int64)

	xx = []
	yy = []
	if abs(A[1]-B[1]) < 1e-8
		yy = collect(range(A[2], B[2], length=factor + 1))
		xx = repeat([A[1]],factor + 1)
	else 
		f(x) = ((B[2] - A[2])/(B[1]-A[1])) * (x-A[1]) + A[2]

		xx = collect(range(A[1], B[1], length=factor + 1))
		yy = f.(xx)
	end

	return xx,yy
end

num = 8
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
n_eigs = 14
nE = 4
for (q1,q2,i) in zip(qx,qy,1:length(qx))
    println(" $(i)-th q of $(length(qx)) q-points")
    Hintra = hamIntra(basis, h, Lat, [q1, q2])
	H = Hinter + Hintra
	@time E, U = eigsolve(H, n_eigs, EigSorter(norm; rev=false); krylovdim=n_eigs + 50);
    sort!(E)
	l = findfirst(x->x>0.,E)
	append!(Eq, E[l-nE:l+nE-1])
end
Eq = reshape(Eq, 2nE, length(qx))
Eq = hcat(Eq, Eq[:,1])
pind = [1, length(qAB[1]), length(qAB[1]) + length(qBC[1]) - 1, length(qAB[1]) + length(qBC[1]) + length(qCD[1])- 2, length(qx)+1]
pname = [L"K", L"K'", L"\Gamma", L"\Gamma", L"K"]
P3 = plot(Eq[1, :], ylims=[-1.1*maximum(abs.(Eq)), 1.1*maximum(abs.(Eq))], ylabel="Energy", guidefontsize=22, color=cols[1], title=L"%$θ^\circ", label="", tickfontsize=20, legendfontsize=20, xticks=(pind, pname),
legend=:topright, grid=:off, box=:on, size=(740, 620), titlefontsize=30, right_margin=3mm, top_margin=3mm, lw = 3)
for i = 2:2nE
	plot!(P3, Eq[i,:],label="", lw = 3)
end
P3 

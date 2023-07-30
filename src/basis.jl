export Basis

"""
rcut : truncation of the reciprocal lattices
Gmax : max index of the reciprocal lattices
G : index of the reciprocal lattices
Gind : sorting index of G
GM : Moire lattices (Θ21*n)
dof : size of Monolayer basis
"""
struct Basis
	rcut::Float64
	Gmax::Int64
	G::Matrix{Int64}
    Gind::Matrix{Int64}
    GM::Matrix{Float64}
	dof::Int64
end

function basisGen(rcut::Float64, Lat::TBLG)
    latM = Lat.latM
	latR = Lat.latR[1]
	
    Gmax = floor(Int, rcut / norm(latR[:, 1]))
    Gmin = -Gmax

	Gt = zeros(Float64, 2)
    Gind = zeros(Int64, 2Gmax + 1, 2Gmax + 1)
	t = zeros(Int64,2)
	G = Int64[]
	l = 1
	for t1 = Gmin : Gmax, t2 = Gmin : Gmax
		t[1] = t1
		t[2] = t2
		mul!(Gt, latR, t)
		if norm(Gt) <= rcut
			append!(G, t)
            @. t = t + Gmax + 1
            Gind[t[1], t[2]] = l

			l += 1
		end
	end

    G = Array(reshape(G, 2, Int(length(G) / 2))')
    GM = G * latM'
    dof = size(G,1)
    println(" rcut = ", rcut, "; Monolayer DOF = ", dof, "; Matrix DOF = ", 2dof * size(Lat.orb, 1))

	return Basis(rcut, Gmax, G, Gind, GM, dof)
end

Basis(rcut::Float64, Lat::TBLG) = basisGen(rcut, Lat)
#-------------------------------------------------------------------------------
#   Assemble truncated momentum space Hamiltonian matrices
#-------------------------------------------------------------------------------
using StaticArrays, SparseArrays

export hamIntra_MS, hamInter_MS, ham_MS

function hamIntra_MS(basis::Basis, h::Hopping, Lat::TBLG, q::Vector{Float64})
	h1 = h.h11.hG
    h2 = h.h22.hG
    GM = basis.GM
	dof = basis.dof

	c1 = sqrt(Lat.latR_UV[1])
    c2 = sqrt(Lat.latR_UV[2])
    orbl = length(Lat.orb)
    N = orbl^2
    j = repeat(collect(1:orbl),orbl)
    i = sort(j)

    indi = Int64[]
    indj = Int64[]
    vals = ComplexF64[]

    G1i = zeros(Int64, N)
    G1j = zeros(Int64, N)
    G2i = zeros(Int64, N)
    G2j = zeros(Int64, N)
    q1 = zeros(Float64, 2)
    q2 = zeros(Float64, 2)
	for k = 1 : dof #loop for the reciprocal lattices
        @views Gk = GM[k, :]
		@. q2 = Gk + q                           
		hv1 = c1 * h1(q2)

        @. q1 = -Gk + q
        hv2 = c2 * h2(q1)

        @. G2i = i + (k - 1) * orbl
        @. G2j = j + (k - 1) * orbl
        @. G1i = i + (k - 1 + dof) * orbl
        @. G1j = j + (k - 1 + dof) * orbl

        append!(indi, G2i, G1i)
        append!(indj, G2j, G1j)
        append!(vals, hv1, hv2)
	end

    return sparse(indi, indj, vals, 2dof * orbl, 2dof * orbl)
end


function hamInter_MS(basis::Basis, h::Hopping, Lat::TBLG, q::Vector{Float64})
	hF12 = h.h12.hFT
    hF21 = h.h21.hFT

    G1 = basis.G1
    G2 = basis.G2
	dof = basis.dof

	c = sqrt(prod(Lat.latR_UV))
	orbl = length(Lat.orb)
    N = orbl^2
	j = repeat(collect(1:orbl),orbl)
    i = sort(j)

    indi = Int64[]
    indj = Int64[]
    vals = ComplexF64[]

	tol = 1e-6
	qmn = zeros(Float64, 2)
    G1i = zeros(Int64, N)
    G1j = zeros(Int64, N)
    G2i = zeros(Int64, N)
    G2j = zeros(Int64, N)
	for k2 = 1:dof # loop for the reciprocal lattices of sheet 2
		for k1 = 1:dof # loop for the reciprocal lattices of sheet1
            @. qmn = q + G1[k1,:] + G2[k2,:]
			v12 = hF12(qmn)
			if norm(v12,Inf) > tol
                @. G2i = i + (k2 - 1) * orbl
                @. G1j = j + (k1 - 1 + dof) * orbl
                @. G1i = i + (k1 - 1 + dof) * orbl
                @. G2j = j + (k2 - 1) * orbl
                v21 = hF21(qmn)

                append!(indi, G2i, G1i)
                append!(indj, G1j, G2j)
                append!(vals, v12, v21)
			end
		end
	end

	@. vals = c * vals

    return sparse(indi, indj, vals, 2dof * orbl, 2dof * orbl)
end

ham_MS(Lat::TBLG, basis::Basis, h::Hopping, q::Vector{Float64}) = hamIntra_MS(basis, h, Lat, q) + hamInter_MS(basis, h, Lat, q)

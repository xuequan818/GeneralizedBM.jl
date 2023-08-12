#-------------------------------------------------------------------------------
#   Assemble truncated momentum space Hamiltonian matrices
#-------------------------------------------------------------------------------
using StaticArrays, SparseArrays

export hamIntra_MS, hamInter_MS, ham_MS, hamInter_MST, ham_MST

function hamIntra_MS(Lat::TBLG, basis::Basis, h::Hopping, q::Vector{Float64})
    h1 = h.h11
    h2 = h.h22
    GM = basis.GM
    dof = basis.dof

    orbl = length(Lat.orb)
    N = orbl^2
    j = repeat(collect(1:orbl), orbl)
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
    hv1 = zeros(ComplexF64, N)
    hv2 = zeros(ComplexF64, N)
    for k = 1:dof #loop for the reciprocal lattices
        @views Gk = GM[k, :]
        @. q2 = Gk + q
        h1(q2, hv1)

        @. q1 = -Gk + q
        h2(q1, hv2)

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


function hamInter_MS(Lat::TBLG, basis::Basis, h::Hopping, q::Vector{Float64})
    hij = h.hij

    G1 = basis.G1
    G2 = basis.G2
    dof = basis.dof

    c = sqrt(prod(Lat.latR_UV))
    orbl = length(Lat.orb)
    N = orbl^2
    j = repeat(collect(1:orbl), orbl)
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
    v12 = zeros(ComplexF64, N)
    v21 = zeros(ComplexF64, N)
    for k2 = 1:dof # loop for the reciprocal lattices of sheet 2
        @. G2i = i + (k2 - 1) * orbl
        @. G2j = j + (k2 - 1) * orbl
        for k1 = 1:dof # loop for the reciprocal lattices of sheet1
            @. qmn = q + G1[k1, :] + G2[k2, :]
            hij(qmn, v12, v21)
            if norm(v12, Inf) > tol
                @. G1j = j + (k1 - 1 + dof) * orbl
                @. G1i = i + (k1 - 1 + dof) * orbl

                append!(indi, G2i, G1i)
                append!(indj, G1j, G2j)
                append!(vals, v12, v21)
            end
        end
    end

    @. vals = c * vals

    return sparse(indi, indj, vals, 2dof * orbl, 2dof * orbl)
end

ham_MS(Lat::TBLG, basis::Basis, h::Hopping, q::Vector{Float64}) = hamIntra_MS(Lat, basis, h, q) + hamInter_MS(Lat, basis, h, q)

function hamInter_MST(Lat::TBLG, basis::Basis, h::Hopping, q::Vector{Float64})
    hij = h.hij
    Btau = h.Bτ
    ltau = size(Btau, 1)

    G = basis.G
    Gind = basis.Gind
    Gmax = basis.Gmax
    G1 = basis.G1
    G2 = basis.G2
    dof = basis.dof

    c = sqrt(prod(Lat.latR_UV))
    orbl = length(Lat.orb)
    N = orbl^2
    j = repeat(collect(1:orbl), orbl)
    i = sort(j)

    indi = Int64[]
    indj = Int64[]
    vals = ComplexF64[]

    qmn = zeros(Float64, 2)
    G1k = zeros(Int64, 2)
    G1i = zeros(Int64, N)
    G1j = zeros(Int64, N)
    G2i = zeros(Int64, N)
    G2j = zeros(Int64, N)
    v12 = zeros(ComplexF64, N)
    v21 = zeros(ComplexF64, N)
    for k2 = 1:dof # loop for the reciprocal lattices of sheet 2
        @views G2k = G[k2, :]
        @. G2i = i + (k2 - 1) * orbl
        @. G2j = j + (k2 - 1) * orbl

        for l = 1:ltau # loop for the corresponding reciprocal lattices of sheet1
            @. G1k = Btau[l, :] - G2k + Gmax + 1
            if maximum(G1k) <= size(Gind, 1) && minimum(G1k) > 0
                k1 = Gind[G1k[1], G1k[2]]
                if k1 > 0
                    @. qmn = q + G1[k1, :] + G2[k2, :]
                    hij(qmn, v12, v21)

                    @. G1j = j + (k1 - 1 + dof) * orbl
                    @. G1i = i + (k1 - 1 + dof) * orbl

                    append!(indi, G2i, G1i)
                    append!(indj, G1j, G2j)
                    append!(vals, v12, v21)
                end
            end
        end
    end

    @. vals = c * vals

    return sparse(indi, indj, vals, 2dof * orbl, 2dof * orbl)
end

ham_MST(Lat::TBLG, basis::Basis, h::Hopping, q::Vector{Float64}) = hamIntra_MS(Lat, basis, h, q) + hamInter_MST(Lat, basis, h, q)
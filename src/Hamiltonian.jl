#-------------------------------------------------------------------------------
#   Assemble truncated momentum space Hamiltonian matrices
#-------------------------------------------------------------------------------
using StaticArrays, SparseArrays

export hamIntra, hamInter, hamiltonian

#---------------------------------------------------------
#   Assemble intarlayer elements
#---------------------------------------------------------
# momentum space for intra
function hamIntra(Lat::TBLG, basis::Basis, h::IntraHopMS, q::Vector{Float64})
    h1 = h.h11
    h2 = h.h22

    G1 = basis.G1
    G2 = basis.G2
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
        @views G1k = G1[k, :]
        @. q2 = Gk + q
        h1(G1k, q2, hv1)

        @. q1 = -Gk + q
        @views G2k = G2[k, :]
        h2(G2k, q1, hv2)

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

# BM model for intra
function hamIntra(Lat::TBLG, basis::Basis, hop::IntraHopGBM, q::Vector{Float64})
    h1 = hop.h11TP
    h2 = hop.h22TP
    Kt = hop.Kt
    q1 = q - Kt[1]
    q2 = q - Kt[2]

    G1 = basis.G1
    G2 = basis.G2
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
    qm1 = zeros(Float64, 2)
    qm2 = zeros(Float64, 2)
    hv1 = zeros(ComplexF64, N)
    hv2 = zeros(ComplexF64, N)
    for k = 1:dof #loop for the reciprocal lattices
        @views GMk = GM[k, :]
        @views G1k = G1[k, :]
        @. qm1 = GMk + q1
        h1(G1k, qm1, hv1)

        @views G2k = G2[k, :]
        @. qm2 = -GMk + q2
        h2(G2k, qm2, hv2)

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

#---------------------------------------------------------
#   Assemble interlayer elements
#---------------------------------------------------------
# momentum space for inter
function hamInter(Lat::TBLG, basis::Basis, h::InterHopMS, q::Vector{Float64})
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

    tol = 1e-14
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
            @views G1k = G1[k1, :]
            @views G2k = G2[k2, :]
            @. qmn = q + G1k + G2k
            hij(G1k, G2k, qmn, v12, v21)
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

# momentum space for inter with hopping truncation
function hamInter(Lat::TBLG, basis::Basis, h::InterHopMST, q::Vector{Float64})
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
        @views G2kr = G2[k2, :]
        @. G2i = i + (k2 - 1) * orbl
        @. G2j = j + (k2 - 1) * orbl

        for l = 1:ltau # loop for the corresponding reciprocal lattices of sheet1
            @. G1k = Btau[l, :] - G2k + Gmax + 1
            if maximum(G1k) <= size(Gind, 1) && minimum(G1k) > 0
                k1 = Gind[G1k[1], G1k[2]]
                if k1 > 0
                    @views G1kr = G1[k1, :]
                    @. qmn = q + G1kr + G2kr
                    hij(G1kr, G2kr, qmn, v12, v21)

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

# BM model for inter 
function hamInter(Lat::TBLG, basis::Basis, hop::InterHopGBM, q::Vector{Float64})
    Kt = hop.Kt
    hij = hop.hijTP
    Btau = hop.Bτ
    ltau = size(Btau, 1)
    qr = q - Kt

    G = basis.G
    GM = basis.GM
    Gind = basis.Gind
    Gmax = basis.Gmax
    dof = basis.dof
    G1 = basis.G1
    G2 = basis.G2

    latM = Lat.latM
    c = sqrt(prod(Lat.latR_UV))
    orbl = length(Lat.orb)
    N = orbl^2
    j = repeat(collect(1:orbl), orbl)
    i = sort(j)

    indi = Int64[]
    indj = Int64[]
    vals = ComplexF64[]

    G1l = zeros(Int64, 2)
    qkt = zeros(Float64, 2)
    qm = zeros(Float64, 2)
    G1i = zeros(Int64, N)
    G1j = zeros(Int64, N)
    G2i = zeros(Int64, N)
    G2j = zeros(Int64, N)
    v12 = zeros(ComplexF64, N)
    v21 = zeros(ComplexF64, N)
    for k2 = 1:dof # loop for the reciprocal lattices of sheet 2       
        @views G2k = G[k2, :]
        @views G2kr = G2[k2, :]
        @. qm = qr + GM[k2, :]
        @. G2i = i + (k2 - 1) * orbl
        @. G2j = j + (k2 - 1) * orbl

        for l = 1:ltau # loop for the corresponding reciprocal lattices of sheet1
            @. G1l = Btau[l, :] - G2k + Gmax + 1
            if maximum(G1l) <= size(Gind, 1) && minimum(G1l) > 0
                k1 = Gind[G1l[1], G1l[2]]
                if k1 > 0
                    @views G1kr = G1[k1, :]
                    hij(G1kr, G2kr, l, qm, v12, v21)

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

hamiltonian(Lat::TBLG, basis::Basis, h::Hopping, q::Vector{Float64}) = hamIntra(Lat, basis, h.intraHop, q) + hamInter(Lat, basis, h.interHop, q)

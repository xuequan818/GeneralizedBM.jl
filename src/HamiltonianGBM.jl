#-------------------------------------------------------------------------------
#   Assemble truncated momentum space Hamiltonian matrices
#-------------------------------------------------------------------------------
export hamIntra_GBM, hamInter_GBM, ham_GBM

function hamIntra_GBM(Lat::TBLG, basis::Basis, hop::Hopping, q::Vector{Float64})
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
        @views GMk = GM[k,:]
        @views G1k = G1[k,:]
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


function hamInter_GBM(Lat::TBLG, basis::Basis, hop::Hopping, q::Vector{Float64})
    Kt = hop.Kt[3]
    hij = hop.hijTP
    Btau = hop.Bτ
    G1tau = hop.G1τ
    G2tau = hop.G2τ
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
                    @. qkt = Kt + G1tau[l, :]
                    @views G1kr = G1[k1, :]
                    hij(G1kr, G2kr, qkt, qm, v12, v21)

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

ham_GBM(Lat::TBLG, basis::Basis, h::Hopping, q::Vector{Float64}) = hamIntra_GBM(Lat, basis, h, q) + hamInter_GBM(Lat, basis, h, q)
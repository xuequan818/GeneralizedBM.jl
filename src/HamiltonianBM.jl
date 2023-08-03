#-------------------------------------------------------------------------------
#   Assemble truncated momentum space Hamiltonian matrices
#-------------------------------------------------------------------------------
export hamIntra_BM, hamInter_BM, ham_BM

function hamIntra_BM(basis::Basis, h::Hopping, Lat::TBLG, q::Vector{Float64})
    hG1 = h.h11.hGv
    dhG1 = h.h11.dhGv
    hG2 = h.h22.hGv
    dhG2 = h.h22.dhGv
    Kt = h.Kt
    q1 = q - Kt[1]
    q2 = q - Kt[2]

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
    GM1 = zeros(Float64, 2)
    GM2 = zeros(Float64, 2)
    dGk1 = zeros(ComplexF64, N)
    dGk2 = zeros(ComplexF64, N)
    hv1 = zeros(ComplexF64, N)
    hv2 = zeros(ComplexF64, N)
    for k = 1:dof #loop for the reciprocal lattices
        @views Gk = GM[k, :]
        @. GM2 = Gk + q1
        mul!(dGk1, dhG1, GM2)
        @. hv1 = hG1 + dGk1

        @. GM1 = -Gk + q2
        mul!(dGk2, dhG2, GM1)
        @. hv2 = hG2 + dGk2

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


function hamInter_BM(basis::Basis, h::Hopping, Lat::TBLG)
    Kt = h.Kt[3]
    hF12 = h.h12.hFT
    hF21 = h.h21.hFT
    Btau = h.Bτ
    G1tau = h.G1τ
    G2tau = h.G2τ
    ltau = size(Btau, 1)

    G = basis.G
    Gind = basis.Gind
    Gmax = basis.Gmax
    dof = basis.dof

    latM = Lat.latM
    c = sqrt(prod(Lat.latR_UV))
    orbl = length(Lat.orb)
    N = orbl^2
    j = repeat(collect(1:orbl), orbl)
    i = sort(j)

    indi = Int64[]
    indj = Int64[]
    vals = ComplexF64[]

    G1 = zeros(Int64, 2)
    qkt = zeros(Float64, 2)
    G1i = zeros(Int64, N)
    G1j = zeros(Int64, N)
    G2i = zeros(Int64, N)
    G2j = zeros(Int64, N)

    for k2 = 1:dof # loop for the reciprocal lattices of sheet 2
        @views G2k = G[k2, :]
        for l = 1:ltau # loop for the corresponding reciprocal lattices of sheet1
            @. G1 = Btau[l] - G2k + Gmax + 1
            if maximum(G1) <= size(Gind, 1) && minimum(G1) > 0
                k1 = Gind[G1[1], G1[2]]
                if k1 > 0
                    @. qkt = Kt + G1tau[l]
                    v12 = hF12(qkt)
                    v21 = hF21(qkt)

                    @. G2i = i + (k2 - 1) * orbl
                    @. G1j = j + (k1 - 1 + dof) * orbl
                    @. G1i = i + (k1 - 1 + dof) * orbl
                    @. G2j = j + (k2 - 1) * orbl

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

ham_BM(Lat::TBLG, basis::Basis, h::Hopping, q::Vector{Float64}) = hamIntra_BM(basis, h, Lat, q) + hamInter_BM(basis, h, Lat)
export hopToyBM, hopTbgBM

function dhGen(h)
    hx(qx, qy) = ForwardDiff.derivative(qx -> h(qx, qy), qx)
    hy(qx, qy) = ForwardDiff.derivative(qy -> h(qx, qy), qy)
    hxx(qx, qy) = ForwardDiff.derivative(qx -> hx(qx, qy), qx)
    hxy(qx, qy) = ForwardDiff.derivative(qy -> hx(qx, qy), qy)
    hyx(qx, qy) = ForwardDiff.derivative(qx -> hy(qx, qy), qx)
    hyy(qx, qy) = ForwardDiff.derivative(qy -> hy(qx, qy), qy)

    return h, hx, hy, hxx, hxy, hyx, hyy
end

function hopping_keep(h_all::Vector, ht_vec; tol=0.007)

    function val_map(f)
        
        #=
        function val_tol(q)
            v = f(q[1], q[2])
            if norm(v) >= tol
                return v
            else
                return 0.0
            end
        end

        [val_tol(ht_vec[i, :]) for i = 1:size(ht_vec, 1)]
        =#
        [f(ht_vec[i, 1], ht_vec[i, 2]) for i = 1:size(ht_vec, 1)]
    end

    
    val = Any[]
    if length(h_all) == 2
        h_full = dhGen(h_all[1])
        g_full = dhGen(h_all[2])
        for (hi, gi) in zip(h_full, g_full)
            fi(qx,qy) = hi(qx,qy) + im * gi(qx,qy)
            append!(val, val_map(fi))
        end
    else
        for hi in dhGen(h_all[1])
            append!(val, val_map(hi))
        end
    end

    l = size(ht_vec, 1)
    mat = reshape(val, l, Int(length(val)/l))
    mat[:,5] .= 0.
    mat[:,6] .= 0.

    mat
end

function inter_bm_tp(orb::Vector{Array{Float64,2}}, P::Int64, htvec::Vector{Vector{Any}}, 
                     G1, G2, q, hval1::Vector{ComplexF64}, hval2::Vector{ComplexF64})

    gval = [htvec[i][1] for i = 1:4]
    for i = 1:4
        for j = 1:P
            hq = htvec[i][2^j : 2^(j+1) - 1]
            if norm(hq) > 0
                for k = 1:j
                    hq = reshape(hq, 2, 2^(j-k))
                    hq = hq' * q
                end
                gval[i] += hq[1] / factorial(j)
            end
        end
    end
            
    hval1[1] = gval[1] * exp(im * (dot(G1, orb[1][:, 1]) - dot(G2, orb[2][:, 1])))
    hval1[2] = gval[2] * exp(im * (dot(G1, orb[1][:, 1]) - dot(G2, orb[2][:, 2])))
    hval1[3] = gval[3] * exp(im * (dot(G1, orb[1][:, 2]) - dot(G2, orb[2][:, 1])))
    hval1[4] = gval[4] * exp(im * (dot(G1, orb[1][:, 2]) - dot(G2, orb[2][:, 2])))
    hval2[1] = conj(hval1[1])
    hval2[2] = conj(hval1[3])
    hval2[3] = conj(hval1[2])
    hval2[4] = conj(hval1[4])

    hval1, hval2
end

function hopToyBM(Lat::TBLG; t=2.6 * 2 / sqrt(3), Pintra=2, Pinter=2, τ=1)
    intraHop = hopToy_intra(Lat, t, Pintra)

    lat = Lat.lat
    latR = Lat.latR
    orb = Lat.orb

    # interlayer hopping
    # h(r, l) = exp(-α*sqrt(r^2+l^2))
    α = 2.0
    hkl(k, l) = (α * exp(-l * sqrt(k[1]^2 + k[2]^2 + α^2)) * (1 + l * sqrt(k[1]^2 + k[2]^2 + α^2)) / (k[1]^2 + k[2]^2 + α^2)^(3 / 2)) / 2pi
    f(l) = Lat.latR_UV[1] * hkl(Lat.KM[1], l) - 0.11
    Lz = find_zeros(f, 0.0, α)[1]
    htoy(qx, qy) = hkl([qx,qy], Lz)
    Bτ, indτ, numτ = BtauGen(τ, Lat.KM[1], latR[1])
    G1τ = Bτ * latR[1]'
    G2τ = Bτ * latR[2]'
    hmat = hopping_keep([htoy], Bτ * Lat.latR[1]' .+ Lat.KM[1]')
    hijTP(G1, G2, qkt, q, hval1, hval2) = inter_bm_tp(orb, Pinter, [hmat[qkt, :] for i = 1:4], G1, G2, q, hval1, hval2)
    interHop = InterHopGBM(Pinter, hijTP, [htoy], Lat.KM[1], τ, numτ, Bτ, G1τ, G2τ)

    return hmat, Hopping(intraHop, interHop)
end

function hopTbgBM(Lat::TBLG; Pintra=nothing, Pinter=nothing, τintra=4, τinter=nothing, nx=100, ny=100, Lrl=8.0, Lft=6.0)
    intraHop = hopTBG_intra(Lat, Pintra, τintra)

    lat = Lat.lat
    latR = Lat.latR
    orb = Lat.orb

    # interlayer hopping
    a = norm(lat[1][:, 1])
    PQv, PQm, J, wt = quad_nodes(Lrl, nx, ny)
    gv = zero(wt)
    phase = zero(wt)
    hftrl, hftim = inter_tbg_ft(a, Lat.θ, PQv, PQm, J, wt, gv, phase)

    Bτ, indτ, numτ = BtauGen(τ, Lat.KM[1], latR[1])
    G1τ = Bτ * latR[1]'
    G2τ = Bτ * latR[2]'
    hvec = Bτ * Lat.latR[1]' .+ Lat.KM[1]'

    function hmatGen()
        hmat = Vector{Any}(undef, 4)
        for i = 1:4
            hirl(qx,qy) = hftrl[i]([qx, qy])
            hiim(qx,qy) = hftim[i]([qx, qy])
            hmat[i] = hopping_keep([hirl, hiim], hvec)
        end
        return hmat
    end
    hmat = hmatGen()

    hijTP(G1, G2, qkt, q, hval1, hval2) = inter_bm_tp(orb, Pinter, [hmat[qkt, :] for i = 1:4], G1, G2, q, hval1, hval2)
    interHop = InterHopGBM(Pinter, hijTP, [hFT], Lat.KM[1], τ, numτ, Bτ, G1τ, G2τ)

    return Hopping(intraHop, interHop)
end
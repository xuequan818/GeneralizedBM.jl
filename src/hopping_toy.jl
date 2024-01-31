export hopToy_intra, hopToy_inter, hopToy

function intraTP(orb::Array{Float64,2}, Frl::Function, Fim::Function, P::Int64, G, qt::Vector{Float64}, q::Vector{Float64}, hval::Vector{ComplexF64})
	
    vrl = Frl(qt)
    vim = Fim(qt)
    for i = 1:P
        # D^i h(qt)*q^i/ i! 
        vrl += TaylorDiff.derivative(Frl, qt, q, i) / factorial(i) # directional derivative of real part
        vim += TaylorDiff.derivative(Fim, qt, q, i) / factorial(i) # directional derivative of imaginary part
    end
	@. hval = 0.0 + 0.0im
	hval[2] = exp(im*dot(G,orb[:,2]-orb[:,1]))*(vrl + im * vim)
	hval[3] = exp(im*dot(G,orb[:,1]-orb[:,2]))*(vrl - im * vim)

	hval
end

function interTP(orb::Vector{Array{Float64,2}}, g::Function, P::Int64, G1, G2, qt::Vector{Float64}, q::Vector{Float64}, hval1::Vector{ComplexF64}, hval2::Vector{ComplexF64})

    gval = g(qt)
    for i = 1:P
        gval += TaylorDiff.derivative(g, qt, q, i) / factorial(i)
    end

    hval1[1] = gval * exp(im * (dot(G1, orb[1][:, 1]) - dot(G2, orb[2][:, 1])))
    hval1[2] = gval * exp(im * (dot(G1, orb[1][:, 1]) - dot(G2, orb[2][:, 2])))
    hval1[3] = gval * exp(im * (dot(G1, orb[1][:, 2]) - dot(G2, orb[2][:, 1])))
    hval1[4] = gval * exp(im * (dot(G1, orb[1][:, 2]) - dot(G2, orb[2][:, 2])))
    hval2[1] = conj(hval1[1])
    hval2[2] = conj(hval1[3])
    hval2[3] = conj(hval1[2])
    hval2[4] = conj(hval1[4])

    hval1, hval2
end

# Bloch transform with orbitals
# hopping sorted as AA, AB, BA, BB
function hopToy_intra(Lat::TBLG, t::Float64, Pintra)
    lat = Lat.lat
    latR = Lat.latR
    orb = Lat.orb

    #F(k, i) = -t * exp(im*dot(k,orb[i][:,2]-orb[i][:,1]) * (1 + exp(-im * dot(k, lat[i][:, 1])) + exp(-im * dot(k, lat[i][:, 2])))
    # intralayer hopping
    Frl(k, i) = -t * (cos(k[1] * (orb[i][1, 1] - orb[i][1, 2]) + k[2] * (orb[i][2, 1] - orb[i][2, 2])) + cos(k[1] * (lat[i][1, 1] + orb[i][1, 1] - orb[i][1, 2]) + k[2] * (lat[i][2, 1] + orb[i][2, 1] - orb[i][2, 2])) + cos(k[1] * (lat[i][1, 2] + orb[i][1, 1] - orb[i][1, 2]) + k[2] * (lat[i][2, 2] + orb[i][2, 1] - orb[i][2, 2])))
    Fim(k, i) = t * (sin(k[1] * (orb[i][1, 1] - orb[i][1, 2]) + k[2] * (orb[i][2, 1] - orb[i][2, 2])) + sin(k[1] * (lat[i][1, 1] + orb[i][1, 1] - orb[i][1, 2]) + k[2] * (lat[i][2, 1] + orb[i][2, 1] - orb[i][2, 2])) + sin(k[1] * (lat[i][1, 2] + orb[i][1, 1] - orb[i][1, 2]) + k[2] * (lat[i][2, 2] + orb[i][2, 1] - orb[i][2, 2])))
	Frl1(k) = Frl(k, 1)
    Frl2(k) = Frl(k, 2)
    Fim1(k) = Fim(k, 1)
    Fim2(k) = Fim(k, 2)
    h11(G, k, hval) = intraTP(orb[1], Frl1, Fim1, 0, G, k, k, hval)
    h22(G, k, hval) = intraTP(orb[2], Frl2, Fim2, 0, G, k, k, hval)
    intraHop = IntraHopMS(h11,h22)
    if Pintra != nothing    
        h11TP(G, q, hval) = intraTP(orb[1], Frl1, Fim1, Pintra, G, Lat.KM[1], q, hval)
        h22TP(G, q, hval) = intraTP(orb[2], Frl2, Fim2, Pintra, G, Lat.KM[2], q, hval)
        intraHop = IntraHopGBM(Pintra, h11TP, h22TP, Lat.KM)
    end

    intraHop
end

function hopToy_inter(Lat::TBLG, τ, Pinter)
    lat = Lat.lat
    latR = Lat.latR
    orb = Lat.orb

    # interlayer hopping
    # h(r, l) = exp(-α*sqrt(r^2+l^2))
    α = 2.0
    hkl(k, l) = (α * exp(-l * sqrt(k[1]^2 + k[2]^2 + α^2)) * (1 + l * sqrt(k[1]^2 + k[2]^2 + α^2)) / (k[1]^2 + k[2]^2 + α^2)^(3 / 2)) / 2pi
    f(l) = Lat.latR_UV[1] * hkl(Lat.KM[1], l) - 0.11
    Lz = find_zeros(f, 0.0, α)[1]
	hFT(k) = hkl(k,Lz) 
    Lat.Lz = Lz
    hij(G1, G2, k, hval1, hval2) = interTP(orb, hFT, 0, G1, G2, k, k, hval1, hval2)
    interHop = InterHopMS([hFT], hij)
    # hopping truncation of interlayer
    if τ != nothing
        Bτ, indτ, numτ = BtauGen(τ, Lat.KM[1], latR[1])
        G1τ = Bτ * latR[1]'
        G2τ = Bτ * latR[2]'
        interHop = InterHopMST([hFT], hij, τ, numτ, Bτ, G1τ, G2τ)
        if Pinter != nothing
            hijTP(G1, G2, qkt, q, hval1, hval2) = interTP(orb, hFT, Pinter, G1, G2, G1τ[qkt, :] + Lat.KM[1], q, hval1, hval2)
            interHop = InterHopGBM(Pinter, hijTP, [hFT], Lat.KM[1], τ, numτ, Bτ, G1τ, G2τ)
        end
    end

    interHop
end

function hopToy(Lat::TBLG; t=2.6 * 2 / sqrt(3), Pintra = nothing, Pinter = nothing, τ = nothing) 
    intraHop = hopToy_intra(Lat, t, Pintra)
    interHop = hopToy_inter(Lat, τ, Pinter)
    
    return Hopping(intraHop, interHop)
end

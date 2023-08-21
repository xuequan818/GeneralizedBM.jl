export Hopping, hopGBM

"""Hopping
Pintra : expansion order of intralayer
Pinter : expansion order of interlayer
hii : intralayer hopping function
hij : interlayer hopping function
hiiTP : Taylor polynomial of intralayer hopping function
hijTP : Taylor polynomial of interlayer hopping function
Kt : Taylor expansion points for h11, h22 and hij, respectively
τ : truncation parameter of interlayer hopping function
Bτ : corresponding truncated basis index
Giτ : corresponding reciprocal lattices of i-sheet
"""
struct Hopping
	Pintra::Int64
	Pinter::Int64
    h11::Function
    h22::Function
    hij::Function
    h11TP::Function
    h22TP::Function
    hijTP::Function
    Kt::Vector{Vector{Float64}}
    τ::Int64
    Bτ::Matrix{Int64}
    G1τ::Matrix{Float64}
    G2τ::Matrix{Float64}
end

function BtauGen(τ::Int64, Kt::Vector{Float64}, latR::Matrix{Float64})

    Bt = zeros(Int64, (2τ + 1)^2, 2)
    l = 1
    for i = -τ:τ, j = -τ:τ
        Bt[l, 1] = i
        Bt[l, 2] = j
        l += 1
    end

    BtG = Bt * latR' .+ Kt'
    dist = round.(sqrt.(BtG[:, 1] .^ 2 + BtG[:, 2] .^ 2); digits=2)
    sp = sortperm(dist)
    sort!(dist)
    dd = unique(dist)
    ind = findfirst(isequal(dd[τ+1]), dist)
    
    return Bt[sp[1:ind-1], :]
end

function intraTP(orb::Array{Float64,2}, Frl::Function, Fim::Function, P::Int64, G, qt::Vector{Float64}, q::Vector{Float64}, hval::Vector{ComplexF64})
	
    vrl = Frl(qt)
    vim = Fim(qt)
    for i = 1:P
        vrl += derivative(Frl, qt, q, i) / factorial(i) # directional derivative of real part
        vim += derivative(Fim, qt, q, i) / factorial(i) # directional derivative of imaginary part
    end
	@. hval = 0.0 + 0.0im
	hval[2] = exp(im*dot(G,orb[:,2]-orb[:,1]))*(vrl + im * vim)
	hval[3] = exp(im*dot(G,orb[:,1]-orb[:,2]))*(vrl - im * vim)

	hval
end

function interTP(orb::Vector{Array{Float64,2}}, g::Function, P::Int64, G1, G2, qt::Vector{Float64}, q::Vector{Float64}, hval1::Vector{ComplexF64}, hval2::Vector{ComplexF64})

    gval = g(qt)
    for i = 1:P
        gval += derivative(g, qt, q, i) / factorial(i)
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

# hopping sorted as AA, AB, BA, BB
function hopGBM(Lat::TBLG, t::Float64, Pintra::Int64, Pinter::Int64, τ::Int64)
    lat = Lat.lat
    latR = Lat.latR
    orb = Lat.orb
    Kt = push!(copy(Lat.KM), Lat.KM[1])

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
    h11TP(G, q, hval) = intraTP(orb[1], Frl1, Fim1, Pintra, G, Kt[1], q, hval)
    h22TP(G, q, hval) = intraTP(orb[2], Frl2, Fim2, Pintra, G, Kt[2], q, hval)

    # interlayer hopping
    # h(r, l) = exp(-α*sqrt(r^2+l^2))
    α = 2.0
    hkl(k, l) = (α * exp(-l * sqrt(k[1]^2 + k[2]^2 + α^2)) * (1 + l * sqrt(k[1]^2 + k[2]^2 + α^2)) / (k[1]^2 + k[2]^2 + α^2)^(3 / 2)) / 2pi
    f(l) = Lat.latR_UV[1] * hkl(Lat.KM[1], l) - 0.11
    Lz = find_zeros(f, 0.0, α)[1]
	hFT(k) = hkl(k,Lz) 
    Lat.Lz = Lz
    hij(G1, G2, k, hval1, hval2) = interTP(orb, hFT, 0, G1, G2, k, k, hval1, hval2)
    hijTP(G1, G2, qt, q, hval1, hval2) = interTP(orb, hFT, Pinter, G1, G2, qt, q, hval1, hval2)

    # hopping truncation of interlayer
    Bτ = BtauGen(τ, Kt[3], latR[1])
    G1τ = Bτ * latR[1]'
    G2τ = Bτ * latR[2]'

    return Hopping(Pintra, Pinter, h11, h22, hij, h11TP, h22TP, hijTP, Kt, τ, Bτ, G1τ, G2τ)
end

hopGBM(Lat::TBLG; t=2.6 * 2 / sqrt(3), Pintra = 1, Pinter = 0, τ = 1) = hopGBM(Lat, t, Pintra, Pinter, τ)

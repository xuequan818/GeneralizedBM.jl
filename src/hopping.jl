#----------------------------------------------------------------------------
# hopping functions: 
# IntraHop for intralayer hopping
# InterHop for interlayer hopping
#----------------------------------------------------------------------------
export Hopping, hopGBM

"""Hopping
Pintra : expansion order of intralayer
Pinter : expansion order of interlayer
hii : intralayer hopping function
hij : interlayer hopping function
hiiTP : Taylor polynomial of intralayer hopping function
hijTP : Taylor polynomial of interlayer hopping function
Kt : Taylor expansion points for h11, h22 and hij, respectively
τ : truncation of hopping function
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
    Bτ::Vector{Vector{Int64}}
    G1τ::Vector{Vector{Float64}}
    G2τ::Vector{Vector{Float64}}
end

# A temporary version that needs to be refined.
function BtauGen(τ::Int64)

    @assert τ == 1

    Btau = [[0, 0], [0, 1], [-1, 0]]
end


function hopTaylor(hrl::Function, him::Function, P::Int64, qt::Vector{Float64}, q::Vector{Float64})
    
	vrl = hrl(qt)
	vim = him(qt)
    for i = 1:P
        vrl += derivative(hrl, qt, q, i) / factorial(i) # directional derivative of real part
        vim += derivative(him, qt, q, i) / factorial(i) # directional derivative of imaginary part
    end

    return vrl, vim
end

function intraTP(Frl::Function, Fim::Function, P::Int64, qt::Vector{Float64}, q::Vector{Float64}, hval::Vector{ComplexF64})
	
	vrl, vim = hopTaylor(Frl, Fim, P, qt, q)
	@. hval = 0.0 + 0.0im
	hval[2] = vrl + im * vim
	hval[3] = vrl - im * vim

	hval
end

function interTP(f::Vector{Function}, g::Function, P::Int64, qt::Vector{Float64}, q::Vector{Float64}, hval1::Vector{ComplexF64}, hval2::Vector{ComplexF64})

	@. hval1 = 0.0 + 0.0im
    @. hval2 = 0.0 + 0.0im

	for i = 1:4
		hrl(x) = g(x) * cos(f[i](x))
        him(x) = g(x) * sin(f[i](x))
        vrl, vim = hopTaylor(hrl, him, P, qt, q)
		hval1[i] = vrl + im * vim
        hval2[i] = vrl - im * vim
	end
	hval2[[2,3]] = hval2[[3,2]]

    hval1, hval2
end

function hopGBM(Lat::TBLG, t::Float64, Kt::Vector{Vector{Float64}}, Pintra::Int64, Pinter::Int64, τ::Int64)
    lat = Lat.lat
    latR = Lat.latR
    orb = Lat.orb

    #F(k, i) = -t * (1 + exp(-im * dot(k, lat[i][:, 1])) + exp(-im * dot(k, lat[i][:, 2])))
    # intralayer hopping
    Frl(k, i) = -t * (1 + cos(k[1] * lat[i][1, 1] + k[2] * lat[i][2, 1]) + cos(k[1] * lat[i][1, 2] + k[2] * lat[i][2, 2]))
    Fim(k, i) = -t * (-sin(k[1] * lat[i][1, 1] + k[2] * lat[i][2, 1]) - sin(k[1] * lat[i][1, 2] + k[2] * lat[i][2, 2]))
	Frl1(k) = Frl(k, 1)
    Frl2(k) = Frl(k, 2)
    Fim1(k) = Fim(k, 1)
    Fim2(k) = Fim(k, 2)
    h11(k, hval) = intraTP(Frl1, Fim1, 0, k, k, hval)
    h22(k, hval) = intraTP(Frl2, Fim2, 0, k, k, hval)
    h11TP(q, hval) = intraTP(Frl1, Fim1, Pintra, Kt[1], q, hval)
    h22TP(q, hval) = intraTP(Frl2, Fim2, Pintra, Kt[2], q, hval)

    # interlayer hopping
    # h(r, l) = exp(-α*sqrt(r^2+l^2))
    α = 2.0
    hkl(k, l) = (α * exp(-l * sqrt(k[1]^2 + k[2]^2 + α^2)) * (1 + l * sqrt(k[1]^2 + k[2]^2 + α^2)) / (k[1]^2 + k[2]^2 + α^2)^(3 / 2)) / 2pi
    f(l) = Lat.latR_UV[1] * hkl(Lat.KM[1], l) - 0.11
    Lz = find_zeros(f, 0.0, α)[1]
	hFT(k) = hkl(k,Lz) 
    Lat.Lz = Lz
	orbf = Vector{Function}(undef, 4)
    orbf[1] = k -> k[1] * (orb[1][1, 1] - orb[2][1, 1]) + k[2] * (orb[1][2, 1] - orb[2][2, 1])
    orbf[2] = k -> k[1] * (orb[1][1, 1] - orb[2][1, 2]) + k[2] * (orb[1][2, 1] - orb[2][2, 2])
    orbf[3] = k -> k[1] * (orb[1][1, 2] - orb[2][1, 1]) + k[2] * (orb[1][2, 2] - orb[2][2, 1])
    orbf[4] = k -> k[1] * (orb[1][1, 2] - orb[2][1, 2]) + k[2] * (orb[1][2, 2] - orb[2][2, 2])
    hij(k, hval1, hval2) = interTP(orbf, hFT, 0, k, k, hval1, hval2)
    hijTP(qt, q, hval1, hval2) = interTP(orbf, hFT, Pinter, qt, q, hval1, hval2)

    Bτ = BtauGen(τ)
    G1τ = map(x -> latR[1] * x, Bτ)
    G2τ = map(x -> latR[2] * x, Bτ)

    return Hopping(Pintra, Pinter, h11, h22, hij, h11TP, h22TP, hijTP, Kt, τ, Bτ, G1τ, G2τ)
end

hopGBM(Lat::TBLG; t=2.6 * 2 / sqrt(3), Kt=push!(copy(Lat.KM), Lat.KM[1]), Pintra = 1, Pinter = 0, τ = 1) = hopGBM(Lat, t, Kt, Pintra, Pinter, τ)

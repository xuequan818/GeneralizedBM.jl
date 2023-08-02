#----------------------------------------------------------------------------
# hopping functions: 
# IntraHop for intralayer hopping
# InterHop for interlayer hopping
#----------------------------------------------------------------------------
export IntraHop, InterHop, Hopping, hopBM

"""IntraHop
Elemets are ordered as AA, AB, BA, BB
hG : Bloch transform of hrl
hGv : value of hG at some point (Dirac point)
dhG : drivative of hG at some point (Dirac point)
"""
struct IntraHop 
	hG::Function
    hGv::Vector{ComplexF64}
    dhGv::Matrix{ComplexF64}
end

function intrahGen(hG::Function, Kt::Vector{Float64})
    hGv = hG(Kt)
	h(k1,k2) = hG([k1,k2])
	dh_dx(k1, k2) = ForwardDiff.derivative(k1 -> h(k1, k2), k1)
    dh_dy(k1, k2) = ForwardDiff.derivative(k2 -> h(k1, k2), k2)
    dhGv = hcat(dh_dx(Kt[1], Kt[2]), dh_dy(Kt[1], Kt[2]))

	IntraHop(hG, hGv, dhGv)
end

"""InterHop
hFT : Fourier transform of h
"""
struct InterHop
    hFT::Function
end

"""Hopping
hii : intralayer hopping function
hij : interlayer hopping function
Kt : taylor expansion points for h11, h22 and hij, respectively
τ : truncation of hopping function
Bτ : corresponding truncated basis index
Giτ : corresponding reciprocal lattices of i-sheet
"""
struct Hopping
    h11::IntraHop
	h22::IntraHop
    h12::InterHop
	h21::InterHop
	Kt::Vector{Vector{Float64}}
	τ::Int64
	Bτ::Vector{Vector{Int64}}
    G1τ::Vector{Vector{Float64}}
    G2τ::Vector{Vector{Float64}}
end

# A temporary version that needs to be refined.
function BtauGen(τ::Int64)
	
	@assert τ == 1 

	Btau = [[0,0], [0,1], [-1,0]]
end

function hopGen(hG11::Function, hG22::Function, hFT12::Function, hFT21::Function, Kt::Vector{Vector{Float64}}, 
	τ::Int64, G1::Matrix{Float64}, G2::Matrix{Float64})

	h11 = intrahGen(hG11, Kt[1])
	h22 = intrahGen(hG22, Kt[2])
	h12 = InterHop(hFT12)
	h21 = InterHop(hFT21)
	
	Bτ = BtauGen(τ)
	G1τ = map(x -> G1*x, Bτ)
	G2τ = map(x -> G2*x, Bτ)

	return Hopping(h11, h22, h12, h21, Kt, τ, Bτ, G1τ, G2τ)
end

function hopBM(Lat::TBLG, t::Float64, Kt::Vector{Vector{Float64}}, α::Float64)
	lat = Lat.lat
	latR = Lat.latR
	orb = Lat.orb
	Lz = Lat.Lz
    c1 = sqrt(Lat.latR_UV[1])
    c2 = sqrt(Lat.latR_UV[2])

	# intralayer hopping
	#F(k,i) = exp(im * dot(k, orb[i][:,2]- orb[i][:,1])) * (1 + exp(-im * dot(k, lat[i][:,1]) + exp(-im * dot(k, lat[i][:,2]))))
    F(k, i) = 1 + exp(-im * dot(k, lat[i][:, 1])) + exp(-im * dot(k, lat[i][:, 2]))
    hG11(k) = -t/c1 .* [0.0, F(k, 1), conj(F(k, 1)), 0.0]
    hG22(k) = -t/c2 .* [0.0, F(k, 2), conj(F(k, 2)), 0.0]

	# interlayer hopping
    # h(r, l) = exp(-α*sqrt(r^2+l^2))
    hFT(k) = (α * exp(-Lz * sqrt(norm(k)^2 + α^2)) * (1 + Lz * sqrt(norm(k)^2 + α^2)) / (norm(k)^2 + α^2)^(3 / 2)) / 2pi
	orb12(k) = [exp(im * dot(k, orb[1][:, 1] - orb[2][:, 1])), exp(im * dot(k, orb[1][:, 1] - orb[2][:, 2])), 
				exp(im * dot(k, orb[1][:, 2] - orb[2][:, 1])), exp(im * dot(k, orb[1][:, 2] - orb[2][:, 2]))]
    hFT12(k) = hFT(k) .* orb12(k)
    orb21(k) = [exp(im * dot(k, orb[2][:, 1] - orb[1][:, 1])), exp(im * dot(k, orb[2][:, 1] - orb[1][:, 2])),
        		exp(im * dot(k, orb[2][:, 2] - orb[1][:, 1])), exp(im * dot(k, orb[2][:, 2] - orb[1][:, 2]))]
    hFT21(k) = hFT(k) .* orb21(k)

	return hopGen(hG11, hG22, hFT12, hFT21, Kt, 1, latR[1], latR[2])
end

hopBM(Lat::TBLG; t = 2 / sqrt(3), Kt = push!(copy(Lat.KM), Lat.KM[1]), α = 1.) = hopBM(Lat, t, Kt, α)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
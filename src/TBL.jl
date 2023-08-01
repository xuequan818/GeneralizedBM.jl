#----------------------------------------------------------------------------
# twisted bilayer structures: 
# TBLG for twisted bilayer graphene
#----------------------------------------------------------------------------

export TwistedBilayer, TBLG

abstract type TwistedBilayer end

"""TBLG
lat : lattice vectors for periodic cells 
latR : lattice vectors for reciprocal lattices
latM : Moire reciprocal lattice
latR_UV : unit cell volume for reciprocal lattice
orb : orbital vectors
KM : positions of Moire dirac points
θ : twist angle (degree)
Lz : distance between two sheets
"""

struct TBLG <: TwistedBilayer
    lat::Vector{Array{Float64,2}}
    latR::Vector{Array{Float64,2}}
    latM::Array{Float64,2}
    latR_UV::Vector{Float64}
    orb::Vector{Array{Float64,2}}
    KM::Vector{Vector{Float64}}
    θ::Float64
    Lz::Float64
end

function TblgGen(latG::Array{Float64,2}, orbG::Array{Float64,2}, KG::Vector{Float64}, θ::Float64, Lz::Float64)
    θ = (2pi * θ) / 360
    R(x) = [cos(x) -sin(x); sin(x) cos(x)]
    R1 = R(-θ/2)
    R2 = R(θ/2)

    lat = [R1 * latG, R2 * latG]
    latRG = 2pi * inv(latG)'
    latR = [R1 * latRG, R2 * latRG]
    latR_UV = repeat([det(latRG)], 2)
    latM = latR[2] - latR[1]
    KM = [R1 * KG, R2 * KG]
    orb = [R1 * orbG, R2 * orbG]

    return TBLG(lat, latR, latM, latR_UV, orb, KM, θ, Lz)
end

TBLG(θ::Float64; Lz = 1., a = 1., latG = a / 2 .* [1. -1.; sqrt(3) sqrt(3)], orbG = [0. 0.; 0. a/sqrt(3)], KG = 4pi/3a .* [1., 0.]) = TblgGen(latG, orbG, KG, θ, Lz)
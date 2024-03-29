export IntraHopping, InterHopping, Hopping, BtauGen

"""IntraHopping
Pintra : expansion order of intralayer
hii : intralayer hopping function
hiiTP : Taylor polynomial of intralayer hopping function
Kt : Taylor expansion points for h11, h22 respectively
"""
abstract type IntraHopping end

struct IntraHopMS <: IntraHopping
    h11::Function
    h22::Function
end

struct IntraHopGBM <: IntraHopping
    Pintra::Int64
    h11TP::Function
    h22TP::Function
    Kt::Vector{Vector{Float64}}
end

"""IntraHopping
Pinter : expansion order of interlayer
hij : interlayer hopping function
hijTP : Taylor polynomial of interlayer hopping function
Kt : Taylor expansion points for hij
τ : truncation parameter of interlayer hopping function
Bτ : corresponding truncated basis index
Giτ : corresponding reciprocal lattices of i-sheet
"""
abstract type InterHopping end

struct InterHopMS <: InterHopping
    interFT::AbstractVector{Any}
    hij::Function
end

struct InterHopMST <: InterHopping
    interFT::AbstractVector{Any}
    hij::Function
    τ::Int64
    numτ::Vector{Int64}
    Bτ::Matrix{Int64}
    G1τ::Matrix{Float64}
    G2τ::Matrix{Float64}
end

struct InterHopGBM <: InterHopping
    Pinter::Int64
    hijTP::Function
    interFT::AbstractVector{Any}
    Kt::Vector{Float64}
    τ::Int64
    numτ::Vector{Int64}
    Bτ::Matrix{Int64}
    G1τ::Matrix{Float64}
    G2τ::Matrix{Float64}
end

struct Hopping
    intraHop::IntraHopping
    interHop::InterHopping
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
    indt = [findfirst(isequal(dd[i]), dist) for i = 1:τ+1]
    numt = indt[2:end] - indt[1:end-1]
    
    return Bt[sp[1:indt[end]-1], :], indt, numt
end
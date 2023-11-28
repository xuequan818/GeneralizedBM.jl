#-------------------------------------------------------------------------------
# hopping functions of TBG
# Bloch transform orbitals
#-------------------------------------------------------------------------------
export hopTBG

#-------------------------------------------------------------------------------
# intra hopping functions setting
#-------------------------------------------------------------------------------
# compF = cos (real) or sin (imaginary)
function graphene_ms(q, compF::Function, lat_hop::Array{Float64,2}, 
                    ind_hop::Vector{Int64}, τ::Int64, t::Vector{Float64})

	phase = compF.(lat_hop*q)

	val = 0.
	for i = 1:τ
		@views phasei = phase[ind_hop[i] : ind_hop[i+1] - 1]
		val += t[i] * sum(phasei)
	end

	return val
end

function intra_tbg_ms(Lat::TBLG; τ = 4)
    lat = Lat.lat
    orb = Lat.orb

    t = [-2.8922, 0.2425, -0.2656, 0.0235, 0.0524, -0.0209, -0.0148, -0.0211]
    hopt, indt = BtauGen(τ, orb[1][:,2], lat[1])
    lat_ht1 = hopt * lat[1]' + ones(size(hopt,1)) * (orb[1][:, 2] - orb[1][:, 1])'
    lat_ht2 = hopt * lat[2]' + ones(size(hopt, 1)) * (orb[2][:, 2] - orb[2][:, 1])'

    F1rl(q) = graphene_ms(q, cos, lat_ht1, indt, τ, t)
    F2rl(q) = graphene_ms(q, cos, lat_ht2, indt, τ, t)
    F1im(q) = graphene_ms(q, sin, lat_ht1, indt, τ, t)
    F2im(q) = graphene_ms(q, sin, lat_ht2, indt, τ, t)

    return F1rl, F2rl, F1im, F2im
end

# Taylor expansion
function intra_tbg_tp(orb::Array{Float64,2}, Frl::Function, Fim::Function, P::Int64, G, qt::Vector{Float64}, q, hval::Vector{ComplexF64})

    vrl = Frl(qt)
    vim = Fim(qt)
    for i = 1:P
        vrl += TaylorDiff.derivative(Frl, qt, q, i) / factorial(i) # directional derivative of real part
        vim += TaylorDiff.derivative(Fim, qt, q, i) / factorial(i) # directional derivative of imaginary part
    end
    @. hval = 0.0 + 0.0im
    hval[2] = exp(im * dot(G, orb[:, 2] - orb[:, 1])) * (vrl + im * vim)
    hval[3] = exp(im * dot(G, orb[:, 1] - orb[:, 2])) * (vrl - im * vim)

    hval
end

#-------------------------------------------------------------------------------
# inter hopping functions setting
#-------------------------------------------------------------------------------
# layer 1 hops to layer 2
# 1 -> orbital A; 2 -> orbital B
function inter_tbg_rl(a::Float64, θ::Float64, r, orb1::Int64, orb2::Int64)
    rd = sqrt(r[1]^2 + r[2]^2 + 0.0000001)#norm(r)
    R(x) = [cos(x) -sin(x); sin(x) cos(x)]
    r1 = (R(-pi / 6 + θ / 2)*r)[1]
	rs = rd / a
    V0 = 0.3155 * exp(-1.7543 * rs^2) * cos(2.001 * rs)
    V3 = -0.0688 * rs^2 * exp(-3.4692 * (rs - 0.5212)^2)
    V6 = -0.0083 * exp(-2.8764 * (rs - 1.5206)^2) * sin(1.5731 * rs)

    ac = acos(r1 / rd) % (2pi/3)
    theta21 = orb1 == 2 ? ac : ac + pi/3
	theta2 = ac + pi - θ
    theta12 = orb2 == 2 ? theta2 : theta2 + pi/3

    r_cut = 8.
    r_cut2 = 7.
	t = 0.
	if rd < r_cut
		t = V0 + V3 * (cos(3 * theta12) + cos(3 * theta21)) + V6 * (cos(6 * theta12) + cos(6 * theta21))
		if rd > r_cut2
            t = t*exp(1/(r_cut2-r_cut)^2-1/(rd-r_cut)^2);
        end
	end

	return t
end

# divide the region into 3 equal parts.
# compute the numerical integration with Gauss-Legendre points.
function quad_nodes(L::Float64, nx::Int64, ny::Int64)
    
    p_x, w_x = gausslegendre(nx)
    p_y, w_y = gausslegendre(ny)

    w_X = [w_x[i] for i = 1:nx for j = 1:ny]
    w_Y = [w_y[j] for i = 1:nx for j = 1:ny]
    wt = @. w_X * w_Y

    xx = [p_x[i] for i = 1:nx for j = 1:ny]
    yy = [p_y[j] for i = 1:nx for j = 1:ny]

	P(x,y) = L*(2x-y+1)/4#sqrt(3) * L * (1 + x) / 4
    Q(x,y) = sqrt(3) * L * (1 + y) / 4#L * ((1 - x) / 4 + y / 2)
	PQv = [[P(x,y),Q(x,y)] for (x,y) in zip(xx,yy)]
	PQm = hcat(P.(xx,yy),Q.(xx,yy))
    J = sqrt(3) * L^2 / 8
    
	return PQv, PQm, J, wt
end

# compF = cos (real) or -sin (imaginary)
function inter_ft(q, compF::Function, hv::Vector{Float64}, 
                PQm::Matrix{Float64}, J::Float64, wt::Vector{Float64}, 
                gv::Vector{Float64}, phase::Vector{Float64})

    q2 = [-q[1] + sqrt(3) * q[2], -sqrt(3) * q[1] - q[2]] / 2
    q3 = [-q[1] - sqrt(3) * q[2], sqrt(3) * q[1] - q[2]] / 2

    @. gv = 0.0
    for qk in (q, q2, q3)
        mul!(phase, PQm, qk)
        @. gv = gv + compF(phase)
    end
    @. gv = gv * hv

    ftv = J * gv' * wt / (2pi)^2

    return ftv
end

function inter_tbg_ft(a::Float64, θ::Float64, PQv::Vector{Vector{Float64}}, 
                    PQm::Matrix{Float64}, J::Float64, wt::Vector{Float64}, gv::Vector{Float64}, phase::Vector{Float64})

	h11v = inter_tbg_rl.(a, θ, PQv, 1, 1)
	h12v = inter_tbg_rl.(a, θ, PQv, 1, 2)
	h21v = inter_tbg_rl.(a, θ, PQv, 2, 1)
	h22v = inter_tbg_rl.(a, θ, PQv, 2, 2)

    h11ftrl(q) = inter_ft(q, cos, h11v, PQm, J, wt, gv, phase)
    h12ftrl(q) = inter_ft(q, cos, h12v, PQm, J, wt, gv, phase)
    h21ftrl(q) = inter_ft(q, cos, h21v, PQm, J, wt, gv, phase)
    h22ftrl(q) = inter_ft(q, cos, h22v, PQm, J, wt, gv, phase)
    h11ftim(q) = -inter_ft(q, sin, h11v, PQm, J, wt, gv, phase)
    h12ftim(q) = -inter_ft(q, sin, h12v, PQm, J, wt, gv, phase)
    h21ftim(q) = -inter_ft(q, sin, h21v, PQm, J, wt, gv, phase)
    h22ftim(q) = -inter_ft(q, sin, h22v, PQm, J, wt, gv, phase)

    return [h11ftrl, h12ftrl, h21ftrl, h22ftrl], [h11ftim, h12ftim, h21ftim, h22ftim]
end

function inter_tbg_spl(hrl::Vector{Function}, him::Vector{Function}; L=8.0)

    xx = collect(-L:0.1:L)
    yy = copy(xx)
    N = length(xx)
    xy = [[x, y] for x in xx for y in yy]
    xy = reshape(xy, N, N)

    hrl11 = @. hrl[1](xy)
    hrl12 = @. hrl[2](xy)
    hrl21 = @. hrl[3](xy)
    hrl22 = @. hrl[4](xy)
    him11 = @. him[1](xy)
    him12 = @. him[2](xy)
    him21 = @. him[3](xy)
    him22 = @. him[4](xy)

    rlspl11 = Spline2D(xx, yy, hrl11)
    rlspl12 = Spline2D(xx, yy, hrl12)
    rlspl21 = Spline2D(xx, yy, hrl21)
    rlspl22 = Spline2D(xx, yy, hrl22)
    imspl11 = Spline2D(xx, yy, him11)
    imspl12 = Spline2D(xx, yy, him12)
    imspl21 = Spline2D(xx, yy, him21)
    imspl22 = Spline2D(xx, yy, him22)

    h11splrl(q) = rlspl11(q[2], q[1])
    h12splrl(q) = rlspl12(q[2], q[1])
    h21splrl(q) = rlspl21(q[2], q[1])
    h22splrl(q) = rlspl22(q[2], q[1])
    h11splim(q) = imspl11(q[2], q[1])
    h12splim(q) = imspl12(q[2], q[1])
    h21splim(q) = imspl21(q[2], q[1])
    h22splim(q) = imspl22(q[2], q[1])

    return [h11splrl, h12splrl, h21splrl, h22splrl], [h11splim, h12splim, h21splim, h22splim]
end

function inter_tbg_ms(orb::Vector{Array{Float64,2}}, hrl::Vector{Function}, him::Vector{Function},
                     G1, G2, q, hval1::Vector{ComplexF64}, hval2::Vector{ComplexF64})

	@. hval1 = 0.0 + 0.0im
	@. hval2 = 0.0 + 0.0im
	for j = 1:4
        hval1[j] = hrl[j](q) + im * him[j](q)
	end

    hval1[1] *= exp(im * (dot(G1, orb[1][:, 1]) - dot(G2, orb[2][:, 1])))
    hval1[2] *= exp(im * (dot(G1, orb[1][:, 1]) - dot(G2, orb[2][:, 2])))
    hval1[3] *= exp(im * (dot(G1, orb[1][:, 2]) - dot(G2, orb[2][:, 1])))
    hval1[4] *= exp(im * (dot(G1, orb[1][:, 2]) - dot(G2, orb[2][:, 2])))
    hval2[1] = conj(hval1[1])
    hval2[2] = conj(hval1[3])
    hval2[3] = conj(hval1[2])
    hval2[4] = conj(hval1[4])

    return hval1, hval2
end

function inter_tbg_tp(hrl::Vector{Function}, him::Vector{Function}, P::Int64, qt::Vector{Float64})

    function tp_func(h::Function, q)

        v = h(qt)
        for i = 1:P
            v += TaylorDiff.derivative(h, qt, q, i) / factorial(i)
        end

        return v
    end

    h11tprl(q) = tp_func(hrl[1], q)
    h12tprl(q) = tp_func(hrl[2], q)
    h21tprl(q) = tp_func(hrl[3], q)
    h22tprl(q) = tp_func(hrl[4], q)
    h11tpim(q) = tp_func(him[1], q)
    h12tpim(q) = tp_func(him[2], q)
    h21tpim(q) = tp_func(him[3], q)
    h22tpim(q) = tp_func(him[4], q)

    return [h11tprl, h12tprl, h21tprl, h22tprl], [h11tpim, h12tpim, h21tpim, h22tpim]
end

function inter_trunc_tp(hrl::Vector{Function}, him::Vector{Function}, P::Int64, Gtau::Matrix{Float64}, Kt; L = 8.)

    qtk = zeros(2)
    hfunc = Matrix{Any}(undef, size(Gtau,1), 2)
    for k = 1:size(Gtau,1)
        @views Gk = Gtau[k,:]
        @. qtk = Gk + Kt

        htprl, htpim = inter_tbg_tp(hrl, him, P, qtk)
        hsplrl, hsplim = inter_tbg_spl(htprl, htpim; L = L)

        hfunc[k,1] = hsplrl
        hfunc[k,2] = hsplim
    end
    
    return hfunc
end

function hopTBG(Lat::TBLG, Pintra::Int64, 
                Pinter::Int64, τintra::Int64, τinter::Int64;
                nx = 100, ny = 100, Lrl = 5.0, Lft = 8.0)
    lat = Lat.lat
    latR = Lat.latR
    orb = Lat.orb
    Kt = push!(copy(Lat.KM), Lat.KM[1])

    # intralayer hopping
    F1rl, F2rl, F1im, F2im = intra_tbg_ms(Lat; τ=τintra)
    h11(G, k, hval) = intra_tbg_tp(orb[1], F1rl, F1im, 0, G, k, k, hval)
    h22(G, k, hval) = intra_tbg_tp(orb[2], F2rl, F2im, 0, G, k, k, hval)
    h11TP(G, q, hval) = intra_tbg_tp(orb[1], F1rl, F1im, Pintra, G, Kt[1], q, hval)
    h22TP(G, q, hval) = intra_tbg_tp(orb[2], F2rl, F2im, Pintra, G, Kt[2], q, hval)

    # interlayer hopping
    a = norm(lat[1][:,1])
    PQv, PQm, J, wt = quad_nodes(Lrl,nx,ny)
    gv = zero(wt)
    phase = zero(wt)
    hftrl, hftim = inter_tbg_ft(a, Lat.θ, PQv, PQm, J, wt, gv, phase)
    hsplrl, hsplim = inter_tbg_spl(hftrl, hftim; L = Lft)
    hij(G1, G2, q, hval1, hval2) = inter_tbg_ms(orb, hsplrl, hsplim, G1, G2, q, hval1, hval2)

    # hopping truncation of interlayer
    Bτ, indτ = BtauGen(τinter, Kt[3], latR[1])
    G1τ = Bτ * latR[1]'
    G2τ = Bτ * latR[2]'
    hspl = inter_trunc_tp(hftrl, hftim, Pinter, G1τ, Kt[3]; L = Lft)
    hijTP(G1, G2, qkt, q, hval1, hval2) = inter_tbg_ms(orb, hspl[qkt, 1], hspl[qkt, 2], G1, G2, q, hval1, hval2)

    return Hopping(Pintra, Pinter, h11, h22, hij, h11TP, h22TP, hijTP, Kt, τinter, Bτ, G1τ, G2τ)
end

hopTBG(Lat::TBLG; Pintra=1, Pinter=0, τintra=4, τinter=1, nx=100, ny=100, Lrl=8., Lft=6.) = hopTBG(Lat, Pintra, Pinter, τintra, τinter; nx=nx, ny=ny, Lrl=Lrl, Lft = Lft)
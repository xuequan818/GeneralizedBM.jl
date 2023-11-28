#-------------------------------------------------------------------------------
# hopping functions of TBG
# Bloch transform without orbitals
#-------------------------------------------------------------------------------
export hopTBG

# intra hopping functions setting
function graphene_ms_rl(q::Vector{Float64}, lat_hop::Array{Float64,2}, ind_hop::Vector{Int64}, 
						τ::Int64, t::Vector{Float64})

	phase = cos.(lat_hop*q)

	val = 0.
	for i = 1:τ
		@views phasei = phase[ind_hop[i] : ind_hop[i+1] - 1]
		val += t[i] * sum(phasei)
	end

	return val
end

function graphene_ms_im(q::Vector{Float64}, lat_hop::Array{Float64,2}, ind_hop::Vector{Int64},
    τ::Int64, t::Vector{Float64})

    phase = sin.(lat_hop * q)

    val = 0.0
    for i = 1:τ
        @views phasei = phase[ind_hop[i]:ind_hop[i+1]-1]
        val += t[i] * sum(phasei)
    end

    return val
end

function intra_tbg_ms(Lat::TBLG; τ = 4)
    lat = Lat.lat
    orb = Lat.orb

    t = [-2.8922, 0.2425, -0.2656, 0.0235, 0.0524, -0.0209, -0.0148, -0.0211]
    hopt, indt = BtauGen(τ, orb[1][:,2], lat[1])
	lat_ht1 = hopt * lat[1]'
    lat_ht2 = hopt * lat[2]'

    F1rl(q) = graphene_ms_rl(q, lat_ht1, indt, τ, t)
    F2rl(q) = graphene_ms_rl(q, lat_ht2, indt, τ, t)
    F1im(q) = graphene_ms_im(q, lat_ht1, indt, τ, t)
    F2im(q) = graphene_ms_im(q, lat_ht2, indt, τ, t)

    return F1rl, F2rl, F1im, F2im
end

# inter hopping functions setting
# layer 1 hops to layer 2
# 1 -> orbital A; 2 -> orbital B
function inter_tbg_rl(a::Float64, θ::Float64, r::Vector{Float64}, orb1::Int64, orb2::Int64)
	rd = norm(r)
	rs = rd / a
    V0 = 0.3155 * exp(-1.7543 * rs^2) * cos(2.001 * rs)
    V3 = -0.0688 * rs^2 * exp(-3.4692 * (rs - 0.5212)^2)
    V6 = -0.0083 * exp(-2.8764 * (rs - 1.5206)^2) * sin(1.5731 * rs)

    ac = atan(r[2], r[1]) + θ / 2 + pi/6
    theta21 = orb1 == 1 ? ac : ac + pi/3
	theta2 = pi + ac - θ
    theta12 = orb2 == 1 ? theta2 : theta2 + pi/3

    r_cut = 8.
    r_cut2 = 7.
	#t = 0.
	#if rd < r_cut
		t = V0 + V3 * (cos(3 * theta12) + cos(3 * theta21)) + V6 * (cos(6 * theta12) + cos(6 * theta21))
		#if rd > r_cut2
        #    t = t*exp(1/(r_cut2-r_cut)^2-1/(rd-r_cut)^2);
        #end
	#end

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

	P(x,y) = sqrt(3) * L * (1 + x) / 4
    Q(x,y) = L * ((1 - x) / 4 + y / 2)
	PQv = [[P(x,y),Q(x,y)] for (x,y) in zip(xx,yy)]
	PQm = hcat(P.(xx,yy),Q.(xx,yy))
    J = sqrt(3) * L^2 / 8
    
	return PQv, PQm, J, wt
end

function inter_ft_rl(q::Vector{Float64}, hv::Vector{Float64}, PQm::Matrix{Float64}, 
					J::Float64, wt::Vector{Float64})
    
	q2 = [-q[1] + sqrt(3) * q[2], -sqrt(3) * q[1] - q[2]] / 2
    q3 = [-q[1] - sqrt(3) * q[2], sqrt(3) * q[1] - q[2]] / 2

    gv = (cos.(PQm * q) .+ cos.(PQm * q2) .+ cos.(PQm * q3)) .* hv
	ftv = J * dot(gv,wt) / (2pi)^2

	return ftv
end

function inter_ft_im(q::Vector{Float64}, hv::Vector{Float64}, PQm::Matrix{Float64},
    J::Float64, wt::Vector{Float64})

    q2 = [-q[1] + sqrt(3) * q[2], -sqrt(3) * q[1] - q[2]] / 2
    q3 = [-q[1] - sqrt(3) * q[2], sqrt(3) * q[1] - q[2]] / 2

    gv = (sin.(PQm * q) .+ sin.(PQm * q2) .+ sin.(PQm * q3)) .* hv
    ftv = - J * dot(gv, wt) / (2pi)^2

    return ftv
end

function inter_tbg_ms(Lat::TBLG, PQv::Vector{Vector{Float64}}, PQm::Matrix{Float64}, J::Float64, wt::Vector{Float64}, L::Float64)
	a = norm(Lat.lat[1][:,1])
	θ = Lat.θ
	L = L * 2pi / a

	h11v = inter_tbg_rl.(a, θ, PQv, 1, 1)
	h12v = inter_tbg_rl.(a, θ, PQv, 1, 2)
	h21v = inter_tbg_rl.(a, θ, PQv, 2, 1)
	h22v = inter_tbg_rl.(a, θ, PQv, 2, 2)

    h11ftrl(q) = inter_ft_rl(q, h11v, PQm, J, wt)
    h12ftrl(q) = inter_ft_rl(q, h12v, PQm, J, wt)
    h21ftrl(q) = inter_ft_rl(q, h21v, PQm, J, wt)
    h22ftrl(q) = inter_ft_rl(q, h22v, PQm, J, wt)
    h11ftim(q) = inter_ft_im(q, h11v, PQm, J, wt)
    h12ftim(q) = inter_ft_im(q, h12v, PQm, J, wt)
    h21ftim(q) = inter_ft_im(q, h21v, PQm, J, wt)
    h22ftim(q) = inter_ft_im(q, h22v, PQm, J, wt)

    return [h11ftrl, h12ftrl, h21ftrl, h22ftrl], [h11ftim, h12ftim, h21ftim, h22ftim]
end

function intra_tbg_tp(Frl::Function, Fim::Function, P::Int64, qt::Vector{Float64}, q::Vector{Float64}, hval::Vector{ComplexF64})

    vrl = Frl(qt)
    vim = Fim(qt)
    for i = 1:P
        vrl += derivative(Frl, qt, q, i) / factorial(i) # directional derivative of real part
        vim += derivative(Fim, qt, q, i) / factorial(i) # directional derivative of imaginary part
    end
    @. hval = 0.0 + 0.0im
    hval[2] = vrl + im * vim
    hval[3] = vrl - im * vim

    hval
end

function inter_tbg_tp(hrl::Vector{Function}, him::Vector{Function}, P::Int64,qt::Vector{Float64}, q::Vector{Float64}, hval1::Vector{ComplexF64}, hval2::Vector{ComplexF64})
	@. hval1 = 0.0 + 0.0im
	@. hval2 = 0.0 + 0.0im
	for j = 1:4
		hrlj = hrl[j]
		himj = him[j]
    	vrl = hrlj(qt)
		vim = himj(qt)
		for i = 1:P
			vrl += derivative(hrlj, qt, q, i) / factorial(i) # directional derivative of real part
			vim += derivative(himj, qt, q, i) / factorial(i) # directional derivative of imaginary part
		end
		hval1[j] = vrl + im * vim
	end

    hval2[1] = conj(hval1[1])
    hval2[2] = conj(hval1[3])
    hval2[3] = conj(hval1[2])
    hval2[4] = conj(hval1[4])

    hval1, hval2
end

function hopTBG(Lat::TBLG, Pintra::Int64, Pinter::Int64, τintra::Int64, τinter::Int64; nx = 5, ny = 5, L = 5.0)
    lat = Lat.lat
    latR = Lat.latR
    orb = Lat.orb
    Kt = push!(copy(Lat.KM), Lat.KM[1])

    # intralayer hopping
    F1rl, F2rl, F1im, F2im = intra_tbg_ms(Lat; τ=τintra)
    h11(k, hval) = intra_tbg_tp(F1rl, F1im, 0,  k, k, hval)
    h22(k, hval) = intra_tbg_tp(F2rl, F2im, 0,  k, k, hval)
    h11TP(q, hval) = intra_tbg_tp(F1rl, F1im, Pintra,  Kt[1], q, hval)
    h22TP(q, hval) = intra_tbg_tp(F2rl, F2im, Pintra, Kt[2], q, hval)

    # interlayer hopping
    PQv, PQm, J, wt = quad_nodes(L,nx,ny)
    hftrl, hftim = inter_tbg_ms(Lat,PQv,PQm,J,wt,L)
    hij(k, hval1, hval2) = inter_tbg_tp(hftrl, hftim, 0, k, k, hval1, hval2)
    hijTP(qt, q, hval1, hval2) = inter_tbg_tp(hftrl, hftim, Pinter, qt, q, hval1, hval2)

    # hopping truncation of interlayer
    Bτ, indτ = BtauGen(τinter, Kt[3], latR[1])
    G1τ = Bτ * latR[1]'
    G2τ = Bτ * latR[2]'

    return Hopping(Pintra, Pinter, h11, h22, hij, h11TP, h22TP, hijTP, Kt, τinter, Bτ, G1τ, G2τ)
end

hopTBG(Lat::TBLG; Pintra=1, Pinter=0, τintra=4, τinter=1, nx=15, ny=15, L=3.) = hopTBG(Lat, Pintra, Pinter, τintra, τinter; nx=nx, ny=ny, L=L)
export path_local, path, band, band_plot

function path_local(A::Vector{Float64}, B::Vector{Float64}, factor::Int64)

    xx = []
    yy = []
    if abs(A[1] - B[1]) < 1e-8
        yy = collect(range(A[2], B[2], length=factor + 1))
        xx = repeat([A[1]], factor + 1)
    else
        f(x) = ((B[2] - A[2]) / (B[1] - A[1])) * (x - A[1]) + A[2]

        xx = collect(range(A[1], B[1], length=factor + 1))
        yy = f.(xx)
    end

    return xx, yy
end

function path(A::Vector{Float64}, B::Vector{Float64}, fac1::Int64, fac2::Int64, direc::Int64, width ::Int64,a,b)  
    @assert direc in [0, 1, 2]
    xx1, yy1 = path_local(A,B,fac1)
    r1 = collect(range(a,b,length=length(xx1)))

    r = copy(r1)
    xx = copy(xx1)
    yy = copy(yy1)
    if direc == 1
        newA = [xx1[1],yy1[1]]
        newB = [xx1[1+width], yy1[1+width]]
        deleteat!(xx1,1:1+width)
        deleteat!(yy1,1:1+width)
        xx2, yy2 = path_local(newA, newB, fac2)
        r2 = collect(range(r1[1], r1[1+width], length=length(xx2)))
        deleteat!(r1, 1:1+width)
        r = vcat(r2, r1)
        xx = vcat(xx2, xx1)
        yy = vcat(yy2, yy1)
    elseif direc == 2
        l = length(xx1)
        newA = [xx1[l-width], yy1[l-width]]
        newB = [xx1[l], yy1[l]]
        deleteat!(xx1, l-width:l)
        deleteat!(yy1, l-width:l)
        xx2, yy2 = path_local(newA, newB, fac2)
        r2 = collect(range(r1[l-width], r1[l], length=length(xx2)))
        deleteat!(r1, l-width:l)
        r = vcat(r1, r2)
        xx = vcat(xx1, xx2)
        yy = vcat(yy1, yy2)
    end

    return [xx, yy], r    
end

function band(H, nE::Int64; fv = 0.01, n_eigs=2nE + 8)
    g(x) = abs(fv - x)

    E, U = eigsolve(H, n_eigs, EigSorter(g; rev=false); krylovdim=n_eigs + 50)
    sort!(E)
    l1 = findfirst(x -> x >= fv, E)
    E1 = abs.([(E[l1-j] + E[l1+j-1]) / 2 for j = 1:nE] .- 2fv)
    E2 = abs.([(E[l1+1-j] + E[l1+j]) / 2 for j = 1:nE] .- 2fv)
    E3 = abs.([(E[l1-j-1] + E[l1+j-2]) / 2 for j = 1:nE] .- 2fv)
    l2 = findmin(sum.([E1, E2, E3]))[2]
    l = [l1, l1 + 1, l1 - 1][l2]

    return E[l-nE:l+nE-1]
end

function band_plot(Lat, basis, hop, fv;num=10, nE=6, width=2)

    # build the symmetric path (K->Γ->M->K)
    A = Lat.KM[1]
    B = [A[1] + norm(Lat.KM[1] - Lat.KM[2]) * sqrt(3) / 2, 0.0]
    C = [A[1], 0.0]
    num = num
    qAB, rAB = path(A, B, Int(round(sqrt(3) * num)), Int(cld(sqrt(3) * num,2)), 2, width, 0, sqrt(3))
    qBC, rBC = path(B, C, 2num, num, 1, width, sqrt(3), sqrt(3)+2)
    qCA, rCA = path(C, A, num, num, 0, width, sqrt(3) + 2, sqrt(3) + 3)
    qx = vcat(qAB[1], qBC[1][2:end], qCA[1][2:end])
    qy = vcat(qAB[2], qBC[2][2:end], qCA[2][2:end])
    rf = vcat(rAB, rBC[2:end], rCA[2:end])

    # generate the band structure
    Eq = []
    nE = nE
    fv = typeof(hop.interHop) == InterHopGBM ? (hop.interHop.Pinter == 0 ? 0 : fv) : fv
    for (q1,q2,i) in zip(qx,qy,1:length(qx))
        println(" $(i)-th q of $(length(qx)) q-points")
        @time H = hamiltonian(Lat, basis, hop, [q1, q2])
        @time E = band(H, nE; fv = fv)
        append!(Eq, E)
    end
    Eq = reshape(Eq, 2nE, length(qx))

    #plot band
    cols = collect(palette(:tab10))
    pind = [rf[1], rf[length(qAB[1])], rf[length(qAB[1])+length(qBC[1])-1], rf[end]]
    pname = [L"K", L"\Gamma", L"M", L"K"]
    P = plot(rf, Eq[1, :], ylims=[-1.1 * maximum(abs.(Eq)), 1.1 * maximum(abs.(Eq))], ylabel="Energy", xticks=(pind, pname), guidefontsize = 22, color=cols[1], label="", tickfontsize=20, legendfontsize=20, legend=:topright, grid=:off, box=:on, size=(740, 620), titlefontsize=30, right_margin=3mm, top_margin=3mm, lw=1.5)
    for i = 2:2nE
        plot!(P, rf, Eq[i,:],label="", lw = 1.5)
    end

    pind, P
end

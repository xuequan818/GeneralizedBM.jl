export path, band

function path(A::Vector{Float64}, B::Vector{Float64}, factor::Int64)

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

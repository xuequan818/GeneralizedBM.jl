export path

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
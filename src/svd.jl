# SVD Algorithms
# --------------

# Truncated SVD using implicitly restarted Lonczos iterations.
function tsvd(A::Matrix{T}, k::Integer) where T
    m, n = size(A)
    F, = svds(A, nsv=k)
    @static if VERSION > v"0.7-"
        U, Σ, V = F.U, F.S, F.V
    else
        U, Σ, V = F[:U], F[:S], F[:V]
    end
    @assert size(U, 2) == length(Σ) == size(V, 2) == k
    return U, Σ, V
end

# H. Li, G. C. Linderman, et al. "Algorithm 971: An Implementation of a
# Randomized Algorithm for Principal Component Analysis", ACM Transactions on
# Mathematical Software (TOMS), 2017
# DOI: https://doi.org/10.1145/3004053
function rsvd(A::Matrix{T}, k::Integer; its=3, l=k+5) where T
    m, n = size(A)
    l = min(l, m, n)
    @assert 0 < k ≤ l ≤ min(m, n)
    Q::Matrix{T} = rand(eltype(T), n, l) .- T(0.5)
    Y = A*Q
    F = lu!(Y)
    for i in 1:its
        Y = @f A'F.L
        F = lu!(Y)
        Y = @f A*F.L
        if i < its
            F = lu!(Y)
        else
            F = qr!(Y)
        end
    end
    Q = Matrix(@f F.Q)
    B = Q'A
    W, Σ, V = svd(B)
    U = Q*W
    return U[:,1:k], Σ[1:k], V[:,1:k]
end

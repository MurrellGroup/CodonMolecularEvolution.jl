function krylov_sqrt_times_vector(A, v; m=10)
    n = length(v)
    m = min(m, n)  # Can't exceed matrix dimension
    
    # Build Krylov subspace and tridiagonal matrix using Lanczos
    V = zeros(n, m+1)
    T = zeros(m, m)
    
    # Initialize
    V[:,1] = v / norm(v)
    w = A * V[:,1]
    alpha = dot(V[:,1], w)
    w = w - alpha * V[:,1]
    
    T[1,1] = alpha
    
    # Lanczos iteration
    for j in 2:m
        beta = norm(w)
        if beta < 1e-12
            # Early termination - subspace invariant
            return norm(v) * V[:,1:j-1] * sqrt(T[1:j-1,1:j-1]) * [1; zeros(j-2)]
        end
        
        V[:,j] = w / beta
        T[j-1,j] = beta
        T[j,j-1] = beta
        
        w = A * V[:,j]
        alpha = dot(V[:,j], w)
        w = w - alpha * V[:,j] - beta * V[:,j-1]
        
        # Optional reorthogonalization for better stability
        for i in 1:j
            proj = dot(V[:,i], w)
            w = w - proj * V[:,i]
        end
        
        T[j,j] = alpha
    end
    
    # Compute sqrt(T) using standard methods (T is small tridiagonal)
    T_sqrt = sqrt(Symmetric(T))
    
    # Return approximation: ||v|| * V_m * sqrt(T_m) * e_1
    return norm(v) * V[:,1:m] * T_sqrt * [1; zeros(m-1)]
end
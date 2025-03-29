using LinearAlgebra
using BenchmarkTools
using Random

"""
    rational_sqrt_times_vector(A, c, order=8)

Compute A^(1/2)*c using a rational approximation of the square root function.
This avoids computing the full matrix square root explicitly.

Parameters:
    A: Symmetric positive definite matrix
    c: Vector
    order: Order of the rational approximation (higher = more accurate)

Returns:
    Vector approximating A^(1/2)*c
"""
function rational_sqrt_times_vector(A, c, order=8)
    # Scale matrix to have eigenvalues in [0,1] for better approximation
    λ_max = opnorm(A, 2)  # Estimate largest eigenvalue
    A_scaled = A / λ_max
    
    # Compute y = (I - A_scaled)c and z = c
    y = c - A_scaled * c
    z = copy(c)
    
    # Initial approximation x₀ = c
    x = copy(c)
    
    # Apply rational approximation iterations
    for k in 1:order
        # Apply (I + A_scaled)/(2) to improve approximation
        x = (z + A_scaled * x) / 2
    end
    
    # Scale result back
    return sqrt(λ_max) * x
end

"""
    chebyshev_sqrt_times_vector(A, c, m=20)

Compute A^(1/2)*c using a Chebyshev polynomial approximation
This is efficient for a one-time calculation

Parameters:
    A: Symmetric positive definite matrix
    c: Vector
    m: Degree of the Chebyshev polynomial approximation

Returns:
    Vector approximating A^(1/2)*c
"""
function chebyshev_sqrt_times_vector(A, c, m=20)
    # Find spectral bounds - we need eigenvalue range for scaling
    # For efficiency, use power iteration to estimate bounds
    n = size(A, 1)
    
    # Estimate λ_min and λ_max using power iteration
    v = randn(n)
    v = v / norm(v)
    
    # Estimate largest eigenvalue
    λ_max_est = 0.0
    for _ in 1:20
        v = A * v
        v = v / norm(v)
        λ_max_est = dot(v, A * v)
    end
    
    # Estimate smallest eigenvalue using inverse iteration
    # For numerical stability, we'll use a shift
    λ_min_est = λ_max_est
    shift = 1e-6 * λ_max_est
    v = randn(n)
    v = v / norm(v)
    
    # Add a small shift to handle potential singularity
    A_shifted = A + shift * I
    
    for _ in 1:20
        v = A_shifted \ v  # Solve linear system
        v = v / norm(v)
        λ_min_est = dot(v, A * v)
    end
    
    # Adjust for the shift
    λ_min_est = max(λ_min_est - shift, 1e-12)
    
    # Map eigenvalue range to [-1,1] for Chebyshev
    a = λ_min_est
    b = λ_max_est
    
    # Scale and shift for mapping [a,b] to [-1,1]
    alpha = (b + a) / 2
    beta = (b - a) / 2
    
    # Initialize Chebyshev iterations
    T_prev = zeros(n)
    T_curr = c
    result = T_curr / 2  # c₀T₀/2
    
    # Compute coefficients for sqrt function on [-1,1]
    coefs = zeros(m+1)
    for i in 0:m
        # Analytical coefficient for sqrt function
        if i == 0
            coefs[i+1] = 2/π
        else
            # Coefficients for even indices are zero
            if i % 2 == 1
                coefs[i+1] = -2/(π * i * (i-1)) * prod(j -> (2*j-1)/(2*j), 1:(i-1)/2)
            end
        end
    end
    
    # Apply Chebyshev recurrence
    for i in 1:m
        # Compute T_i(x)
        T_next = 2 * ((A - alpha*I)/beta) * T_curr - T_prev
        
        # Update result with coefficient
        result += coefs[i+1] * T_next
        
        # Update for next iteration
        T_prev = T_curr
        T_curr = T_next
    end
    
    # Scale by √β to account for change of variables
    return sqrt(beta) * result
end

"""
    krylov_sqrt_times_vector(A, v, m=30)

Compute A^(1/2)*v using the Krylov subspace method
This is generally the most efficient for a one-time calculation

Parameters:
    A: Symmetric positive definite matrix
    v: Vector
    m: Dimension of Krylov subspace (smaller than matrix size)

Returns:
    Vector approximating A^(1/2)*v
"""
function krylov_sqrt_times_vector(A, v, m=15)
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
   @simd for j in 2:m
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

"""
    native_sqrt_times_vector(A, c)

Calculate A^(1/2)*c using Julia's built-in sqrt function for reference.

Parameters:
    A: Symmetric positive definite matrix
    c: Vector

Returns:
    Vector A^(1/2)*c
"""
function native_sqrt_times_vector(A, c)
    return sqrt(A) * c
end

# Run accuracy and performance tests
function run_tests(n=400)
    println("\n=== Testing with $n × $n matrix ===")
    
    # Create a random SPD matrix
    Random.seed!(42)
    A = rand(n, n)
    A = A' * A + 4I  # Make it SPD with good conditioning
    c = rand(n)
    
    # Normalize c for easier comparison
    c = c / norm(c)
    
    # Reference solution
    println("\nComputing reference solution...")
    reference = native_sqrt_times_vector(A, c)
    
    # Test Krylov method
    println("\nTesting Krylov method...")
    krylov_m = min(20, n)  # Krylov subspace dimension
    krylov_result = krylov_sqrt_times_vector(A, c, krylov_m)
    krylov_error = norm(krylov_result - reference) / norm(reference)
    println("Krylov method relative error (m=$krylov_m): $krylov_error")
    
    # Time comparison
    println("\nPerformance comparison:")
    
    println("Native sqrt(A)*v:")
    @btime native_result = native_sqrt_times_vector($A, $c)
    
    println("Krylov method:")
    @btime krylov_result = krylov_sqrt_times_vector($A, $c, $krylov_m)
    
    # println("Chebyshev method:")
    # @time cheb_result = chebyshev_sqrt_times_vector(A, c)
    # cheb_error = norm(cheb_result - reference) / norm(reference)
    # println("Chebyshev method relative error: $cheb_error")
    
    # println("Rational approximation:")
    # @time rational_result = rational_sqrt_times_vector(A, c)
    # rational_error = norm(rational_result - reference) / norm(reference)
    # println("Rational approximation relative error: $rational_error")
    
    return nothing
end


run_tests(400)

    
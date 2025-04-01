function FUBAR_fitEM(con_lik_matrix, iters, conc; verbosity=1)
    verbosity > 0 && println("Step 4: Model fitting.")
    L = size(con_lik_matrix,1)
    LDAθ = weightEM(con_lik_matrix, ones(L)./L, conc = conc, iters = iters);
    return LDAθ
end
#Need to make this match the smoothFUBAR setup with the first argument controlling the method (EM, Gibbs, etc)
function FUBAR(f::FUBARgrid, outpath;
    pos_thresh=0.95, verbosity=1, exports=true, code=MolecularEvolution.universal_code, optimize_branch_lengths=false,
    method = (sampler = :DirichletEM, concentration = 0.5, iterations = 2500), plots = true)
    exports && init_path(outpath)
    θ = FUBAR_fitEM(f.cond_lik_matrix, method.iterations, method.concentration, verbosity=verbosity)
    df_results = FUBAR_tabulate_from_θ(θ, f, outpath, posterior_threshold = pos_thresh, verbosity = verbosity, plots = plots)
    return df_results, (θ = θ, )
end

function alternating_maximize(f, a_bounds, b_bounds; final_tol = 1e-20, max_iters = 50)
    a = sum(a_bounds)/2
    b = sum(b_bounds)/2
    m = -f(a,b)
    for i in 1:3
        a = brents_method_minimize(x -> -f(x,b), a_bounds[1], a_bounds[2], identity, 1e-5)
        b = brents_method_minimize(x -> -f(a,x), b_bounds[1], b_bounds[2], identity, 1e-5)
    end
    m = -f(a,b)
    next_m = -Inf
    for i in 1:max_iters
        a = brents_method_minimize(x -> -f(x,b), a_bounds[1], a_bounds[2], identity, final_tol)
        b = brents_method_minimize(x -> -f(a,x), b_bounds[1], b_bounds[2], identity, final_tol)
        next_m = -f(a,b)
        if abs(next_m - m) < final_tol
            break;
        end
        m = next_m
    end
    return a,b,-m
end

#Generalize this to work with any grid dimensions!
function interpolating_LRS(grid)
    itp = interpolate(grid, BSpline(Cubic(Line(OnGrid()))))
    
    #Null model:
    ab = brents_method_minimize(x -> -itp(x,x), 1, 20, identity, 1e-20)
    LL_null = itp(ab,ab)

    #Alt model:
    a,b,LL_alt = alternating_maximize(itp, (1.0,20.0), (1.0,20.0))
    LRS = 2(LL_alt-LL_null)
    p = 1-cdf(Chisq(1), LRS)
    return p, (p_value = p, LRS = LRS, LL_alt = LL_alt, α_alt = a, β_alt = b, LL_null = LL_null, αβ_null = ab, sel = ifelse(p<0.05,ifelse(b>a, "Positive", "Purifying"), ""))
end

#Frequentist interpolating fixed-effects FUBAR
#This needs to have exports like regular FUBAR - plots etc
#We could plot LL surfaces for each site under selection, like this:
#=
    itp = interpolate((LLmatrix[:,:,190]'), BSpline(Cubic(Line(OnGrid()))));
    x = range(1, 20, length=200)
    y = range(1, 20, length=200)
    z = @. itp(x', y)
    contour(x, y, z, levels=30, color=:turbo, fill=true, linewidth = 0, colorbar = false, size = (400,400))
=#
function FIFE(f::FUBARgrid, outpath; verbosity=1, exports=true)
    exports && init_path(outpath)
    LLmatrix = reshape(f.LL_matrix, length(f.grid_values),length(f.grid_values),:) 
    #Note: dim1 = beta, dim2 = alpha, so we transpose going in:
    stats = [interpolating_LRS(LLmatrix[:,:,i]') for i in 1:size(LLmatrix, 3)]
    df_results = DataFrame([s[2] for s in stats])
    df_results.site = 1:size(LLmatrix, 3)
    df_results.α_alt .= f.grid_function.(df_results.α_alt)
    df_results.β_alt .= f.grid_function.(df_results.β_alt)
    df_results.αβ_null .= f.grid_function.(df_results.αβ_null)
    CSV.write(outpath * "_results.csv", df_results)
    return df_results
end

export FIFE

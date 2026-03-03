# ========== shapley_mfm_standard_no_splitmerge.jl ==========
using BayesianMixtures
const B = BayesianMixtures
using DelimitedFiles, Statistics, Printf, Random

# Paths relative to the repo root
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DATA_IN  = joinpath(ROOT, "data", "raw", "shapley")
const DATA_OUT = joinpath(ROOT, "data", "processed", "julia_run", "MFM")

mkpath(DATA_OUT)  # make sure output dir exists

# --------------------- Data ---------------------
raw = readdlm(joinpath(DATA_IN, "Shapley_galaxy.dat"))
x_all = vec(raw[:, 4]) ./ 1000.0
cut = quantile(x_all, 0.995) # trim top 0.5%
x_all = x_all[x_all .<= cut]
nmax = length(x_all)
# sample_sizes = [40, 50, 100, 200, 400, 500, 1000, 2000, 4000, nmax]
sample_sizes = [40, 120, 400, 1200, 4000]

# Randomly permute once, so subsets are nested but random
seed = 20251014
Random.seed!(seed)

perm = randperm(nmax)
x_perm = x_all[perm]

open(joinpath(DATA_OUT, "x_perm.csv"), "w") do io
    println(io, "x_perm")
    for val in x_perm
        @printf(io, "%.10f\n", val)
    end
end

# --------------------- Prior / MCMC controls ---------------------
# MFM
log_pk_str = "k -> (k <= 30 ? -log(30.0) : -Inf)"

gamma_dir  = 1.0          # Dirichlet weights parameter 
t_max_cap  = 40           # safe cap for occupied clusters
n_total = 200000         # total sweeps
n_burn  = 100000         # burn-in
# n_keep  = 1000           # kept sweeps
n_keep  = n_total - n_burn 

k_max = 30


# --------------------- Main loop ---------------------
for n in sample_sizes
    @printf("\n=== MFM, n = %d ===\n", n)
    x_n = view(x_perm, 1:n)

    opt = B.options(
        "Normal", "MFM", x_n, n_total;
        n_keep = n_keep,
        n_burn = n_burn,
        verbose = true,
        use_hyperprior = true,
        t_max = t_max_cap,
        gamma = gamma_dir,
        log_pk = log_pk_str
    )

    res = B.run_sampler(opt)

    # Posterior on K (true components)
    pk = B.k_posterior(res; upto = t_max_cap)
    ks = collect(1:length(pk))
    last_nz = findlast(p -> p > 0, pk)
    if last_nz !== nothing
        ks, pk = ks[1:last_nz], pk[1:last_nz]
    end

    prefix = joinpath(DATA_OUT, "shapley_mfm_n$(n)")

    open(prefix * "_k_posterior.csv", "w") do io
        println(io, "k,prob")
        for (kk, p) in zip(ks, pk)
            @printf(io, "%d,%.10g\n", kk, p)
        end
    end

    # --------------------- Posterior draws of K (for boxplots) ---------------------
    # Use post–burn-in draws of t (# occupied clusters)
    iters = (n_burn+1):n_total
    t_draws = @view res.t[iters]

    # under prior support K<=30, we should never see t > 30
    maximum(t_draws) > k_max && error("Observed t > $k_max; cannot sample K with support 1..$k_max.")

    # Build p_kt[k,t] = P(K=k | t) under the same (log_pk, gamma) used in the run
    lpk = eval(Meta.parse(log_pk_str))
    log_pk_fn(k) = Base.invokelatest(lpk, k)
    p_kt = B.MFM.p_kt(log_pk_fn, gamma_dir, n, k_max, t_max_cap)

    # Sample K^(s) ~ P(K | t^(s)) for each saved iteration
    K_draws = similar(t_draws)
    for idx in eachindex(t_draws)
        tt = t_draws[idx]
        K_draws[idx] = B.rand_categorical(view(p_kt, :, tt))
    end

    # Write posterior draws of K 
    open(prefix * "_k_draws.csv", "w") do io
        println(io, "iteration,K")
        for (i, kk) in zip(iters, K_draws)
            @printf(io, "%d,%d\n", i, kk)
        end
    end
end
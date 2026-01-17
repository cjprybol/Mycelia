const EXP_FAMILY_PCA_ENABLED = false  # ExpFamilyPCA removed due to circular deps and poor perf

"""
    sanity_check_matrix(M::AbstractMatrix)

Checks matrix shape, value types, and distributional properties.
Suggests the most appropriate ePCA function and distance metric. When
`EXP_FAMILY_PCA_ENABLED == false`, `suggested_epca` will be `nothing`
for ExpFamilyPCA-backed methods.

Returns a NamedTuple with fields:
- `n_features`, `n_samples`
- `value_type`
- `range`
- `is_binary`
- `is_integer`
- `is_nonnegative`
- `is_strictly_positive`
- `is_in_01`
- `is_centered`
- `is_overdispersed`
- `suggested_epca`
- `suggested_distance`
"""
function sanity_check_matrix(M::AbstractMatrix)
    # Check dimensions
    ndims(M) == 2 || throw(ArgumentError("Input must be a 2D matrix"))
    n_features, n_samples = size(M)
    summary = Dict{Symbol,Any}()
    summary[:n_features] = n_features
    summary[:n_samples] = n_samples

    # Value checks
    elty = eltype(M)
    summary[:value_type] = elty
    summary[:is_integer] = elty <: Integer
    summary[:is_binary] = all(x -> x == 0 || x == 1, M)
    summary[:is_nonnegative] = all(x -> x >= 0, M)
    summary[:is_strictly_positive] = all(x -> x > 0, M)
    summary[:is_in_01] = all(x -> 0 < x < 1, M)
    summary[:range] = (minimum(M), maximum(M))

    # Probability vector check: each column non-negative, sums to 1 (within tolerance)
    function is_probability_vector_matrix(M)
        all(x -> x >= 0, M) && all(abs.(sum(M, dims=1) .- 1) .< 1e-8)
    end
    summary[:is_probability_vector] = is_probability_vector_matrix(M)

    # Centering
    feature_means = mapslices(Statistics.mean, M; dims=2)
    centered = all(abs.(feature_means) .< 1e-6)
    summary[:is_centered] = centered

    # Overdispersion (for count data)
    feature_vars = mapslices(Statistics.var, M; dims=2)
    feature_means_vec = vec(feature_means)
    feature_vars_vec = vec(feature_vars)
    overdispersion = summary[:is_integer] && summary[:is_nonnegative] &&
                     Statistics.mean(feature_vars_vec .- feature_means_vec) > 1.0
    summary[:is_overdispersed] = overdispersion

    # Suggest ePCA and distance
    if summary[:is_binary]
        suggested_epca = :bernoulli_pca_epca
        suggested_distance = :jaccard_distance
    elseif summary[:is_integer] && summary[:is_nonnegative]
        if overdispersion
            suggested_epca = :negbin_pca_epca
            suggested_distance = :bray_curtis_distance
        else
            suggested_epca = :poisson_pca_epca
            suggested_distance = :bray_curtis_distance
        end
    elseif summary[:is_probability_vector]
        suggested_epca = nothing  # No direct EPCA for probability vectors
        suggested_distance = :jensen_shannon_divergence
    elseif summary[:is_in_01]
        suggested_epca = :contbernoulli_pca_epca
        suggested_distance = :cosine_distance
    elseif summary[:is_strictly_positive]
        suggested_epca = :gamma_pca_epca
        suggested_distance = :cosine_distance
    elseif centered
        suggested_epca = :gaussian_pca_epca
        suggested_distance = :euclidean_distance
    else
        suggested_epca = :pca_transform
        suggested_distance = :euclidean_distance
    end
    summary[:suggested_epca_raw] = suggested_epca
    summary[:suggested_epca] = EXP_FAMILY_PCA_ENABLED || suggested_epca == :pca_transform ? suggested_epca : nothing
    summary[:exp_family_pca_enabled] = EXP_FAMILY_PCA_ENABLED
    summary[:suggested_distance] = suggested_distance

    # Print warnings for assumption violations
    if suggested_epca == :poisson_pca_epca && overdispersion
        @warn "Data appears overdispersed (variance > mean); consider using negbin_pca_epca."
    elseif suggested_epca == :negbin_pca_epca && !overdispersion
        @warn "Data does not appear overdispersed (variance ≈ mean); consider using poisson_pca_epca."
    elseif suggested_epca in [:gaussian_pca_epca, :pca_transform] && !centered
        @warn "Data is not centered (mean ≠ 0); consider centering before PCA."
    end

    return summary
end

"""
    pca_transform(
      M::AbstractMatrix{<:Real};
      k::Int = 0,
      var_prop::Float64 = 1.0
    )

Perform standard PCA on `M` (features × samples), returning enough PCs
to either:

- match a user‐supplied `k > 0`, or  
- explain at least `var_prop` of the total variance (0 < var_prop ≤ 1).  

By default (`k=0, var_prop=1.0`), this will capture **all** variance,
i.e. use `min(n_samples-1, n_features)` components.

# When to use
Use for real-valued, continuous, and approximately Gaussian data. PCA is most suitable when features are linearly related and data is centered and scaled. Not ideal for count, binary, or highly skewed data.

# Returns
A NamedTuple with fields
- `model`    : the fitted `MultivariateStats.PCA` object  
- `scores`   : k×n_samples matrix of PC scores  
- `loadings` : k×n_features matrix of PC loadings  
- `chosen_k` : the number of components actually used
"""
function pca_transform(
  M::AbstractMatrix{<:Real};
  k::Int = 0,
  var_prop::Float64 = 1.0
)
  if any(!isfinite, M)
    throw(ArgumentError("PCA input contains non-finite values (NaN or Inf)."))
  end
  n_feats, n_samps = size(M)
  # Warn if not centered
  feature_means = mapslices(Statistics.mean, M; dims=2)
  if any(abs.(feature_means) .> 1e-6)
    @warn "PCA assumes centered data (mean ≈ 0 for each feature); consider centering before PCA."
  end
  # max possible PCs = full rank of X (columns = samples)
  rank_max = min(n_samps - 1, n_feats)

  # user‐supplied k takes priority
  if k > 0
    chosen_k = min(k, rank_max)
    model = MultivariateStats.fit(
      MultivariateStats.PCA, M; maxoutdim=chosen_k
    )
    Z = MultivariateStats.transform(model, M)           # (n_samples × chosen_k)

  else
    # need to auto‐select by variance proportion
    # fit full‐rank PCA to get all eigenvariances
    full = MultivariateStats.fit(
      MultivariateStats.PCA, M; maxoutdim=rank_max
    )
    vars = MultivariateStats.principalvars(full)        # length = rank_max
    total = MultivariateStats.tvar(full)               # = sum(vars) :contentReference[oaicite:1]{index=1}
    # find smallest k so cum‐var / total ≥ var_prop
    if var_prop < 1.0
      cum = cumsum(vars) ./ total
      idx = findfirst(x -> x ≥ var_prop, cum)
      chosen_k = idx === nothing ? rank_max : idx
    else
      chosen_k = rank_max
    end
    # slice out the first chosen_k components
    Z = MultivariateStats.transform(full, M)[:, 1:chosen_k]  # (n_samples × chosen_k)
    model = full
  end

  return (
    model    = model,
    scores   = transpose(Z),                       # (chosen_k × n_samples)
    loadings = MultivariateStats.projection(model)'[:, 1:chosen_k],  # (chosen_k × n_features)
    chosen_k = chosen_k
  )
end

#=
ExpFamilyPCA-backed transforms are disabled due to dependency issues and
performance concerns. The original implementations are left here for
reference and can be re-enabled once the dependency is restored.

"""
    logistic_pca_epca(M::AbstractMatrix{Bool}; k::Int=0)

Synonym for `bernoulli_pca_epca(M; k=k)`.
"""
function logistic_pca_epca(M::AbstractMatrix{Bool}; k::Int=0)
    bernoulli_pca_epca(M; k=k)
end

"""
  bernoulli_pca_epca(M::AbstractMatrix{Bool}; k::Int=0)

Perform Bernoulli (logistic) EPCA on a 0/1 matrix `M` (features × samples).

# When to use
Use for binary (0/1) data, such as presence/absence or yes/no features.

# Returns
A NamedTuple with
- `model`    : the fitted `ExpFamilyPCA.BernoulliEPCA` object  
- `scores`   : k×n_samples matrix of low‐dimensional sample scores  
- `loadings` : k×n_features matrix of feature loadings  
"""
function bernoulli_pca_epca(M::AbstractMatrix{Bool}; k::Int=0)
  # Assert all values are 0 or 1
  if !all(x -> x == 0 || x == 1, M)
    throw(ArgumentError("Bernoulli EPCA requires all entries to be 0 or 1."))
  end
  n_features, n_samples = size(M)
  if k < 1
    k = min(min(n_samples-1, n_features), 10)
  end
  X = transpose(M)                                 # samples × features
  model = ExpFamilyPCA.BernoulliEPCA(n_features, k)
  A     = ExpFamilyPCA.fit!(model, X)              # returns (n_samples×k)
  scores   = transpose(A)                          # k×n_samples
  loadings = model.V                               # k×n_features
  return (model=model, scores=scores, loadings=loadings)
end

"""
  poisson_pca_epca(M::AbstractMatrix{<:Integer}; k::Int=0)

Perform Poisson EPCA on a count matrix `M` (features × samples).

# When to use
Use for non-negative integer count data, such as raw event or read counts.

# Returns
A NamedTuple with
- `model`    : the fitted `ExpFamilyPCA.PoissonEPCA` object  
- `scores`   : k×n_samples matrix of low‐dimensional sample scores  
- `loadings` : k×n_features matrix of feature loadings  
"""
function poisson_pca_epca(M::AbstractMatrix{<:Integer}; k::Int=0)
  if any(M .< 0)
    throw(ArgumentError("Poisson EPCA requires non-negative integer counts."))
  end
  n_features, n_samples = size(M)
  # Warn if overdispersed
  feature_means = mapslices(Statistics.mean, M; dims=2)
  feature_vars = mapslices(Statistics.var, M; dims=2)
  feature_means_vec = vec(feature_means)
  feature_vars_vec = vec(feature_vars)
  overdispersion = Statistics.mean(feature_vars_vec .- feature_means_vec) > 1.0
  if overdispersion
    @warn "Poisson EPCA assumes variance ≈ mean; data appears overdispersed (variance > mean). Consider using negbin_pca_epca."
  end
  if k < 1
    k = min(min(n_samples-1, n_features), 10)
  end
  X = transpose(M)
  model = ExpFamilyPCA.PoissonEPCA(n_features, k)
  A     = ExpFamilyPCA.fit!(model, X)
  scores   = transpose(A)
  loadings = model.V
  return (model=model, scores=scores, loadings=loadings)
end

# ── 1. Negative‐Binomial EPCA ────────────────────────────────────────────────

"""
    negbin_pca_epca(M::AbstractMatrix{<:Integer};
                   k::Int=0,
                   r::Int=1)

Perform Negative-Binomial EPCA on a count matrix `M` (features × samples).

# When to use
Use for overdispersed count data (variance > mean), such as RNA-seq or metagenomic counts.

# Keyword arguments
- `k` : desired number of latent dimensions; if `k<1` defaults to `min(n_samples-1, n_features, 10)`
- `r` : known NB “number of successes” parameter

# Returns
NamedTuple with fields
- `model`    : the fitted `ExpFamilyPCA.NegativeBinomialEPCA` object  
- `scores`   : k×n_samples matrix of sample scores  
- `loadings` : k×n_features matrix of feature loadings  
"""
function negbin_pca_epca(
    M::AbstractMatrix{<:Integer};
    k::Int = 0,
    r::Int = 10
)
    if any(M .< 0)
        throw(ArgumentError("Negative Binomial EPCA requires non-negative integer counts."))
    end
    n_feats, n_samps = size(M)
    # Warn if not overdispersed
    feature_means = mapslices(Statistics.mean, M; dims=2)
    feature_vars = mapslices(Statistics.var, M; dims=2)
    feature_means_vec = vec(feature_means)
    feature_vars_vec = vec(feature_vars)
    overdispersion = Statistics.mean(feature_vars_vec .- feature_means_vec) > 1.0
    if !overdispersion
        @warn "Negative Binomial EPCA assumes overdispersed data (variance > mean); data does not appear overdispersed. Consider using poisson_pca_epca."
    end
    if k < 1
        k = min(min(n_samps-1, n_feats), 10)
    end

    # transpose so each row is a sample
    X = transpose(M)  # (n_samples × n_features)

    # construct and fit NB-EPCA
    model = ExpFamilyPCA.NegativeBinomialEPCA(n_feats, k, r)
    A     = ExpFamilyPCA.fit!(model, X)  # returns (n_samples × k)

    return (
      model    = model,
      scores   = transpose(A),  # (k × n_samples)
      loadings = model.V        # (k × n_features)
    )
end
=#


# ── 4. PCoA from a distance matrix ───────────────────────────────────────────

struct MDSWarningFilterLogger <: Logging.AbstractLogger
    parent::Logging.AbstractLogger
end

Logging.min_enabled_level(logger::MDSWarningFilterLogger) =
    Logging.min_enabled_level(logger.parent)

Logging.shouldlog(logger::MDSWarningFilterLogger, level, _module, group, id) =
    Logging.shouldlog(logger.parent, level, _module, group, id)

Logging.catch_exceptions(logger::MDSWarningFilterLogger) =
    Logging.catch_exceptions(logger.parent)

function Logging.handle_message(
    logger::MDSWarningFilterLogger,
    level,
    message,
    _module,
    group,
    id,
    file,
    line;
    kwargs...
)
    msg = message isa AbstractString ? message : string(message)
    if level == Logging.Warn && occursin("degenerate with", msg)
        return
    end
    return Logging.handle_message(logger.parent, level, message, _module, group, id, file, line; kwargs...)
end

"""
    pcoa_from_dist(D::AbstractMatrix{<:Real}; maxoutdim::Int = 2)

Perform Principal Coordinates Analysis directly from a precomputed
distance matrix `D` (n_samples×n_samples).

# Keyword arguments
- `maxoutdim` : target embedding dimension (default=2)

# Returns
NamedTuple with fields
- `model`      : the fitted `MultivariateStats.MDS` model  
- `coordinates`: maxoutdim×n_samples matrix of embedded coordinates  
"""
function pcoa_from_dist(
    D::AbstractMatrix{<:Real};
    maxoutdim::Int = 3
)
    @assert size(D, 1) == size(D,2) "size(D,1) != size(D,2) $(size(D))"
    model = Logging.with_logger(MDSWarningFilterLogger(Logging.current_logger())) do
        MultivariateStats.fit(
            MultivariateStats.MDS, D; distances=true, maxoutdim=maxoutdim
        )
    end
    Y = MultivariateStats.predict(model)

    return (
      model       = model,
      coordinates = Y
    )
end


"""
$(DocStringExtensions.TYPEDSIGNATURES)

    umap_embed(scores::AbstractMatrix{<:Real};
               n_neighbors::Int=15,
               min_dist::Float64=0.1,
               n_components::Int=2)

Embed your PC/EPCA scores (k×n_samples) into `n_components` via UMAP.

# When to use
Use for visualizing high-dimensional data in 2 or 3 dimensions, especially when the data may have nonlinear structure. UMAP is suitable for both continuous and discrete data, and is robust to non-Gaussian distributions. Input should be a matrix of features or dimensionally-reduced scores (e.g., from PCA or EPCA).

# Arguments
- `scores`      : (components × samples) matrix  
- `n_neighbors` : UMAP neighborhood size  
- `min_dist`    : UMAP min_dist  
- `n_components`: output dimension (2 or 3)

# Returns
- `model`        : the trained UMAP.UMAP model  
"""
function umap_embed(X::AbstractMatrix{<:Real};
                    n_neighbors::Int=Mycelia.nearest_prime(Int(round(log(maximum(size(X)))))),
                    min_dist::Float64=0.1,
                    n_components::Int=2)
    if any(!isfinite, X)
        throw(ArgumentError("UMAP input contains non-finite values (NaN or Inf)."))
    end

    if size(X, 2) < 4096
        dist_matrix = Distances.pairwise(Distances.Euclidean(), X; dims=2)
        model = UMAP.UMAP_(dist_matrix, n_components;
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=:precomputed
        )
    else
        model = UMAP.UMAP_(X, n_components, n_neighbors=n_neighbors, min_dist=min_dist)
    end

    return model
end

#=
Additional ExpFamilyPCA transforms (binomial, continuous Bernoulli,
gamma, Gaussian) were here prior to disabling the dependency. Re-enable
when ExpFamilyPCA is reinstated.

# ── Binomial EPCA ──────────────────────────────────────────────────────────
"""
    binomial_pca_epca(M::AbstractMatrix{<:Integer}; k::Int=0, ntrials::Int=1)

Perform Binomial EPCA on a count matrix `M` (features × samples).

# When to use
Use for integer count data representing the number of successes out of a fixed number of trials (e.g., number of mutated alleles out of total alleles).

# Keyword arguments
- `k` : desired number of latent dimensions; if `k<1` defaults to `min(n_samples-1, n_features, 10)`
- `ntrials` : number of trials for the Binomial distribution (default=1)

# Returns
NamedTuple with fields
- `model`    : the fitted `ExpFamilyPCA.BinomialEPCA` object  
- `scores`   : k×n_samples matrix of sample scores  
- `loadings` : k×n_features matrix of feature loadings  
"""
function binomial_pca_epca(
    M::AbstractMatrix{<:Integer};
    k::Int = 0,
    ntrials::Int = 1
)
    if any(M .< 0) || any(M .> ntrials)
        throw(ArgumentError("Binomial EPCA requires integer counts in 0:ntrials for each entry."))
    end
    n_feats, n_samps = size(M)
    if k < 1
        k = min(min(n_samps-1, n_feats), 10)
    end
    X = transpose(M)  # (n_samples × n_features)
    model = ExpFamilyPCA.BinomialEPCA(n_feats, k, ntrials)
    A     = ExpFamilyPCA.fit!(model, X)
    return (
      model    = model,
      scores   = transpose(A),  # (k × n_samples)
      loadings = model.V        # (k × n_features)
    )
end

# ── Continuous Bernoulli EPCA ───────────────────────────────────────────────
"""
    contbernoulli_pca_epca(M::AbstractMatrix{<:Real}; k::Int=0)

Perform Continuous Bernoulli EPCA on a matrix `M` (features × samples).

# When to use
Use for continuous data in the open interval (0, 1), such as probabilities or normalized intensities.

# Keyword arguments
- `k` : desired number of latent dimensions; if `k<1` defaults to `min(n_samples-1, n_features, 10)`

# Returns
NamedTuple with fields
- `model`    : the fitted `ExpFamilyPCA.ContinuousBernoulliEPCA` object  
- `scores`   : k×n_samples matrix of sample scores  
- `loadings` : k×n_features matrix of feature loadings  
"""
function contbernoulli_pca_epca(
    M::AbstractMatrix{<:Real};
    k::Int = 0
)
    if any(M .<= 0) || any(M .>= 1)
        throw(ArgumentError("Continuous Bernoulli EPCA requires all entries strictly in (0, 1)."))
    end
    n_feats, n_samps = size(M)
    if k < 1
        k = min(min(n_samps-1, n_feats), 10)
    end
    X = transpose(M)  # (n_samples × n_features)
    model = ExpFamilyPCA.ContinuousBernoulliEPCA(n_feats, k)
    A     = ExpFamilyPCA.fit!(model, X)
    return (
      model    = model,
      scores   = transpose(A),  # (k × n_samples)
      loadings = model.V        # (k × n_features)
    )
end

# ── Gamma EPCA ──────────────────────────────────────────────────────────────
"""
    gamma_pca_epca(M::AbstractMatrix{<:Real}; k::Int=0)

Perform Gamma EPCA on a matrix `M` (features × samples).

# When to use
Use for positive continuous data, such as rates, times, or strictly positive measurements.

# Keyword arguments
- `k` : desired number of latent dimensions; if `k<1` defaults to `min(n_samples-1, n_features, 10)`

# Returns
NamedTuple with fields
- `model`    : the fitted `ExpFamilyPCA.GammaEPCA` object  
- `scores`   : k×n_samples matrix of sample scores  
- `loadings` : k×n_features matrix of feature loadings  
"""
function gamma_pca_epca(
    M::AbstractMatrix{<:Real};
    k::Int = 0
)
    if any(M .<= 0)
        throw(ArgumentError("Gamma EPCA requires all entries to be strictly positive."))
    end
    n_feats, n_samps = size(M)
    if k < 1
        k = min(min(n_samps-1, n_feats), 10)
    end
    X = transpose(M)  # (n_samples × n_features)
    model = ExpFamilyPCA.GammaEPCA(n_feats, k)
    A     = ExpFamilyPCA.fit!(model, X)
    return (
      model    = model,
      scores   = transpose(A),  # (k × n_samples)
      loadings = model.V        # (k × n_features)
    )
end

# ── Gaussian EPCA ───────────────────────────────────────────────────────────
"""
    gaussian_pca_epca(M::AbstractMatrix{<:Real}; k::Int=0)

Perform Gaussian EPCA on a matrix `M` (features × samples).

# When to use
Use for real-valued continuous data (centered, can be negative or positive), such as normalized or standardized measurements.

# Keyword arguments
- `k` : desired number of latent dimensions; if `k<1` defaults to `min(n_samples-1, n_features, 10)`

# Returns
NamedTuple with fields
- `model`    : the fitted `ExpFamilyPCA.GaussianEPCA` object  
- `scores`   : k×n_samples matrix of sample scores  
- `loadings` : k×n_features matrix of feature loadings  
"""
function gaussian_pca_epca(
    M::AbstractMatrix{<:Real};
    k::Int = 0
)
    if any(!isfinite, M)
        throw(ArgumentError("Gaussian EPCA input contains non-finite values (NaN or Inf)."))
    end
    n_feats, n_samps = size(M)
    # Warn if not centered
    feature_means = mapslices(Statistics.mean, M; dims=2)
    if any(abs.(feature_means) .> 1e-6)
        @warn "Gaussian EPCA assumes centered data (mean ≈ 0 for each feature); consider centering before use."
    end
    if k < 1
        k = min(min(n_samps-1, n_feats), 10)
    end
    X = transpose(M)  # (n_samples × n_features)
    model = ExpFamilyPCA.GaussianEPCA(n_feats, k)
    A     = ExpFamilyPCA.fit!(model, X)
    return (
      model    = model,
      scores   = transpose(A),  # (k × n_samples)
      loadings = model.V        # (k × n_features)
    )
end
=#

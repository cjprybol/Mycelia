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
  n_feats, n_samps = size(M)
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

"""
    logistic_pca_epca(M::AbstractMatrix{<:Integer}, k::Int)

Perform Bernoulli (logistic) EPCA on a 0/1 matrix `M` (features × samples).

# Returns
A NamedTuple with
- `model`    : the fitted `ExpFamilyPCA.BernoulliEPCA` object  
- `scores`   : k×n_samples matrix of low‐dimensional sample scores  
- `loadings` : k×n_features matrix of feature loadings  
"""
function logistic_pca_epca(M::AbstractMatrix{<:Integer}; k::Int=0)
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
    glm_pca_epca(M::AbstractMatrix{<:Integer}, k::Int)

Perform Poisson EPCA on a count matrix `M` (features × samples).

# Returns
A NamedTuple with
- `model`    : the fitted `ExpFamilyPCA.PoissonEPCA` object  
- `scores`   : k×n_samples matrix of low‐dimensional sample scores  
- `loadings` : k×n_features matrix of feature loadings  
"""
function glm_pca_epca(M::AbstractMatrix{<:Integer}; k::Int=0)
    n_features, n_samples = size(M)
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
    r::Int = 1
)
    n_feats, n_samps = size(M)
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


# ── 4. PCoA from a distance matrix ───────────────────────────────────────────

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
    maxoutdim::Int = 2
)
    @assert size(D, 1) == size(D,2) "size(D,1) != size(D,2) $(size(D))"
    model = MultivariateStats.fit(
      MultivariateStats.MDS, D; distances=true, maxoutdim=maxoutdim
    )
    Y = MultivariateStats.predict(model)

    return (
      model       = model,
      coordinates = Y
    )
end


"""
    umap_embed(scores::AbstractMatrix{<:Real};
               n_neighbors::Int=15,
               min_dist::Float64=0.1,
               n_components::Int=2)

Embed your PC/EPCA scores (k×n_samples) into `n_components` via UMAP.

# Arguments
- `scores`      : (components × samples) matrix  
- `n_neighbors` : UMAP neighborhood size  
- `min_dist`    : UMAP min_dist  
- `n_components`: output dimension (2 or 3)

# Returns
- `embedding` : n_components×n_samples matrix  
- `um`        : the trained UMAP.UMAP model  
"""
function umap_embed(scores::AbstractMatrix{<:Real};
                    n_neighbors::Int=15,
                    min_dist::Float64=0.1,
                    n_components::Int=2)

    # transpose to (samples × components)
    X = transpose(scores)

    # build & fit UMAP
    um = UMAP.UMAP(n_neighbors=n_neighbors,
                   min_dist=min_dist,
                   n_components=n_components)
    embedding = UMAP.fit_transform(um, X)  # returns samples×n_components

    # return embedding in components×samples orientation
    return transpose(embedding), um
end
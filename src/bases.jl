import LinearAlgebra: norm_sqr

@inline lj_potential_uij(
    r::Float64; ϵ = 1.0, σ = 1.0
) = let r² = r ^ 2, σ² = σ ^ 2
    @fastmath 4 * ϵ * ((σ² ^ 6) / (r² ^ 6) - (σ² ^ 3) / (r² ^ 3))
end

@inline lj_potential_uij(
    r; ϵ = 1.0, σ = 1.0
) = let r² = norm_sqr(r), σ² = σ ^ 2
    @fastmath 4 * ϵ * ((σ² ^ 6) / (r² ^ 6) - (σ² ^ 3) / (r² ^ 3))
end

@inline lj_potential_uij(
    r₁, r₂;
    ϵ = 1.0, σ = 1.0,
    dist = default_distance
) = lj_potential_uij(dist(r₁, r₂), ϵ = ϵ, σ = σ)

_lj_potential_wij(
    r², ϵ = 1.0, σ² = 1.0
) = @fastmath 24 * ϵ * (2 * (σ² ^ 6) / (r² ^ 6) - (σ² ^ 3) / (r² ^ 3))

@inline function lj_potential_fij(
    r; ϵ = 1.0, σ = 1.0
)
    r² = norm_sqr(r)
    wij = _lj_potential_wij(r², ϵ, σ ^ 2)
    fij = r * wij / r²
end

@inline lj_potential_fij(
    r₁, r₂;
    ϵ = 1.0, σ = 1.0,
    dist = default_distance
) = lj_potential_fij(dist(r₁, r₂), ϵ = ϵ, σ = σ)

@inline lj_potential_uij_cutoff(
    r;
    ϵ = 1.0, σ = 1.0,
    r_cutoff = 2.5σ
) = norm_sqr(r) ≥ r_cutoff ^ 2 ? 0.0 :
    lj_potential_uij(r, ϵ = ϵ, σ = σ) - lj_potential_uij(r_cutoff, ϵ = ϵ, σ = σ)

@inline lj_potential_uij_cutoff(
    r₁, r₂;
    ϵ = 1.0, σ = 1.0,
    dist = default_distance,
    r_cutoff = 2.5σ
) = lj_potential_uij_cutoff(dist(r₁, r₂), ϵ = ϵ, σ = σ, r_cutoff = r_cutoff)

@inline lj_potential_fij_cutoff(
    r;
    ϵ = 1.0, σ = 1.0,
    r_cutoff = 2.5σ
) = norm_sqr(r) ≥ r_cutoff ^ 2 ? zero(r) :
    lj_potential_fij(r, ϵ = ϵ, σ = σ)

@inline lj_potential_fij_cutoff(
    r₁, r₂;
    ϵ = 1.0, σ = 1.0,
    dist = default_distance,
    r_cutoff = 2.5σ
) = lj_potential_fij_cutoff(dist(r₁, r₂), ϵ = ϵ, σ = σ, r_cutoff = r_cutoff)

# Weeks–Chandler–Anderson potential function, a repulsive potential.

@inline lj_potential_wca_uij(
    r; ϵ = 1.0, σ = 1.0
) = norm_sqr(r) ≥ (R_MIN * σ) ^ 2 ? 0.0 :
    (lj_potential_uij(r, ϵ = ϵ, σ = σ) + ϵ)

@inline lj_potential_wca_uij(
    r₁, r₂;
    ϵ = 1.0, σ = 1.0,
    dist = default_distance
) = lj_potential_wca_uij(dist(r₁, r₂), ϵ = ϵ, σ = σ)

@inline lj_potential_wca_fij(
    r; ϵ = 1.0, σ = 1.0
) = norm_sqr(r) ≥ (R_MIN * σ) ^ 2 ? zero(r) :
    lj_potential_fij(r, ϵ = ϵ, σ = σ)

@inline lj_potential_wca_fij(
    r₁, r₂;
    ϵ = 1.0, σ = 1.0,
    dist = default_distance
) = lj_potential_wca_fij(dist(r₁, r₂), ϵ = ϵ, σ = σ)

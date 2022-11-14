using Test
using LennardJones
using LinearAlgebra

@testset "Test" begin

@test LennardJones.lj_potential_uij(LennardJones.R_MIN) ≈ -1.0

r₁ = [0.916914, 0.111325, 0.594877]
r₂ = [0.454572, 0.979098, 0.988925]
dr = r₂ - r₁
r = norm(dr)
ϵ = 0.8
σ = 0.9
dist = (a, b) -> a - b

u = LennardJones.lj_potential_uij(r, ϵ = ϵ, σ = σ)
@test u ≈ -0.7509440391688087
@test LennardJones.lj_potential_uij(dr, ϵ = ϵ, σ = σ) ≈ u
@test LennardJones.lj_potential_uij(r₁, r₂, ϵ = ϵ, σ = σ, dist = dist) ≈ u

f = LennardJones.lj_potential_fij(dr, ϵ = ϵ, σ = σ)
@test f ≈ [0.7369705345615358, -1.3832252568186918, -0.6281102846008022]
@test LennardJones.lj_potential_fij(r₁, r₂, ϵ = ϵ, σ = σ, dist = dist) ≈ -f

cutoff = 2.5
u₂ = LennardJones.lj_potential_uij_cutoff(dr, ϵ = ϵ, σ = σ, r_cutoff = cutoff)
@test u₂ ≈ -0.7439934985138913
@test LennardJones.lj_potential_uij_cutoff(r₁, r₂, ϵ = ϵ, σ = σ, r_cutoff = cutoff, dist = dist) ≈ u₂

f₂ = LennardJones.lj_potential_fij_cutoff(dr, ϵ = ϵ, σ = σ, r_cutoff = cutoff)
@test f₂ ≈ f
@test LennardJones.lj_potential_fij_cutoff(r₁, r₂, ϵ = ϵ, σ = σ, r_cutoff = cutoff, dist = dist) ≈ -f₂

r₃ = [0.669256, 0.300412, 0.671123]
dr₂ = r₃ - r₁
u₃ = LennardJones.lj_potential_wca_uij(dr₂, ϵ = ϵ, σ = σ)
@test u₃ ≈ 759679.8217248165
@test LennardJones.lj_potential_wca_uij(r₁, r₃, ϵ = ϵ, σ = σ, dist = dist) ≈ u₃

f₃ = LennardJones.lj_potential_wca_fij(dr₂, ϵ = ϵ, σ = σ)
@test f₃ ≈ [-2.1962740885915317e7, 1.6768563042159226e7, 6.761627492701628e6]
@test LennardJones.lj_potential_wca_fij(r₁, r₃, ϵ = ϵ, σ = σ, dist = dist) ≈ -f₃

end

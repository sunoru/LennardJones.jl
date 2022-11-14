# LennardJones.jl

[![CI](https://github.com/sunoru/LennardJones.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sunoru/LennardJones.jl/actions/workflows/ci.yml)

A simple package of functions for Lennard-Jones potential energies.

## Usage

`LennardJones` provides `lj_potential_uij` and `lj_potential_fij`, and their cutoff versions `lj_potential_uij_cutoff` and `lj_potential_fij_cutoff`, for calculating the potential energy and force between two particles.

A repulsive potential called Weeks–Chandler–Anderson potential is also implemented here. See `lj_potential_wca_uij` and `lj_potential_wca_fij`.

You can use keyword arguments to set the LJ parameters. The default values are `ϵ = 1.0`, `σ = 1.0`, `r_cutoff = 2.5σ`,
and a distance function `dist(r₁, r₂)` can be passed to specify the way to calculate the distance between two particles.

## License

The [MIT License](https://sunoru.mit-license.org/).

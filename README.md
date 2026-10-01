# MagBerry.jl

`MagBerry.jl` is a Julia package built on [`Sunny.jl`](https://github.com/SunnySuite/Sunny.jl) for topological magnonics, thermal transport, and related bosonic excitations.
It provides practical workflows for computing Berry-curvature-based responses, Boltzmann transport, spin and orbital effects, and several visualization tools for magnon, phonon, and their hybrid models.

It currently focuses on:

- magnon Berry curvature, Chern numbers
- intrinsic thermal Hall effect and symmetry analysis tools
- spin Berry curvature, spin Nernst effect
- magnon orbital moment Berry curvature, orbital Nernst effect
- Boltzmann transport, spin polarization, and nonlinear thermal Hall / extended BCD calculations
- group velocity, chirality, and edge / surface state calculations
- phonon models, magnon-phonon hybrid models, hybrid Berry curvature, and hybrid thermal Hall effect

## Installation

For local development:

```julia
using Pkg
Pkg.activate("path/to/MagBerry")
```

Then load the package with:

```julia
using MagBerry
```

Recommended environment:

- Julia `1.10`

## Quick Start

```julia
using MagBerry, Sunny

units = Units(:meV, :angstrom)
a=b=c=5
latvecs = lattice_vectors(a, b, c, 90, 90, 90)   
positions = [[0.0, 0.0, 0.0]]    
cryst = Crystal(latvecs, positions, 136; types=["Mn"])   
moments = [1 => Moment(s=5/2, g=2)]    
sys = System(cryst, moments, :dipole)
J1, J2, J3, d, D= 1.0, -0.3, -0.2, 0.1,0.1
set_exchange!(sys, J1 * I + d * dmvec([1 / sqrt(2), -1 / sqrt(2), 0]), Bond(2, 1, [0, 0, 0]))
set_exchange!(sys, J2, Bond(1, 1, [0, 0, 1]))
set_exchange!(sys, J2, Bond(1, 1, [1, 0, 0]))
set_exchange!(sys, J3 + D, Bond(1, 1, [1, 1, 0]))
set_exchange!(sys, J3 - D, Bond(1, 1, [1, -1, 0]))
θ = π / 3
ϕ = π / 3
set_dipole!(sys, [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)], (1, 1, 1, 1))
set_dipole!(sys, [-sin(θ) * cos(ϕ), -sin(θ) * sin(ϕ), -cos(θ)], (1, 1, 1, 2))
swt = SpinWaveTheory(sys; measure=nothing)
grid = q_space_grid(
    cryst,
    [1, 0, 0], range(-0.5, 0.5, 101),
    [0, 1, 0], range(-0.5, 0.5, 101),
    [0, 0, 1], range(-0.5, 0.5, 101)
)
bc = magnon_berry_curvature(swt, grid)
kyz = thermal_hall_conductivity(bc; temperature=20.0, component=:yz, units=units)
```

## Notes

`MagBerry.jl` relies on `Sunny.jl` for crystal handling, magnetic structures, and spin-wave calculations. Please also refer to the `Sunny.jl` documentation for the underlying modeling workflow.

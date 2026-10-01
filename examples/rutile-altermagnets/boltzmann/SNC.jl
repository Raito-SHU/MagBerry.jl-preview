using Pkg
Pkg.activate("d:/Trae/MapBerry/MagBerry")
using MagBerry,Sunny
using GLMakie
using LinearAlgebra


include("../system.jl")

swt = SpinWaveTheory(sys; measure=nothing)
n_neel = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]

# 1. Define 3D Grid for Spin Berry Curvature
# 50x50x50 grid
nx,ny,nz=50,50,50
grid = q_space_grid(cryst, 
    [1, 0, 0], range(-0.5, 0.5, nx),
    [0, 1, 0], range(-0.5, 0.5, ny),
    [0, 0, 1], range(-0.5, 0.5, nz)
)


# 2. Compute Spin Berry Curvature (projected on Néel axis)
# Note: spin_magnon_berry_curvature calculates the tensor for all bands.
# It doesn't take an axis argument for projection during calculation in the current API?
# Checking API: spin_magnon_berry_curvature(swt, grid; ...) returns rank-4 tensor (band, spin_comp, dim1, dim2).
# Wait, I previously used `axis=n_neel`. Let me check if that is supported.
# Reading source code: spin_magnon_berry_curvature does NOT take `axis`. 
# It computes the full tensor. We need to project it later or use specific components.
# The `spin_nernst_conductivity` function takes `component` which can specify spin axis.
println("Computing 3D Spin Berry Curvature...")
sbc = spin_magnon_berry_curvature(swt, grid)

# 3. Calculate Spin Nernst Conductivity vs Temperature
temps = range(0.1, 100, length=50)

# User requested "xzy" and "xxz" components.
# Interpreting as (Spin polarization, Current direction, Gradient direction)
# or standard indices. 
# Based on `_resolve_spin_nernst_component`:
# Tuple: (spin, current, gradient)
# "xzy" -> Spin X, Current Z, Gradient Y? Or indices 1, 3, 2?
# "xxz" -> Spin X, Current X, Gradient Z? Or indices 1, 1, 3?
# Since user previously asked for Néel axis quantization, maybe they want:
# Spin along Néel, Current X, Gradient Y?
# But here they explicitly said "xzy, xxz". I will assume (Spin, Current, Gradient) axes.

# However, the user might mean tensor indices indices $\alpha_{ij}^k$.
# Common notation $\alpha_{ij}^k = J_i^k / \nabla_j T$. (Current i, Gradient j, Spin k).
# The library uses (spin, current, gradient).
# So I will calculate:
# 1. Spin X, Current Z, Gradient Y (xzy)
# 2. Spin X, Current X, Gradient Z (xxz)

println("Calculating Spin Nernst Conductivity...")

# Component 1: Spin X, Current Z, Gradient Y
# We pass component as a tuple of symbols (:x, :z, :y)
axy_1 = [spin_nernst_conductivity(sbc; temperature=T, component=(:x, :z, :y)) for T in temps]

# Component 2: Spin X, Current X, Gradient Z
# We pass component as a tuple of symbols (:x, :x, :z)
axy_2 = [spin_nernst_conductivity(sbc; temperature=T, component=(:x, :x, :z)) for T in temps]

# 4. Plot
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1], xlabel="Temperature (K)", ylabel="α (J/(K·m))", title="Spin Nernst Conductivity")

lines!(ax, temps, axy_1, label="α_xzy (Spin X, J_z, ∇_y T)")
lines!(ax, temps, axy_2, label="α_xxz (Spin X, J_x, ∇_z T)")

axislegend(ax)
save("SNC_components.png", fig)
println("Saved SNC_components.png")

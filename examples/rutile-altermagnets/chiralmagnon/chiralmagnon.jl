using MagBerry
using Sunny
using GLMakie
using LinearAlgebra

# Include system definition
include("system.jl")

# 1. Determine Néel vector (Quantization Axis)
# Based on system.jl:
# Site 1: [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
# Site 2: [-sin(θ)*cos(ϕ), -sin(θ)*sin(ϕ), -cos(θ)]
# Néel vector is along (θ, ϕ).
n_neel = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
println("Néel vector: ", n_neel)

# 2. Setup SpinWaveTheory
# measure=nothing speeds up initialization if we don't need intensities immediately
swt = SpinWaveTheory(sys; measure=nothing)

# 3. Define k-path
# Points from user request: R2, G, U2, S2, G, S1, U1, G, R1
# R2=[-0.5,-0.5,0.5], G=[0,0,0], U2=[0,-0.5,0], S2=[0.5,-0.5,0]
# G=[0,0,0], S1=[0.5,0.5,0], U1=[0.5,0,0], G=[0,0,0], R1=[0.5,0.5,0.5]
points = [
    [0.5, -0.5, 0.5], # R2
    [0.0, 0.0, 0.0],   # G
    [0.0, -0.5, 0.0],  # U2
    [0.5, -0.5, 0.0],  # S2
    [0.0, 0.0, 0.0],   # G
    [0.5, 0.5, 0.0],   # S1
    [0.5, 0.0, 0.0],   # U1
    [0.0, 0.0, 0.0],   # G
    [0.5, 0.5, 0.5]    # R1
]
labels = ["R2", "Γ", "U2", "S2", "Γ", "S1", "U1", "Γ", "R1"]
path = q_space_path(cryst, points, 4000; labels=labels)

# 4. Calculate Chirality
# axis=n_neel projects the chirality onto the Néel order
res = intensities_chirality(swt, path; axis=n_neel)

# 5. Plot
fig = plot_intensities(res; title="Chiral Magnon",allpositive=false)
save("chiral_magnon.png", fig)
println("Saved chiral_magnon.png")

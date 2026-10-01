using MagBerry
using Sunny
using GLMakie
using LinearAlgebra

include("../system.jl")

# 1. Setup
θ,ϕ=0,0
set_dipole!(sys, [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)], (1, 1, 1, 1))
set_dipole!(sys, [-sin(θ) * cos(ϕ), -sin(θ) * sin(ϕ), -cos(θ)], (1, 1, 1, 2))
swt = SpinWaveTheory(sys; measure=nothing)
n_neel = [sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]

# 2. Define Integration Grid (3D)
# 50x50x50 grid for accurate 3D integration
nx,ny,nz=50,50,50
grid = q_space_grid(cryst, 
    [1, 0, 0], range(-0.5, 0.5, nx),
    [0, 1, 0], range(-0.5, 0.5, ny),
    [0, 0, 1], range(-0.5, 0.5, nz)
)


# 3. Parameters
# Temperature range 0.001 to 50.001 K
temps = range(0.1, 100, length=50)
tau = 5e-13# Relaxation time (s)

# 4. Calculate Conductivity vs Temperature
# We calculate for each temperature and extract xx, xy, zz components for spin up/down/total
kxx_up, kxy_up, kzz_up = Float64[], Float64[], Float64[]
kxx_dn, kxy_dn, kzz_dn = Float64[], Float64[], Float64[]
ls_xx, ls_xy, ls_zz = Float64[], Float64[], Float64[]

println("Calculating Boltzmann conductivity over temperature range...")

for t in temps
    # spin_axis=n_neel allows separating spin-up/down contributions
    res = magnon_boltzmann_conductivity(swt, grid; temperature=t, relaxation_time=tau, periodic=true, spin_axis=n_neel)
    
    # Extract spin-up components
    push!(kxx_up, res.spin_up[1,1])
    push!(kxy_up, res.spin_up[1,2])
    push!(kzz_up, res.spin_up[3,3])
    
    # Extract spin-down components
    push!(kxx_dn, res.spin_down[1,1])
    push!(kxy_dn, res.spin_down[1,2])
    push!(kzz_dn, res.spin_down[3,3])

    # Extract Spin Seebeck components
    push!(ls_xx, res.spin_seebeck[1,1])
    push!(ls_xy, res.spin_seebeck[1,2])
    push!(ls_zz, res.spin_seebeck[3,3])
end

# 5. Plot
fig = Figure(resolution = (1000, 800))

# Plot xx component
ax1 = Axis(fig[1, 1], xlabel="Temperature (K)", ylabel="κ_xx (arb. units)", title="Longitudinal (xx)")
lines!(ax1, temps, kxx_up, label="Spin Up")
lines!(ax1, temps, kxx_dn, label="Spin Down", linestyle=:dash)
axislegend(ax1)

# Plot xy component
ax2 = Axis(fig[1, 2], xlabel="Temperature (K)", ylabel="κ_xy (arb. units)", title="Transverse (xy)")
lines!(ax2, temps, kxy_up, label="Spin Up")
lines!(ax2, temps, kxy_dn, label="Spin Down", linestyle=:dash)
axislegend(ax2)

# Plot zz component
ax3 = Axis(fig[2, 1], xlabel="Temperature (K)", ylabel="κ_zz (arb. units)", title="Longitudinal (zz)")
lines!(ax3, temps, kzz_up, label="Spin Up")
lines!(ax3, temps, kzz_dn, label="Spin Down", linestyle=:dash)
axislegend(ax3)

save("boltzmann_spin_components.png", fig)
println("Saved boltzmann_spin_components.png")

# Plot Spin Seebeck
fig_sse = Figure(resolution = (800, 600))
ax_sse = Axis(fig_sse[1, 1], xlabel="Temperature (K)", ylabel="L_sT (arb. units)", title="Spin Seebeck Coefficient")
lines!(ax_sse, temps, ls_xx, label="xx", color=:red)
lines!(ax_sse, temps, ls_xy, label="xy", color=:blue)
lines!(ax_sse, temps, ls_zz, label="zz", color=:green)
axislegend(ax_sse)
save("boltzmann_spin_seebeck.png", fig_sse)
println("Saved boltzmann_spin_seebeck.png")

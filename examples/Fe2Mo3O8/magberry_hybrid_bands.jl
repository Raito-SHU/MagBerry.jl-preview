using Pkg
Pkg.activate("d:/Trae/MapBerry/MagBerry")
using MagBerry,Sunny
using GLMakie
using LinearAlgebra
using Statistics

const units = Units(:meV, :angstrom)
const a = 5.773
const c = 20.0
const S = 2.0
const Dmp = 0.8

latvecs = lattice_vectors(a, a, c, 90, 90, 120)
positions = [
    [2 / 3, 1 / 3, 1 / 2],
    [1 / 3, 2 / 3, 1 / 2],
]

cryst = Crystal(latvecs, positions, 1; types=["Fe1", "Fe2"])
sys = System(cryst, [1 => Moment(s=S, g=2), 2 => Moment(s=S, g=2)], :dipole)

J1 = 0.5742*0.89
J21 = -0.06522*0.89
J22 = -0.01386*0.89
J3 = -0.2113*0.89
Δ1 = -3.005*0.89
Δ2 = -2.250*0.89
V1 = 36.5
V21 = 9.1
V22 = 7.6
V3 = 15.5
V4 = 5.5

nn_bonds = [
    Bond(1, 2, [0, 0, 0]),
    Bond(1, 2, [0, -1, 0]),
    Bond(1, 2, [1, 0, 0]),
]

nna_bonds = [
    Bond(1, 1, [1, 0, 0]),
    Bond(1, 1, [0, 1, 0]),
    Bond(1, 1, [1, 1, 0]),
]

nnb_bonds = [
    Bond(2, 2, [1, 0, 0]),
    Bond(2, 2, [0, 1, 0]),
    Bond(2, 2, [1, 1, 0]),
]

third_bonds = [
    Bond(1, 2, [1, 1, 0]),
    Bond(2, 1, [1, 1, 0]),
    Bond(1, 2, [1, -1, 0]),
]

fourth_bonds = [
    Bond(1, 2, [2, 1, 0]),
    Bond(2, 1, [1, 2, 0]),
    Bond(1, 2, [2, 0, 0]),
    Bond(1, 2, [0, -2, 0]),
    Bond(1, 2, [0, 1, 0]),
    Bond(1, 2, [-1, 0, 0]),
]

for bond in nn_bonds
    set_exchange!(sys, J1, bond)
end

for bond in nna_bonds
    set_exchange!(sys, J21, bond)
end

for bond in nnb_bonds
    set_exchange!(sys, J22, bond)
end

for bond in third_bonds
    set_exchange!(sys, J3, bond)
end

set_onsite_coupling!(sys, s -> Δ1 * s[3]^2, 1)
set_onsite_coupling!(sys, s -> Δ2 * s[3]^2, 2)

set_dipole!(sys, [0.0, 0.0, -2.0], (1, 1, 1, 1))
set_dipole!(sys, [0.0, 0.0, 2.0], (1, 1, 1, 2))

psys = PhononSystem(cryst, [1.0, 1.0])

function bond_vector(cryst::Crystal, bond::Bond)
    ri = cryst.positions[bond.i]
    rj = cryst.positions[bond.j] + bond.n
    return cryst.latvecs * (rj - ri)
end

function spring_matrix(k::Real, bond::Bond)
    r = bond_vector(cryst, bond)
    r̂ = r / norm(r)
    Φ = zeros(Float64, 3, 3)
    Φ .-= k .* (r̂ * r̂')
    return Φ
end

for bond in nn_bonds
    set_force!(psys, spring_matrix(V1, bond), bond)
end

for bond in nna_bonds
    set_force!(psys, spring_matrix(V21, bond), bond)
end

for bond in nnb_bonds
    set_force!(psys, spring_matrix(V22, bond), bond)
end

for bond in third_bonds
    set_force!(psys, spring_matrix(V3, bond), bond)
end

for bond in fourth_bonds
    set_force!(psys, spring_matrix(V4, bond), bond)
end

enforce_asr!(psys)

dmi_term = DMPhononCoupling(:D, Dmp, nn_bonds; directions=1:2)
hsys = HybridSpinWaveTheory(sys, psys; measure=nothing, dmi_terms=[dmi_term])

kpoints = [
    [0.0, 0.0, 0.0],
    [1 / 2, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [-1 / 3, 2 / 3, 0.0],
    [1 / 3, 1 / 3, 0.0],
    [0.0, 0.0, 0.0],
]

labels = ["Γ", "M", "Γ", "K₁", "K₂", "Γ"]
kpath = q_space_path(cryst, kpoints, 4000; labels=labels)
projection_data = hybrid_band_projections(hsys, kpath.qs; axis=(0.0, 0.0, 1.0))
bands = projection_data.energies
selected_bands = collect(3:6)
berry_plot_bands = [5, 6]
x = 1:length(kpath.qs)
phonon_threshold = 0.5

fig = Figure(size=(1300, 600))
ax1 = Axis(
    fig[1, 1];
    xlabel="Momentum",
    ylabel="Energy (meV)",
    title="Magnon-Phonon Character",
    xticks=kpath.xticks,
)

for band in selected_bands
    lines!(ax1, x, bands[band, :], color=(:black, 0.25), linewidth=1.0)
    scatter!(ax1, x, bands[band, :];
        color=projection_data.magnon_phonon_character[band, :],
        colorrange=(-1.0, 1.0),
        colormap=:cool,
        markersize=7,
    )
end

Colorbar(fig[1, 2], limits=(-1.0, 1.0), colormap=:cool, label="+1 Magnon   0 Hybrid   -1 Phonon")

ax2 = Axis(
    fig[1, 3];
    xlabel="Momentum",
    ylabel="Energy (meV)",
    title="Phonon Circular Polarization",
    xticks=kpath.xticks,
)

for band in selected_bands
    lines!(ax2, x, bands[band, :], color=(:black, 0.25), linewidth=1.0)
    phonon_like = findall(projection_data.phonon_weight[band, :] .>= phonon_threshold)
    scatter!(ax2, x[phonon_like], bands[band, phonon_like];
        color=projection_data.phonon_circular_polarization[band, phonon_like],
        colorrange=(-1.0, 1.0),
        colormap=:cool,
        markersize=7,
    )
end

Colorbar(fig[1, 4], limits=(-1.0, 1.0), colormap=:cool, label="Circular Polarization (phonon-like only)")

save("magnon-polaron-bands-magberry.png", fig)
display(fig)

berry_grid = q_space_grid(
    cryst,
    [1, 0, 0], range(-0.5, 0.5, 301),
    [0, 1, 0], range(-0.5, 0.5, 301),
)
kubo_fdstep = 1e-5
kubo_degeneracy_tol = 1e-5

berry_data_fukui = hybrid_berry_curvature(hsys, berry_grid; bands=selected_bands, method=:fukui)
chern_numbers_fukui = hybrid_chern_numbers(berry_data_fukui)

berry_data_kubo = hybrid_berry_curvature(
    hsys,
    berry_grid;
    bands=selected_bands,
    method=:kubo,
    fdstep=kubo_fdstep,
    degeneracy_tol=kubo_degeneracy_tol,
)
chern_numbers_kubo = hybrid_chern_numbers(berry_data_kubo)

println("Fukui Chern numbers:")
for (idx, band) in enumerate(selected_bands)
    println("hybrid band ", band, " chern = ", chern_numbers_fukui[idx])
end

println("\nKubo Chern numbers:")
for (idx, band) in enumerate(selected_bands)
    println("hybrid band ", band, " chern = ", chern_numbers_kubo[idx])
end

v1 = cryst.recipvecs * berry_grid.axes[1]
v2 = cryst.recipvecs * berry_grid.axes[2]
e1 = normalize(v1)
v2_orth = v2 - dot(v2, e1) * e1
e2 = normalize(v2_orth)
X = zeros(Float64, size(berry_grid.qs))
Y = zeros(Float64, size(berry_grid.qs))
for idx in CartesianIndices(berry_grid.qs)
    q_cart = cryst.recipvecs * berry_grid.qs[idx]
    X[idx] = dot(q_cart, e1)
    Y[idx] = dot(q_cart, e2)
end
b1_2d = [dot(v1, e1), dot(v1, e2)]
b2_2d = [dot(v2, e1), dot(v2, e2)]
bzcell_2d = Sunny.Brillouin.cartesianize!(Sunny.Brillouin.wignerseitz(eachcol(hcat(b1_2d, b2_2d))))
segments = GLMakie.Point2f[]
for face in bzcell_2d
    append!(segments, face)
    push!(segments, face[1])
    push!(segments, GLMakie.Point2f(NaN))
end

function plot_hybrid_berry_curvature_maps(data, filename; title_prefix="")
    plot_band_indices = [findfirst(==(band), data.bands) for band in berry_plot_bands]
    omega_samples = Float64[]
    for band_idx in plot_band_indices
        band_slice = abs.(@view data.curvature[band_idx, 3, :, :])
        append!(omega_samples, vec(band_slice))
    end
    filter!(isfinite, omega_samples)
    omega_max = if isempty(omega_samples)
        1.0
    else
        clipped = quantile(omega_samples, 0.995)
        clipped > 0 ? clipped : maximum(omega_samples)
    end
    omega_max = omega_max > 0 ? omega_max : 1.0
    levels = collect(range(-omega_max, omega_max; length=101))
    ncols = min(3, length(berry_plot_bands))
    nrows = cld(length(berry_plot_bands), ncols)
    berry_fig = Figure(size=(420 * ncols, 340 * nrows))
    last_plot = nothing

    for (plot_idx, band) in enumerate(berry_plot_bands)
        row = cld(plot_idx, ncols)
        col = mod1(plot_idx, ncols)
        ax = Axis(
            berry_fig[row, col];
            xlabel="kx (Å⁻¹)",
            ylabel="ky (Å⁻¹)",
            title=isempty(title_prefix) ? "Band $band" : "$title_prefix Band $band",
            aspect=DataAspect(),
        )
        band_idx = findfirst(==(band), data.bands)
        contour = contourf!(
            ax,
            X,
            Y,
            data.curvature[band_idx, 3, :, :];
            colormap=:balance,
            levels=levels,
        )
        lines!(ax, segments; color=:black, linewidth=1.5, inspectable=false)
        last_plot = contour
    end

    Colorbar(berry_fig[1:nrows, ncols + 1], last_plot, label="Ωz")
    save(filename, berry_fig)
    display(berry_fig)
    return berry_fig
end

plot_hybrid_berry_curvature_maps(
    berry_data_fukui,
    "magnon-polaron-berry-magberry-fukui.png";
    title_prefix="Fukui",
)

plot_hybrid_berry_curvature_maps(
    berry_data_kubo,
    "magnon-polaron-berry-magberry-kubo.png";
    title_prefix="Kubo",
)

println("\nKubo settings:")
println("fdstep = ", kubo_fdstep)
println("degeneracy_tol = ", kubo_degeneracy_tol)

thermal_hall_temperatures = collect(range(0.1, 300.0; length=50))
kappa_xy_fukui = [
    thermal_hall_conductivity(berry_data_fukui; temperature=T, component=:xy, units=units)
    for T in thermal_hall_temperatures
]
kappa_xy_kubo = [
    thermal_hall_conductivity(berry_data_kubo; temperature=T, component=:xy, units=units)
    for T in thermal_hall_temperatures
]

thermal_hall_fig = Figure(size=(760, 480))
thermal_hall_ax = Axis(
    thermal_hall_fig[1, 1];
    xlabel="Temperature (K)",
    ylabel="κxy (W/K)",
    title="Hybrid Thermal Hall Conductivity",
)
lines!(thermal_hall_ax, thermal_hall_temperatures, kappa_xy_fukui; linewidth=3, label="Fukui")
lines!(thermal_hall_ax, thermal_hall_temperatures, kappa_xy_kubo; linewidth=3, linestyle=:dash, label="Kubo")
axislegend(thermal_hall_ax; position=:rb)
save("magnon-polaron-thermal-hall-magberry.png", thermal_hall_fig)
display(thermal_hall_fig)

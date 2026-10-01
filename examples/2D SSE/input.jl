using Pkg
Pkg.activate("d:/Trae/MapBerry/MagBerry")
using Sunny, MagBerry, LinearAlgebra, GLMakie, DelimitedFiles
units = Units(:meV, :angstrom)

# CTO model
latvecs = lattice_vectors(4.12, 4.12, 20, 90.000000, 90.000000, 90.000000)
positions = [[0.0, 0.0, 0.5], [0.5, 0.5, 0.5]]
cryst = Crystal(latvecs, positions, 1; types=["Cr1", "Cr2"])
moments = [1 => Moment(s=3/2, g=2), 2 => Moment(s=3/2, g=2)]
sys = System(cryst, moments, :dipole)

J1=3.90
J2a=-7.9
J2b=-1.21
K=-0.73

set_exchange!(sys, J1, Bond(2,1,[1,1,0]))
set_exchange!(sys, J1, Bond(1,2,[-1,0,0]))
set_exchange!(sys, J1, Bond(1,2,[0,-1,0]))
set_exchange!(sys, J1, Bond(1,2,[0,0,0]))

set_exchange!(sys, J2a, Bond(1,1,[1,0,0]))
set_exchange!(sys, J2b, Bond(1,1,[0,1,0]))

set_exchange!(sys, J2a, Bond(2,2,[0,1,0]))
set_exchange!(sys, J2b, Bond(2,2,[1,0,0]))

set_onsite_coupling!(sys, S -> K*S[3]^2, 1)
set_onsite_coupling!(sys, S -> K*S[3]^2, 2)

randomize_spins!(sys)
minimize_energy!(sys)

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))

kpoints = [[0.0, 0.0, 0.0],[0.5, 0.5, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.5], [0.0, 0.0, 0.0]]
kpath = q_space_path(cryst, kpoints, 5000)
chi = intensities_chirality(swt,kpath;axis=Vec3((0, 0, 1)))

fig = plot_intensities(chi; units, allpositive=false)  
  
# Save the figure  
save("chiralmagnon.png", fig)

grid_2d = q_space_grid(cryst,
    [1, 0, 0], range(-0.5, 0.5, length=81),
    [0, 1, 0], range(-0.5, 0.5, length=81)
)

temps = range(0.1, 100.0, length=60)
tau = 1e-11

results_sse = [
    magnon_boltzmann_conductivity(
        swt,
        grid_2d;
        temperature=T,
        relaxation_time=tau,
        energy_unit=:meV,
    )
    for T in temps
]

fig_sse, ax_sse, values_sse = plot_boltzmann_conductivity(
    first(results_sse).spin_seebeck;
    temperatures=temps,
    results=results_sse,
    data_kind=:spin_seebeck,
    spin_axis=[0, 0, 1],
    components=(:xx, :xy, :yy),
    linewidth=3,
)

ax_sse.title = "2D Spin Seebeck"

save("2d_spin_seebeck.png", fig_sse)
writedlm("2d_spin_seebeck.csv", hcat(collect(temps), permutedims(values_sse)), ',')

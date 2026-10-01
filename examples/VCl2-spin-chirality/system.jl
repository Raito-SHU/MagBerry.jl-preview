using Sunny, MagBerry

# 本文件只负责 VCl2 的系统设置：晶体、磁模型、声子模型和杂化系统。
# 具体的能带、Berry 曲率、热霍尔和自适应网格计算放在 input.jl。

units = Units(:meV, :angstrom)

poscar_unitcell_path = joinpath(@__DIR__, "POSCAR-unitcell")
poscar_supercell_path = joinpath(@__DIR__, "POSCAR")
force_constants_path = joinpath(@__DIR__, "new-FORCE_CONSTANTS")

S = 3 / 2
g = 2.0
J = [5.955, -0.004, 0.026] ./ S^2
dJdr = 0.5

cryst = read_poscar_crystal(poscar_unitcell_path)
mag_cryst = subcrystal(cryst, "V")

sys = System(mag_cryst, [site => Moment(s=S, g=g) for site in 1:Sunny.natoms(mag_cryst)], :dipole)

j1_bonds = Bond[
    Bond(2, 3, [0, 1, 0]),
    Bond(1, 2, [0, 0, 0]),
    Bond(3, 1, [1, 1, 0]),
    Bond(2, 1, [0, 1, 0]),
    Bond(3, 1, [1, 0, 0]),
    Bond(2, 3, [0, 0, 0]),
    Bond(3, 2, [1, 0, 0]),
    Bond(1, 3, [0, 0, 0]),
    Bond(2, 1, [1, 1, 0]),
]

j2_bonds = Bond[
    Bond(2, 2, [0, 1, 0]),
    Bond(1, 1, [0, 1, 0]),
    Bond(3, 3, [0, 1, 0]),
    Bond(2, 2, [1, 1, 0]),
    Bond(3, 3, [1, 1, 0]),
    Bond(1, 1, [1, 1, 0]),
    Bond(3, 3, [1, 0, 0]),
    Bond(1, 1, [1, 0, 0]),
    Bond(2, 2, [1, 0, 0]),
]

j3_bonds = Bond[
    Bond(1, 3, [0, 1, 0]),
    Bond(2, 1, [1, 2, 0]),
    Bond(3, 2, [1, 1, 0]),
    Bond(3, 2, [1, -1, 0]),
    Bond(2, 1, [1, 0, 0]),
    Bond(3, 1, [0, 1, 0]),
    Bond(1, 2, [1, 0, 0]),
    Bond(3, 1, [2, 1, 0]),
    Bond(2, 3, [1, 1, 0]),
]

for (value, bonds) in zip(J, [j1_bonds, j2_bonds, j3_bonds])
    for bond in bonds
        set_exchange!(sys, value, bond)
    end
end

# 120 度顺时针自旋手性。若要比较相反手性，把 angle 的符号改成正号即可。
for site in 1:Sunny.natoms(mag_cryst)
    angle = -2π * (site - 1) / 3
    set_dipole!(sys, S .* [cos(angle), sin(angle), 0.0], (1, 1, 1, site))
end

set_onsite_coupling!(sys, S -> 0.0001 * S[3]^2, 1)
set_onsite_coupling!(sys, S -> 0.0001 * S[3]^2, 2)
set_onsite_coupling!(sys, S -> 0.0001 * S[3]^2, 3)
set_field!(sys, [0, 0, 0] * units.meV)

minimize_energy!(sys; maxiters=10000)

mag_moments = [magnetic_moment(sys, site) for site in eachsite(sys)]
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))

psys = read_phonopy_phonon(cryst, poscar_supercell_path, force_constants_path)
# 如需强制声子和力常数满足对称性，可按需要打开下面三行。
# symmetrize_phonon!(psys; symmetry=:spacegroup)
# enforce_permutation_symmetry!(psys)
# enforce_asr!(psys)

ex_term = EXPhononCoupling(:J1, dJdr, j1_bonds; directions=1:3)

# hsys 是带 magnon-phonon exchange-displacement 耦合的主系统。
hsys = HybridSpinWaveTheory(swt, psys; ex_terms=[ex_term])

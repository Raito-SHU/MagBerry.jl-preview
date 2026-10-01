using Pkg
Pkg.activate("d:/Trae/MapBerry/MagBerry")
using MagBerry,Sunny
using GLMakie
using LinearAlgebra

#RuO2(MnF2) 磁结构
units = Units(:meV, :angstrom)
a = 5    #without loss of generality, we consider a cubic unit cell with a lattice constant of a = 0.5 nm
b = 5
c = 5
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



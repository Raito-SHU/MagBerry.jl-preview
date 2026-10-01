using Pkg
Pkg.activate("d:/Trae/MapBerry/MagBerry")
using MagBerry,Sunny
using GLMakie
using LinearAlgebra
using Printf

include("../system.jl")
θ = 0
ϕ = π / 2
set_dipole!(sys, [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)], (1, 1, 1, 1))
set_dipole!(sys, [-sin(θ) * cos(ϕ), -sin(θ) * sin(ϕ), -cos(θ)], (1, 1, 1, 2))


swt = SpinWaveTheory(sys; measure=nothing)

# 1. Define 3D Grid for Berry Curvature
# 50x50x50 grid for accurate 3D integration
nx,ny,nz=200,200,200
grid = q_space_grid(cryst, 
    [1, 0, 0], range(-0.5, 0.5, nx),
    [0, 1, 0], range(-0.5, 0.5, ny),
    [0, 0, 1], range(-0.5, 0.5, nz)
)


# 2. Compute Berry Curvature (3D)
println("Computing 3D Berry Curvature...")
bc = magnon_berry_curvature(swt, grid)

# 3. Identify magnetic symmetry and symbolic thermal Hall tensor constraints
println("Identifying magnetic symmetry...")
magnetic_analysis = identify_magnetic_symmetry(swt)
tensor_report_symbolic = thermal_hall_tensor_report(swt)

println("\nMagnetic symmetry summary:")
display(magnetic_analysis)

println("\nSymbolic thermal Hall tensor report:")
display(tensor_report_symbolic)

open("THC_symmetry_report.txt", "w") do io
    println(io, "Magnetic symmetry operations")
    println(io, "Parent crystal: ", magnetic_analysis.crystal.sg.label)
    println(io, "Supercell dims: ", magnetic_analysis.supercell_dims)
    println(io, "Number of magnetic operations: ", length(magnetic_analysis.magnetic_symops))
    println(io, "Number of unitary operations: ", length(magnetic_analysis.unitary_symops))
    println(io, "Number of antiunitary operations: ", length(magnetic_analysis.antiunitary_symops))
    println(io)
    println(io, "Operations:")
    for op in magnetic_analysis.operation_strings
        println(io, "  ", op)
    end
    println(io)
    println(io, "Independent tensor elements: ", isempty(tensor_report_symbolic.independent_tensor_elements) ? "none" : join(tensor_report_symbolic.independent_tensor_elements, ", "))
    println(io, "Zero tensor elements: ", isempty(tensor_report_symbolic.zero_tensor_elements) ? "none" : join(tensor_report_symbolic.zero_tensor_elements, ", "))
    println(io)
    println(io, "Component relations:")
    for rel in tensor_report_symbolic.component_relations
        println(io, "  ", rel)
    end
    println(io)
    println(io, "Tensor relations:")
    for rel in tensor_report_symbolic.tensor_relations
        println(io, "  ", rel)
    end
    println(io)
    println(io, "Symbolic tensor:")
    for i in 1:3
        println(io, join(tensor_report_symbolic.symbolic_tensor[i, :], "\t"))
    end
end
println("Saved THC_symmetry_report.txt")

# 4. Calculate Thermal Hall Conductivity vs Temperature
# Range 0 to 100K. Avoid exactly 0 to prevent division by zero in Bose function if not handled.
# Implementation usually handles T->0, but let's use small epsilon.
temps = range(0.1, 100, length=50)

println("Calculating Thermal Hall Conductivity...")
# Calculate xy, yz, zx components
kxy = [thermal_hall_conductivity(bc; temperature=T, component=:xy) for T in temps]
kyz = [thermal_hall_conductivity(bc; temperature=T, component=:yz) for T in temps]
kzx = [thermal_hall_conductivity(bc; temperature=T, component=:zx) for T in temps]
tensor_reports = [thermal_hall_tensor_report(bc, swt; temperature=T) for T in temps]
tensors = [report.numeric_tensor for report in tensor_reports]

open("THC_tensor_vs_T.dat", "w") do io
    println(io, "# T(K) k_xx k_xy k_xz k_yx k_yy k_yz k_zx k_zy k_zz")
    for (T, tensor) in zip(temps, tensors)
        @printf(
            io,
            "%.6f %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
            T,
            tensor[1, 1], tensor[1, 2], tensor[1, 3],
            tensor[2, 1], tensor[2, 2], tensor[2, 3],
            tensor[3, 1], tensor[3, 2], tensor[3, 3],
        )
    end
end
println("Saved THC_tensor_vs_T.dat")

# 5. Plot
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1], xlabel="Temperature (K)", ylabel="κ (W/(K·m))", title="Thermal Hall Conductivity")

lines!(ax, temps, kxy, label="κ_xy")
lines!(ax, temps, kyz, label="κ_yz")
lines!(ax, temps, kzx, label="κ_zx")

axislegend(ax)
save("THC_components.png", fig)
println("Saved THC_components.png")

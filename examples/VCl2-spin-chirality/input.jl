using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using Sunny, MagBerry, GLMakie

include(joinpath(@__DIR__, "system.jl"))

# 本文件只放计算流程：
# 1. 杂化能带和 Berry 曲率投影图；
# 2. 普通粗网格 Berry 曲率和热霍尔；
# 3. Wannier90 式自适应网格热霍尔。

const HALL_TEMPERATURE = 30.0
const DEGENERACY_TOL = 1e-5
const SPREAD_TOL = DEGENERACY_TOL
const FDSTEP = 1e-5

# 自适应网格参数：当粗网格点的 |Omega_xy| 超过 threshold 时，
# 用该点周围 kmesh^D 个子网格点的平均值替换粗点贡献。
const ADAPTIVE_REFINEMENT = (kmesh=3, threshold=100.0)

kpoints = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1 / 3, 1 / 3, 0.0], [0.0, 0.0, 0.0]]
kpath = q_space_path(mag_cryst, kpoints, 600; labels=["G", "M", "K", "G"])

println("绘制 VCl2 杂化能带和 Berry 曲率投影图。")

plot_hybrid_bands(hsys, kpath;
    color=:character,
    filename=joinpath(@__DIR__, "vcl2_hybrid_bands.png"),
    title="VCl2 hybrid magnon-phonon bands",
)

plot_hybrid_bands(hsys, kpath;
    color=:phonon_chirality,
    axis=(0.0, 0.0, 1.0),
    filename=joinpath(@__DIR__, "vcl2_hybrid_phonon_chirality.png"),
    title="VCl2 hybrid phonon chirality",
)

plot_bc_bands(hsys, kpath;
    component=:z,
    projection=:magnon,
    degeneracy_handling=:hybrid,
    degeneracy_tol=1e-4,
    smearing=:adaptive,
    smearing_factor=0.5,
    smearing_max=0.05,
    filename=joinpath(@__DIR__, "vcl2_hybrid_bc_magnon_projection.png"),
    title="VCl2 hybrid Berry curvature, magnon projection",
)

plot_bc_bands(hsys, kpath;
    bands=bc_plot_bands,
    component=:z,
    projection=:phonon,
    degeneracy_handling=:hybrid,
    degeneracy_tol=1e-4,
    smearing=:adaptive,
    smearing_factor=0.5,
    smearing_max=0.05,
    filename=joinpath(@__DIR__, "vcl2_hybrid_bc_phonon_projection.png"),
    title="VCl2 hybrid Berry curvature, phonon projection",
)

plot_phonon(psys, kpath;
    color=:chirality,
    axis=(1.0, 0.0, 0.0),
    filename=joinpath(@__DIR__, "vcl2_phonon_chirality.png"),
    title="VCl2 phonon chirality",
)

plot_phonon(psys, kpath;
    color=:atom_projection,
    elements="V",
    filename=joinpath(@__DIR__, "vcl2_phonon_v_projection.png"),
    title="VCl2 V-site phonon projection",
)

println("构造热霍尔二维积分网格。")
hall_grid = q_space_grid(mag_cryst,
    [1, 0, 0], range(-0.5, 0.5; length=50),
    [0, 1, 0], range(-0.5, 0.5; length=50),
)
temperatures = collect(range(0.0001, 100.0; length=40))

println("计算 magnon / phonon / hybrid 的粗网格 Berry 曲率。")
magnon_bc = magnon_berry_curvature(swt, hall_grid;
    method=:kubo,
    fdstep=:auto,
    degeneracy_handling=:hybrid,
    degeneracy_tol=DEGENERACY_TOL,
)

hybrid_bc = hybrid_berry_curvature(hsys, hall_grid;
    method=:kubo,
    fdstep=FDSTEP,
    degeneracy_handling=:hybrid,
    degeneracy_tol=DEGENERACY_TOL,
    bloch_embedding=:basis_positions,
)

phonon_bc = phonon_berry_curvature(psys, hall_grid;
    method=:kubo,
    fdstep=FDSTEP,
    degeneracy_tol=DEGENERACY_TOL,
    bloch_embedding=:basis_positions,
)

println("绘制粗网格 Berry 曲率平面图。")
fig_mbc_plane = plot_berry_curvature_plane(magnon_bc;
    view=:subspace,
    sum_bands=true,
    component=:z,
    size=(720, 620),
)
save(joinpath(@__DIR__, "vcl2_magnon_bc_plane.png"), fig_mbc_plane)

fig_hbc_plane = plot_berry_curvature_plane(hybrid_bc;
    view=:subspace,
    sum_bands=true,
    component=:z,
    size=(720, 620),
)
save(joinpath(@__DIR__, "vcl2_hybrid_bc_plane.png"), fig_hbc_plane)

fig_pbc_plane = plot_berry_curvature_plane(phonon_bc;
    sum_bands=true,
    component=:z,
    size=(720, 620),
)
save(joinpath(@__DIR__, "vcl2_phonon_bc_plane.png"), fig_pbc_plane)

println("绘制粗网格热霍尔随温度曲线。")
fig_mthe = plot_thermal_hall_conductivity(magnon_bc;
    temperatures,
    units,
    components=(:xy,),
    axis_kwargs=(title="VCl2 magnon thermal Hall",),
    size=(760, 480),
)
save(joinpath(@__DIR__, "vcl2_magnon_thermal_hall.png"), fig_mthe)

fig_hthe = plot_thermal_hall_conductivity(hybrid_bc;
    temperatures,
    units,
    components=(:xy,),
    axis_kwargs=(title="VCl2 hybrid thermal Hall",),
    size=(760, 480),
)
save(joinpath(@__DIR__, "vcl2_hybrid_thermal_hall.png"), fig_hthe)

fig_pthe = plot_thermal_hall_conductivity(phonon_bc;
    temperatures,
    units,
    components=(:xy,),
    axis_kwargs=(title="VCl2 phonon thermal Hall",),
    size=(760, 480),
)
save(joinpath(@__DIR__, "vcl2_phonon_thermal_hall.png"), fig_pthe)

println("计算 T=$(HALL_TEMPERATURE) K 的粗网格热霍尔。")
kappa_magnon = thermal_hall_conductivity(magnon_bc;
    temperature=HALL_TEMPERATURE,
    component=:xy,
    units=units,
    spread_tol=SPREAD_TOL,
)
kappa_phonon = thermal_hall_conductivity(phonon_bc;
    temperature=HALL_TEMPERATURE,
    component=:xy,
    units=units,
)
kappa_hybrid_coarse = hybrid_thermal_hall_conductivity(hybrid_bc;
    temperature=HALL_TEMPERATURE,
    component=:xy,
    units=units,
    spread_tol=SPREAD_TOL,
)

println("用 Wannier90 式自适应网格重新计算 hybrid 热霍尔。")
kappa_hybrid_adaptive = hybrid_thermal_hall_conductivity(hsys, hall_grid;
    method=:kubo,
    fdstep=FDSTEP,
    degeneracy_handling=:hybrid,
    degeneracy_tol=DEGENERACY_TOL,
    spread_tol=SPREAD_TOL,
    bloch_embedding=:basis_positions,
    temperature=HALL_TEMPERATURE,
    component=:xy,
    units=units,
    adaptive_refinement=ADAPTIVE_REFINEMENT,
    return_diagnostics=true,
)

println()
println("===== VCl2 热霍尔结果汇总 =====")
println("T = $(HALL_TEMPERATURE) K")
println("magnon coarse κ_xy = ", kappa_magnon)
println("phonon coarse κ_xy = ", kappa_phonon)
println("hybrid coarse κ_xy = ", kappa_hybrid_coarse)
println("hybrid adaptive κ_xy = ", kappa_hybrid_adaptive.value)
println("hybrid adaptive refined points = ",
    kappa_hybrid_adaptive.refined_points, " / ", kappa_hybrid_adaptive.total_points)

open(joinpath(@__DIR__, "vcl2_thermal_hall_summary.csv"), "w") do io
    println(io, "quantity,value")
    println(io, "temperature_K,$HALL_TEMPERATURE")
    println(io, "magnon_coarse_kxy,$kappa_magnon")
    println(io, "phonon_coarse_kxy,$kappa_phonon")
    println(io, "hybrid_coarse_kxy,$kappa_hybrid_coarse")
    println(io, "hybrid_adaptive_kxy,$(kappa_hybrid_adaptive.value)")
    println(io, "adaptive_kmesh,$(kappa_hybrid_adaptive.kmesh)")
    println(io, "adaptive_threshold,$(kappa_hybrid_adaptive.threshold)")
    println(io, "hybrid_refined_points,$(kappa_hybrid_adaptive.refined_points)")
    println(io, "total_points,$(kappa_hybrid_adaptive.total_points)")
end

println("已保存 vcl2_thermal_hall_summary.csv。")

using Pkg
Pkg.activate("d:/Trae/MapBerry/MagBerry")
using Sunny,MagBerry,LinearAlgebra,GLMakie
units = Units(:meV, :angstrom)
# ACS nano dual topology
# ================= 参数设置 =================  
S = 1.0       # 自旋大小  
J = -103       # 最近邻交换作用 (铁磁)  
A = -0.68      # 单离子各向异性  
Dz = 1.0      # 最近邻 DMI  
a = 1.0       # 晶格常数  
  
# ================= 晶格设置 =================  
latvecs = lattice_vectors(a, a, 10*a, 90, 90, 120)  
positions = [[0, 0, 0.5],[1/3,2/3,0.5]]  # 单原子基元  
cryst = Crystal(latvecs, positions,1)  
  
# ================= 系统创建 =================  
# 创建自旋系统  
sys = System(cryst, [1 => Moment(s=S, g=2),2 => Moment(s=S, g=2)], :dipole)  
  
# ================= 相互作用设置 =================  

set_exchange!(sys, J, Bond(1, 2, [0, 0, 0]))   
set_exchange!(sys, J, Bond(2, 1, [1, 1, 0]))   
set_exchange!(sys, J, Bond(1, 2, [0, -1, 0]))  

set_exchange!(sys, dmvec([0, 0, -Dz]), Bond(1, 1, [-1, 0, 0]))
set_exchange!(sys, dmvec([0, 0, -Dz]), Bond(1, 1, [0, -1, 0]))
set_exchange!(sys, dmvec([0, 0, -Dz]), Bond(1, 1, [1, 1, 0]))

set_exchange!(sys, dmvec([0, 0, Dz]), Bond(2, 2, [-1, 0, 0]))
set_exchange!(sys, dmvec([0, 0, Dz]), Bond(2, 2, [0, -1, 0]))
set_exchange!(sys, dmvec([0, 0, Dz]), Bond(2, 2, [1, 1, 0]))

# 单离子各向异性  
set_onsite_coupling!(sys, S -> A*S[3]^2, 1)
set_onsite_coupling!(sys, S -> A*S[3]^2, 2)

swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
kpoints = [  
    [0, 0, 0],      # Γ'  
    [1/3, 1/3, 0],  # K  
    [1/2, 0, 0],    # M  
    [0, 0, 0]       # Γ'  
]
kpath = q_space_path(cryst, kpoints, 5000)
res = intensities_bands(swt, kpath)
fig1 = plot_intensities(res; units)
save("bulk_bands.png", fig1)

# ================= 边界态计算 (MagBerry) =================
# 我们沿着一个方向（例如 ky 即 periodic_dim=2）计算色散
# 并在另一个方向（例如 kx 即 finite_dim=1）进行截断，构建条带(Strip)结构
k_path_edge = range(-0.5, 0.5, length=300)

# 计算条带的能带
# finite_dim=1 (沿着a1方向截断)，periodic_dim=2 (沿着a2方向周期性)
# width=40 指条带宽度为 40 个原胞
edge_bands = magnon_strip_bands(swt, k_path_edge; 
                                width=30, 
                                finite_dim=2, 
                                periodic_dim=1, 
                                n_sampling=100)


fig_edge_bands = plot_strip_bands(edge_bands; 
                                  title="Magnon Strip Bands (Edge States)", 
                                  colormap=:RdBu)
save("edge_bands_discrete.png", fig_edge_bands)

fig_edge_spectral = plot_spectral_function(edge_bands; 
                                           broadening=10.0,   
                                           projection=:left, # 可选 :left, :right, :bulk, :total
                                           colormap=:magma,
                                           y_res=500)
save("edge_spectral_density.png", fig_edge_spectral)

# ================= 热霍尔效应计算 (MagBerry) =================
# 使用 2D q-space 网格计算贝利曲率
nx,ny=100,100
grid_2d = q_space_grid(cryst, 
    [1, 0, 0], range(-0.5, 0.5, nx),
    [0, 1, 0], range(-0.5, 0.5, ny)
)

temps = range(0.1, 300, length=30)

bc_data = magnon_berry_curvature(swt, grid_2d)

for b in 1:2
    fig = plot_berry_curvature_plane(bc_data; band=b, component=:z, colormap=:cool, 
                                     bz_color=:black, bz_linewidth=1.0, bz_linestyle=:dash, levels=200)
    save("2D_BC_Band_$(b).png", fig)
end

fig_intrinsic_the = plot_thermal_hall_conductivity(bc_data; 
    temperatures=temps, 
    units=units, 
    components=(:xy,),
    axis_kwargs=(title="Intrinsic Thermal Hall Conductivity",)
)

save("intrinsic_thermal_hall_conductivity.png", fig_intrinsic_the)
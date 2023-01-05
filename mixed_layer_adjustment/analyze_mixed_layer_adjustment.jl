using Oceananigans
using Oceananigans.Units
using GLMakie
using Statistics

filename = "mixed_layer_baroclinic_adjustment"
#filename = "mixed_layer_baroclinic_equilibrium"

bt_t = FieldTimeSeries(filename * "_top_slice.jld2", "b")
ζt_t = FieldTimeSeries(filename * "_top_slice.jld2", "ζ")
bb_t = FieldTimeSeries(filename * "_bottom_slice.jld2", "b")
ζb_t = FieldTimeSeries(filename * "_bottom_slice.jld2", "ζ")
B_t = FieldTimeSeries(filename * "_zonal_average.jld2", "b")
U_t = FieldTimeSeries(filename * "_zonal_average.jld2", "u")
K_t = FieldTimeSeries(filename * "_zonal_average.jld2", "k")

grid = K_t.grid
xc, yc, zc = nodes(bt_t)
xf, yf, zf = nodes(ζt_t)
zf = znodes(Face, grid)[2:grid.Nz]

times = bt_t.times
Nt = length(times)

fig = Figure(resolution = (1800, 1200))
slider = Slider(fig[1, 1:2], range=1:Nt, startvalue=1)
n = slider.value

axbt = Axis(fig[2, 1])
btn = @lift interior(bt_t[$n], :, :, 1)
heatmap!(axbt, xc, yc, btn, colormap=:thermal)

axbb = Axis(fig[3, 1])
bbn = @lift interior(bb_t[$n], :, :, 1)
heatmap!(axbb, xc, yc, bbn, colormap=:thermal)

ζlim = 1e-4
axζt = Axis(fig[2, 2])
ζtn = @lift interior(ζt_t[$n], :, :, 1)
hm = heatmap!(axζt, xf, yf, ζtn, colormap=:balance, colorrange=(-ζlim, ζlim))
Colorbar(fig[2, 3], hm, label="Surface vorticity (s⁻¹)")

axζb = Axis(fig[3, 2])
ζbn = @lift interior(ζb_t[$n], :, :, 1)
hm = heatmap!(axζb, xf, yf, ζbn, colormap=:balance, colorrange=(-ζlim, ζlim))
Colorbar(fig[3, 3], hm, label="Bottom vorticity (s⁻¹)")

Un = @lift interior(U_t[$n], 1, :, :)
Bn = @lift interior(B_t[$n], 1, :, :)
Kn = @lift interior(K_t[$n], 1, :, :)

N² = @lift begin
    N²_op = ∂z(B_t[$n])
    N² = compute!(Field(N²_op))
    interior(N², 1, :, 2:grid.Nz)
end

Uzn = @lift begin
    Uz_op = ∂z(U_t[$n])
    Uz = compute!(Field(Uz_op))
    interior(Uz, 1, :, 2:grid.Nz)
end

Kzn = @lift begin
    Kz_op = ∂z(K_t[$n])
    Kz = compute!(Field(Kz_op))
    interior(Kz, 1, :, 2:grid.Nz)
end

Uzlim = 1e-3
Ulim = 0.2
Klim = 0.2
Kzlim = 1e-3

axu = Axis(fig[4, 2])
hm = heatmap!(axu, yc, zc, Un, colormap=:balance, colorrange=(-Ulim, Ulim))
Colorbar(fig[4, 3], hm, label="Zonal velocity (m s⁻¹)")
contour!(axu, yc, zc, Bn; levels=15, linewidth=2, color=:black)

axs = Axis(fig[5, 2])
hm = heatmap!(axs, yc, zf, Uzn, colormap=:balance, colorrange=(-Uzlim, Uzlim))
Colorbar(fig[5, 3], hm, label="Zonal shear (s⁻¹)")
contour!(axs, yc, zc, Bn; levels=15, linewidth=2, color=:black)

axN = Axis(fig[6, 2])
hm = heatmap!(axN, yc, zf, N², colorrange=(0, 5e-5))
Colorbar(fig[6, 3], hm, label="Buoyancy frequency (s⁻²)")
contour!(axN, yc, zc, Bn; levels=15, linewidth=2, color=:black)

axk = Axis(fig[4, 1])
hm = heatmap!(axk, yc, zc, Kn, colormap=:thermal, colorrange=(0, Klim))
Colorbar(fig[4, 0], hm, label="Eddy kinetic energy (m² s⁻²)", flipaxis=false)
contour!(axu, yc, zc, Bn; levels=15, linewidth=2, color=:black)

axkz = Axis(fig[5, 1])
hm = heatmap!(axkz, yc, zf, Kzn, colormap=:balance, colorrange=(-Kzlim, Kzlim))
Colorbar(fig[5, 0], hm, label="EKE vertical gradient (m s⁻²)", flipaxis=false)
contour!(axs, yc, zc, Bn; levels=15, linewidth=2, color=:black)

display(fig)

record(fig, filename * ".mp4", 1:Nt, framerate=12) do nn
    n[] = nn
end

#=
Nx, Ny, Nz = size(grid)
K_zt = zeros(Nz, Nt)

for n = 1:Nt
    K_zt[:, n] .= interior(mean(K_t[n], dims=2), 1, 1, :)
end

Klim = maximum(abs, K_zt) / 2

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time (days)", ylabel="z (m)")
ax2 = Axis(fig[1, 2], xlabel="Mesoscale kinetic energy (m² s⁻²)", ylabel="z (m)")
heatmap!(ax, times ./ day, zc, permutedims(K_zt, (2, 1)), colorrange=(0, Klim))

ns = [80, 120, 150, 160]

for n in ns 
    lines!(ax2, K_zt[:, n], zc, label="t = " * prettytime(times[n]))
end

axislegend(ax2, position=:rb)

ylims!(ax2, -2000, 0)

display(fig)
=#

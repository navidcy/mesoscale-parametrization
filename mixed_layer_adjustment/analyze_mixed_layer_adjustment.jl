using Oceananigans
using Oceananigans.Units
using GLMakie

filename = "mixed_layer_baroclinic_adjustment"

fig = Figure(resolution = (2400, 1600))

bt_t = FieldTimeSeries(filename * "_top_slice.jld2", "b")
ζt_t = FieldTimeSeries(filename * "_top_slice.jld2", "ζ")
bb_t = FieldTimeSeries(filename * "_bottom_slice.jld2", "b")
ζb_t = FieldTimeSeries(filename * "_bottom_slice.jld2", "ζ")
B_t = FieldTimeSeries(filename * "_zonal_average.jld2", "b")
U_t = FieldTimeSeries(filename * "_zonal_average.jld2", "u")

xc, yc, zc = nodes(b_t)
xf, yf, zf = nodes(ζ_t)
zf = znodes(Face, grid)[2:grid.Nz]

times = b_t.times
Nt = length(times)
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

Uzlim = 5e-4
Ulim = 0.05

axu = Axis(fig[4, 1:2])
hm = heatmap!(axu, yc, zc, Un, colormap=:balance, colorrange=(-Ulim, Ulim))
Colorbar(fig[4, 3], hm, label="Zonal velocity (m s⁻¹)")
contour!(axu, yc, zc, Bn; levels=15, linewidth=2, color=:black)

axs = Axis(fig[5, 1:2])
hm = heatmap!(axs, yc, zf, Uzn, colormap=:balance, colorrange=(-Uzlim, Uzlim))
Colorbar(fig[5, 3], hm, label="Zonal shear (s⁻¹)")
contour!(axs, yc, zc, Bn; levels=15, linewidth=2, color=:black)

axN = Axis(fig[6, 1:2])
hm = heatmap!(axN, yc, zf, N², colorrange=(0, 1e-5))
Colorbar(fig[6, 3], hm, label="Buoyancy frequency (s⁻²)")
contour!(axN, yc, zc, Bn; levels=15, linewidth=2, color=:black)

display(fig)


using Oceananigans
using Oceananigans.Units
using GLMakie

filename = "mixed_layer_baroclinic_adjustment"

fig = Figure(resolution = (800, 800))

b_t = FieldTimeSeries(filename * "_top_slice.jld2", "b")
ζ_t = FieldTimeSeries(filename * "_top_slice.jld2", "ζ")
B_t = FieldTimeSeries(filename * "_zonal_average.jld2", "b")
U_t = FieldTimeSeries(filename * "_zonal_average.jld2", "u")

xc, yc, zc = nodes(b_t)
xf, yf, zf = nodes(ζ_t)
zf = znodes(Face, grid)[2:grid.Nz]

times = b_t.times
Nt = length(times)
slider = Slider(fig[1, 1:2], range=1:Nt, startvalue=1)
n = slider.value

axb = Axis(fig[2, 1])
bn = @lift interior(b_t[$n], :, :, 1)
heatmap!(axb, xc, yc, bn)

ζlim = 1e-4
axζ = Axis(fig[2, 2])
ζn = @lift interior(ζ_t[$n], :, :, 1)
heatmap!(axζ, xf, yf, ζn, colormap=:balance, colorrange=(-ζlim, ζlim))

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

Uzlim = 1e-4
Ulim = 0.03

axu = Axis(fig[3, 1:2])
hm = heatmap!(axu, yc, zc, Un, colormap=:balance, colorrange=(-Ulim, Ulim))
Colorbar(fig[3, 3], hm)
contour!(axu, yc, zc, Bn; levels=15, linewidth=2, color=:black)

axs = Axis(fig[4, 1:2])
hm = heatmap!(axs, yc, zf, Uzn, colormap=:balance, colorrange=(-Uzlim, Uzlim))
Colorbar(fig[4, 3], hm)
contour!(axs, yc, zc, Bn; levels=15, linewidth=2, color=:black)

axN = Axis(fig[5, 1:2])
hm = heatmap!(axN, yc, zf, N², colorrange=(0, 2e-6))
Colorbar(fig[5, 3], hm)
contour!(axN, yc, zc, Bn; levels=15, linewidth=2, color=:black)

display(fig)



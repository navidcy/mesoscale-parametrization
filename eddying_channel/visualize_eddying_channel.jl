using Oceananigans
using GLMakie

fileprefix = "eddying_channel_catke"

xy_filename = fileprefix * "_top_slice.jld2"
xy_bottom_filename = fileprefix * "_bottom_slice.jld2"
yz_filename = fileprefix * "_east_slice.jld2"

et_xy = FieldTimeSeries(xy_filename, "e")
et_yz = FieldTimeSeries(yz_filename, "e")
bt_xy = FieldTimeSeries(xy_filename, "b")
bt_yz = FieldTimeSeries(yz_filename, "b")
Nt_yz = FieldTimeSeries(yz_filename, "N²")
ut_xy = FieldTimeSeries(xy_filename, "u")
ut_yz = FieldTimeSeries(yz_filename, "u")
vt_xy = FieldTimeSeries(xy_filename, "v")
vt_yz = FieldTimeSeries(yz_filename, "v")
ζt_xy = FieldTimeSeries(xy_filename, "ζ")
kt_xy = FieldTimeSeries(xy_filename, "k")
kt_yz = FieldTimeSeries(yz_filename, "k")
ζb_xy = FieldTimeSeries(xy_bottom_filename, "ζ")
kb_xy = FieldTimeSeries(xy_bottom_filename, "k")

times = et_xy.times
grid = et_xy.grid
x, y, z = nodes(et_xy)
xζ, yζ, zζ = nodes(ζt_xy)
zf = znodes(Face, grid)

x = x ./ 1e3
y = y ./ 1e3

xζ = xζ ./ 1e3
yζ = yζ ./ 1e3

Nt = length(times)

zlim = -450
klim = 0.5
ζlim = 2e-5
elim = 2e-4
Nlim = 5e-5
ulim = 0.5

fig = Figure(resolution=(1600, 1200))
slider = Slider(fig[1, 1:6], range=1:Nt, startvalue=1)
n = slider.value

titlestr = @lift string("Eddying channel spinup at t = ", prettytime(times[$n]))
Label(fig[2, 1:6], titlestr)

axe = Axis(fig[4, 1], xlabel="x (km)", ylabel="y (km)", aspect=1)
en_xy = @lift interior(et_xy[$n], :, :, 1)
hm = heatmap!(axe, x, y, en_xy, colorrange=(0, 1.2elim))
Colorbar(fig[3, 1], hm, flipaxis=true, vertical=false, label="Turbulent kinetic energy (m² s⁻²)", ticks=[1e-4, 2e-4])

axζ = Axis(fig[4, 2], title="Top", xlabel="x (km)", ylabel="y (km)", aspect=1)
ζn_xy = @lift interior(ζt_xy[$n], :, :, 1)
hm = heatmap!(axζ, x, y, ζn_xy, colormap=:balance, colorrange=(-ζlim, ζlim))
Colorbar(fig[3, 2:3], hm, flipaxis=true, vertical=false, label="Vertical vorticity (s⁻¹)")

axζ = Axis(fig[4, 3], title="Bottom", xlabel="x (km)", ylabel="y (km)", aspect=1)
ζn_xy = @lift interior(ζb_xy[$n], :, :, 1)
heatmap!(axζ, x, y, ζn_xy, colormap=:balance, colorrange=(-ζlim, ζlim))

axk = Axis(fig[6, 1], title="Top", xlabel="x (km)", ylabel="y (km)", aspect=1)
kn_xy = @lift interior(kt_xy[$n], :, :, 1)
hm = heatmap!(axk, x, y, kn_xy, colormap=:solar, colorrange=(0, klim))
Colorbar(fig[5, 1:2], hm, flipaxis=true, vertical=false, label="Eddy kinetic energy (m s⁻¹)")

axk = Axis(fig[6, 2], title="Bottom", xlabel="x (km)", ylabel="y (km)", aspect=1)
kn_xy = @lift interior(kb_xy[$n], :, :, 1)
heatmap!(axk, x, y, kn_xy, colormap=:solar, colorrange=(0, klim))

axb = Axis(fig[6, 3], title="Buoyancy (m s⁻²)", xlabel="x (km)", ylabel="y (km)", aspect=1)
bn_xy = @lift interior(bt_xy[$n], :, :, 1)
hm = heatmap!(axb, x, y, bn_xy, colormap=:thermal, colorrange=(-1e-4, 1e-3))
Colorbar(fig[5, 3], hm, flipaxis=true, vertical=false)

axe = Axis(fig[7, 1:3], title="Turbulent kinetic energy (m² s⁻²)", xlabel="y (km)", ylabel="z (m)")
ylims!(axe, zlim, 0)
en_yz = @lift interior(et_yz[$n], 1, :, :)
hm = heatmap!(axe, y, z, en_yz, colorrange=(0, elim))
Colorbar(fig[7, 4], hm)

axN = Axis(fig[8, 1:3], title="Buoyancy frequency (s⁻²)", xlabel="y (km)", ylabel="z (m)")
#ylims!(axN, zlim, 0)
Nn_yz = @lift interior(Nt_yz[$n], 1, :, :) #interior(compute!(Field(∂z(bt_yz[$n]))), 1, :, :)
hm = heatmap!(axN, y, zf, Nn_yz, colorrange=(0, Nlim))
Colorbar(fig[8, 4], hm)

axk = Axis(fig[9, 1:3], title="K (m² s⁻²)", xlabel="y (km)", ylabel="z (m)")
kn_yz = @lift interior(kt_yz[$n], 1, :, :)
hm = heatmap!(axv, y, z, kn_yz, colormap=:balance, colorrange=(0, elim))
Colorbar(fig[9, 4], hm, flipaxis=true)

#=
axu = Axis(fig[8, 1:6], title="u", xlabel="y (km)", ylabel="z (m)")
#ylims!(axu, zlim, 0)
un_yz = @lift interior(ut_yz[$n], 1, :, :)
hm = heatmap!(axu, y, z, un_yz, colormap=:balance, colorrange=(-ulim, ulim))
Colorbar(fig[8, 7], hm)
=#

#=
axκ = Axis(fig[5, 1], title="TKE diffusivity (m s⁻²)", xlabel="y (km)", ylabel="z (m)")
ylims!(axκ, zlim, 0)
κn_yz = @lift interior(κet_yz[$n], 1, :, :)
hm = heatmap!(axκ, y, zf, κn_yz, colorrange=(0, 3))
Colorbar(fig[5, 0], hm, flipaxis=false)

axℓ = Axis(fig[5, 2], title="Mixing length for TKE (m)", xlabel="y (km)", ylabel="z (m)")
ylims!(axℓ, zlim, 0)
ℓen_yz = @lift interior(ℓet_yz[$n], 1, :, :)
hm = heatmap!(axℓ, y, zf, ℓen_yz, colorrange=(0, 200))
Colorbar(fig[5, 3], hm)
=#

display(fig)

#record(fig, "eddying_channel_spinup.mp4", 1:Nt, framerate=6) do nn
#    n[] = nn
#end


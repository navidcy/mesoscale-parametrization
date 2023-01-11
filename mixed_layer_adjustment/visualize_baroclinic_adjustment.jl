using Oceananigans
using GLMakie

fileprefix = "simple_baroclinic_adjustment"

xy_filename = fileprefix * "_top_slice.jld2"
xy_bottom_filename = fileprefix * "_bottom_slice.jld2"
yz_filename = fileprefix * "_east_slice.jld2"
zonal_average_filename = fileprefix * "_zonal_average.jld2"

Bt = FieldTimeSeries(zonal_average_filename, "b")
Ut = FieldTimeSeries(zonal_average_filename, "u")
St = FieldTimeSeries(zonal_average_filename, "uz")
Kt = FieldTimeSeries(zonal_average_filename, "k")

ζt_xy = FieldTimeSeries(xy_filename, "ζ")

times = ζt_xy.times
grid = ζt_xy.grid
xζ, yζ, zζ = nodes(ζt_xy)
x, y, z = nodes(Bt)
xz, yz, zz = nodes(St)

x = x ./ 1e3
y = y ./ 1e3

xζ = xζ ./ 1e3
yζ = yζ ./ 1e3

Nt = length(times)

fig = Figure(resolution=(2400, 1200))
aspect = 8
slider = Slider(fig[1, 1], range=1:Nt, startvalue=1)
n = slider.value

titlestr = @lift string("Eddying channel spinup at t = ", prettytime(times[$n]))
Label(fig[2, 1], titlestr, tellwidth=false)

axz = Axis(fig[3, 1]; xlabel="x (km)", ylabel="y (km)", aspect)
ζn_xy = @lift interior(ζt_xy[$n], :, :, 1)
hm = heatmap!(axz, xζ, yζ, ζn_xy, colormap=:balance, colorrange=(-1e-4, 1e-4))
Colorbar(fig[3, 2], hm, vertical=true, label="Vertical vorticity (s⁻¹)")

bn = @lift interior(Bt[$n], 1, :, :)

axu = Axis(fig[4, 1]; xlabel="x (km)", ylabel="y (km)", aspect)
Un = @lift interior(Ut[$n], 1, :, :)
hm = heatmap!(axu, y, z, Un, colormap=:balance, colorrange=(-0.5, 0.5))
contour!(axu, y, z, bn, levels=25, color=:black, linewidth=2)
Colorbar(fig[4, 2], hm, label="Zonally-averaged zonal velocity (m s⁻¹)")

axs = Axis(fig[5, 1]; xlabel="x (km)", ylabel="y (km)", aspect)
Sn = @lift interior(St[$n], 1, :, :)
hm = heatmap!(axs, yz, zz, Sn, colormap=:balance, colorrange=(-1e-3, 1e-3))
contour!(axs, y, z, bn, levels=25, color=:black, linewidth=2)
Colorbar(fig[5, 2], hm, label="Zonally-averaged shear (s⁻¹)")

axk = Axis(fig[6, 1]; xlabel="x (km)", ylabel="y (km)", aspect)
Kn = @lift interior(Kt[$n], 1, :, :)
hm = heatmap!(axk, y, z, Kn, colormap=:solar, colorrange=(0.1, 0.5))
contour!(axk, y, z, bn, levels=25, color=:black, linewidth=2)
Colorbar(fig[6, 2], hm, label="Zonally-averaged eddy kinetic energy (m² s⁻²)")

display(fig)

record(fig, "eddying_channel_spinup.mp4", 1:Nt, framerate=6) do nn
    n[] = nn
end


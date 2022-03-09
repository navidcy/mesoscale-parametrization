using Oceananigans
using GLMakie
using JLD2

filename = "eddying_channel_catke"

zonal_file = jldopen(filename * "_zonal_average.jld2")

B = FieldTimeSeries(filename * "_zonal_average.jld2", "b")
U = FieldTimeSeries(filename * "_zonal_average.jld2", "u")
V = FieldTimeSeries(filename * "_zonal_average.jld2", "v")
bb = FieldTimeSeries(filename * "_zonal_average.jld2", "bb")
wb = FieldTimeSeries(filename * "_zonal_average.jld2", "wb")
vb = FieldTimeSeries(filename * "_zonal_average.jld2", "vb")

fig = Figure(resolution=(1200, 1600))

axbb = Axis(fig[1, 1])
axwb = Axis(fig[2, 1])
axvb = Axis(fig[3, 1])
axN = Axis(fig[4, 1])
axS = Axis(fig[5, 1])

Nt = length(bb.times)

slider = Slider(fig[6, :], range=1:Nt, startvalue=1)
n = slider.value

x, y, z = nodes(bb)

bbn = @lift interior(bb[$n])[1, :, :]
wbn = @lift interior(wb[$n])[1, :, :]
vbn = @lift interior(vb[$n])[1, :, :]
Bn = @lift interior(B[$n])[1, :, :]
Nsq = @lift begin
    Bn = B[$n]
    Nsq = Field(∂z(Bn))
    compute!(Nsq)
    interior(Nsq)[1, :, :]  
end

Sh = @lift begin
    Un = U[$n]
    # Sh = Field(sqrt(∂z(Un)^2))
    Vn = V[$n]
    Sh = Field(sqrt(∂z(Un)^2 + ∂z(Vn)^2))
    compute!(Sh)
    interior(Sh)[1, :, :]
end

buoyancycontours!(ax) = contour!(ax, y * 1e-3, z, Bn; levels=25, color = (:black, 0.6), linewidth = 3)

heatmap!(axbb, y * 1e-3, z, bbn, colormap=:speed)
buoyancycontours!(axbb)

heatmap!(axwb, y * 1e-3, z, wbn, colormap=:curl)
buoyancycontours!(axwb)

heatmap!(axvb, y * 1e-3, z, vbn, colormap=:curl)
buoyancycontours!(axvb)

heatmap!(axN, y * 1e-3, z, Nsq, colormap=:thermal, colorrange=(1e-6, 1e-5))
buoyancycontours!(axN)

heatmap!(axS, y * 1e-3, z, Sh, colormap=:solar, colorrange=(0, 1e-4))
buoyancycontours!(axS)

display(fig)

#=
record(fig, "eddy_fluxes.mp4", 1:Nt; framerate=12) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end
=#

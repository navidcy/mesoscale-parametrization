# using Pkg
# pkg"add Oceananigans CairoMakie JLD2"
# pushfirst!(LOAD_PATH, @__DIR__)
# pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..", "..")) # add Oceananigans

using Printf
using Statistics
using JLD2

using Oceananigans
using Oceananigans.Units
using Oceananigans.OutputReaders: FieldTimeSeries
using Oceananigans.Grids: xnode, ynode, znode
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity

using Random
Random.seed!(1234)

arch = GPU()

filename = "eddying_channel_catke"

# Domain
const Lx = 2000kilometers # zonal domain length [m]
const Ly = 2000kilometers # meridional domain length [m]
const Lz = 2kilometers    # depth [m]

# number of grid points
Nx = 128
Ny = 128
Nz = 60

save_fields_interval = 7days
stop_time = 20years + 1day
Δt₀ = 5minutes

# stretched grid

# we implement here a linearly streched grid in which the top grid cell has Δzₜₒₚ
# and every other cell is bigger by a factor σ, e.g.,
# Δzₜₒₚ, Δzₜₒₚ * σ, Δzₜₒₚ * σ², ..., Δzₜₒₚ * σᴺᶻ⁻¹,
# so that the sum of all cell heights is Lz

# Given Lz and stretching factor σ > 1 the top cell height is Δzₜₒₚ = Lz * (σ - 1) / σ^(Nz - 1)

σ = 1.04 # linear stretching factor
Δz_center_linear(k) = Lz * (σ - 1) * σ^(Nz - k) / (σ^Nz - 1) # k=1 is the bottom-most cell, k=Nz is the top cell
linearly_spaced_faces(k) = k==1 ? -Lz : - Lz + sum(Δz_center_linear.(1:k-1))

refinement = 2 # controls spacing near surface (higher means finer spaced)
stretching = 4  # controls rate of stretching at bottom

# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz

# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement

# Bottom-intensified stretching function
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(arch;
                       topology = (Periodic, Bounded, Bounded),
                       size = (Nx, Ny, Nz),
                       halo = (3, 3, 3),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = (-Lz, 0)) # z_faces)

# The vertical spacing versus depth for the prescribed grid
# using GLMakie
# plot(grid.Δzᵃᵃᶜ[1:Nz], grid.zᵃᵃᶜ[1:Nz], marker = :circle,
#      axis=(xlabel = "Vertical spacing (m)",
#            ylabel = "Depth (m)"))

@info "Built a grid: $grid."

#####
##### Boundary conditions
#####

α  = 2e-4     # [K⁻¹] thermal expansion coefficient
g  = 9.8061   # [m s⁻²] gravitational constant
cᵖ = 3994.0   # [J K⁻¹] heat capacity
ρ  = 1024.0   # [kg m⁻³] reference density

parameters = (Ly = Ly,
              Lz = Lz,
              Qᵇ = 10 / (ρ * cᵖ) * α * g,          # buoyancy flux magnitude [m² s⁻³]    
              y_shutoff = 5/6 * Ly,                # shutoff location for buoyancy flux [m]
              τ = 0.15/ρ,                          # surface kinematic wind stress [m² s⁻²]
              μ = 1 / 30days,                      # bottom drag damping time-scale [s⁻¹]
              ΔB = 8 * α * g,                      # surface vertical buoyancy gradient [s⁻²]
              H = Lz,                              # domain depth [m]
              h = 1000.0,                          # exponential decay scale of stable stratification [m]
              y_sponge = 19/20 * Ly,               # southern boundary of sponge layer [m]
              λt = 7days                           # relaxation time scale [s]
)

@inline function buoyancy_flux(i, j, grid, clock, model_fields, p)
    y = ynode(Center(), j, grid)
    return ifelse(y < p.y_shutoff, p.Qᵇ * cos(3π * y / p.Ly), 0.0)
end

buoyancy_flux_bc = FluxBoundaryCondition(buoyancy_flux, discrete_form=true, parameters=parameters)

@inline function u_stress(i, j, grid, clock, model_fields, p)
    y = ynode(Center(), j, grid)
    return - p.τ * sin(π * y / p.Ly)
end

u_stress_bc = FluxBoundaryCondition(u_stress, discrete_form=true, parameters=parameters)

@inline u_drag(i, j, grid, clock, model_fields, p) = @inbounds - p.μ * p.Lz * model_fields.u[i, j, 1] 
@inline v_drag(i, j, grid, clock, model_fields, p) = @inbounds - p.μ * p.Lz * model_fields.v[i, j, 1]

u_drag_bc = FluxBoundaryCondition(u_drag, discrete_form=true, parameters=parameters)
v_drag_bc = FluxBoundaryCondition(v_drag, discrete_form=true, parameters=parameters)

b_bcs = FieldBoundaryConditions(top = buoyancy_flux_bc)

u_bcs = FieldBoundaryConditions(top = u_stress_bc, bottom = u_drag_bc)
v_bcs = FieldBoundaryConditions(bottom = v_drag_bc)

#####
##### Coriolis
#####

const f = -1e-4     # [s⁻¹]
const β =  1e-11    # [m⁻¹ s⁻¹]
coriolis = BetaPlane(f₀ = f, β = β)

#####
##### Forcing and initial condition
#####

@inline initial_buoyancy(z, p) = p.ΔB * (exp(z / p.h) - exp(-p.Lz / p.h)) / (1 - exp(-p.Lz / p.h))
@inline mask(y, p) = max(0.0, y - p.y_sponge) / (Ly - p.y_sponge)

@inline function buoyancy_relaxation(i, j, k, grid, clock, model_fields, p)
    timescale = p.λt
    y = ynode(Center(), j, grid)
    z = znode(Center(), k, grid)
    target_b = initial_buoyancy(z, p)
    b = @inbounds model_fields.b[i, j, k]
    return - 1 / timescale  * mask(y, p) * (b - target_b)
end

Fb = Forcing(buoyancy_relaxation, discrete_form = true, parameters = parameters)

# Turbulence closures

κh = 0.5e-5 # [m²/s] horizontal diffusivity
νh = 30.0   # [m²/s] horizontal viscocity
κz = 0.5e-5 # [m²/s] vertical diffusivity
νz = 3e-4   # [m²/s] vertical viscocity

vertical_diffusive_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν = νz, κ = κz)

horizontal_diffusive_closure = HorizontalScalarDiffusivity(ν = νh, κ = κh)

catke = CATKEVerticalDiffusivity()

#####
##### Model building
#####

@info "Building a model..."

model = HydrostaticFreeSurfaceModel(grid = grid,
                                    free_surface = ImplicitFreeSurface(),
                                    momentum_advection = WENO5(),
                                    tracer_advection = WENO5(),
                                    buoyancy = BuoyancyTracer(),
                                    coriolis = coriolis,
                                    closure = (catke, vertical_diffusive_closure, horizontal_diffusive_closure),
                                    tracers = (:b, :c, :e),
                                    boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs),
                                    forcing = (; b=Fb))

@info "Built $model."

#####
##### Initial conditions
#####

# resting initial condition
ε(σ) = σ * randn()
bᵢ(x, y, z) = parameters.ΔB * ( exp(z / parameters.h) - exp(-Lz / parameters.h) ) / (1 - exp(-Lz / parameters.h)) + ε(1e-8)
uᵢ(x, y, z) = ε(1e-8)
vᵢ(x, y, z) = ε(1e-8)
wᵢ(x, y, z) = ε(1e-8)

Δy = 100kilometers
Δz = 100
Δc = 2Δy
cᵢ(x, y, z) = exp(-(y - Ly/2)^2 / 2Δc^2) * exp(-(z + Lz/4)^2 / 2Δz^2)

set!(model, b=bᵢ, u=uᵢ, v=vᵢ, w=wᵢ, c=cᵢ)

#####
##### Simulation building
#####

simulation = Simulation(model, Δt=Δt₀, stop_time=stop_time)

# add timestep wizard callback
wizard = TimeStepWizard(cfl=0.1, max_change=1.1, max_Δt=20minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(20))

# add progress callback
wall_clock = [time_ns()]

function print_progress(sim)
    @printf("[%05.2f%%] i: %d, t: %s, wall time: %s, max(u): (%6.3e, %6.3e, %6.3e) m/s, next Δt: %s\n",
            100 * (sim.model.clock.time / sim.stop_time),
            sim.model.clock.iteration,
            prettytime(sim.model.clock.time),
            prettytime(1e-9 * (time_ns() - wall_clock[1])),
            maximum(abs, sim.model.velocities.u),
            maximum(abs, sim.model.velocities.v),
            maximum(abs, sim.model.velocities.w),
            prettytime(sim.Δt))

    wall_clock[1] = time_ns()
    
    return nothing
end

simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(20))


#####
##### Diagnostics
#####

u, v, w = model.velocities
b, c = model.tracers.b, model.tracers.c
η = model.free_surface.η

ζ = Field(∂x(v) - ∂y(u))

B = Field(Average(b, dims=1))
C = Field(Average(c, dims=1))
U = Field(Average(u, dims=1))
η̄ = Field(Average(η, dims=1))
V = Field(Average(v, dims=1))
W = Field(Average(w, dims=1))

b′ = b - B
u′ = u - U
v′ = v - V
w′ = w - W
c′ = c - C

tke_op = @at (Center, Center, Center) (u′ * u′ + v′ * v′ + w′ * w′) / 2
tke = Field(Average(tke_op, dims=1))

uv_op = @at (Center, Center, Center) u′ * v′
vw_op = @at (Center, Center, Center) v′ * w′
uw_op = @at (Center, Center, Center) u′ * w′

u′v′ = Field(Average(uv_op, dims=1))
v′w′ = Field(Average(vw_op, dims=1))
u′w′ = Field(Average(uw_op, dims=1))

b′b′ = Field(Average(b′ * b′, dims=1))
v′b′ = Field(Average(b′ * v′ , dims=1))
w′b′ = Field(Average(b′ * w′ , dims=1))

c′c′ = Field(Average(c′ * c′, dims=1))
v′c′ = Field(Average(c′ * v′, dims=1))
w′c′ = Field(Average(c′ * w′, dims=1))

outputs = (; b, c, ζ, u, v, w)

zonally_averaged_outputs = (b=B, u=U, v=V, w=W, c=C, η=η̄,
                            vb=v′b′, wb=w′b′, vc=v′c′, wc=w′c′, bb=b′b′,
                            tke=tke, uv=u′v′, vw=v′w′, uw=u′w′, cc=c′c′)

#####
##### Build checkpointer and output writer
#####

simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                        schedule = TimeInterval(10years),
                                                        prefix = filename,
                                                        force = true)

slicers = (west = (1, :, :),
           east = (grid.Nx, :, :),
           south = (:, 1, :),
           north = (:, grid.Ny, :),
           bottom = (:, :, 1),
           top = (:, :, grid.Nz))

for side in keys(slicers)
    indices = slicers[side]

    simulation.output_writers[side] = JLD2OutputWriter(model, outputs;
                                                       schedule = TimeInterval(save_fields_interval),
                                                       indices,
                                                       prefix = filename * "_$(side)_slice",
                                                       force = true)
end

simulation.output_writers[:zonal] = JLD2OutputWriter(model, zonally_averaged_outputs;
                                                     schedule = TimeInterval(save_fields_interval),
                                                     prefix = filename * "_zonal_average",
                                                     force = true)

#=
simulation.output_writers[:zonal] = JLD2OutputWriter(model, zonally_averaged_outputs,
                                                     schedule = AveragedTimeInterval(60days),
                                                     prefix = filename * "_zonal_time_average",
                                                     force = true)
=#

#=
simulation.output_writers[:averages] = JLD2OutputWriter(model, averaged_outputs,
                                                        schedule = AveragedTimeInterval(1days, window=1days, stride=1),
                                                        prefix = "eddying_channel_averages",
                                                        verbose = true,
                                                        force = true)
=#

@info "Running the simulation..."

run!(simulation, pickup=false)

# simulation.stop_time += 61days
# run!(simulation, pickup=true)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)


#####
##### Visualization
#####

# ENV["GKSwstype"] = "100"
# using CairoMakie

#=

using GLMakie
using JLD2

fig = Figure(resolution = (2300, 1400))

ax_b = fig[1, 1] = LScene(fig)
ax_ζ = fig[1, 2] = LScene(fig)

axis_rotation_angles = (π/24, -π/6, 0)

# Extract surfaces on all 6 boundaries

iter = Observable(0)

filename = "eddying_channel_catke"

zonal_file = jldopen(filename * "_zonal_average.jld2")

slicers = (west = FieldSlicer(i=1),
           east = FieldSlicer(i=grid.Nx),
           south = FieldSlicer(j=1),
           north = FieldSlicer(j=grid.Ny),
           bottom = FieldSlicer(k=1),
           top = FieldSlicer(k=grid.Nz))

sides = keys(slicers)

slice_files = NamedTuple(side => jldopen(filename * "_$(side)_slice.jld2") for side in sides)

# Build coordinates, rescaling the vertical coordinate

xζ, yζ, zζ = nodes((Face, Face, Center), grid)
xu, yu, zu = nodes((Face, Center, Center), grid)
xv, yv, zv = nodes((Center, Face, Center), grid)
xw, yw, zw = nodes((Center, Center, Face), grid)
xb, yb, zb = nodes((Center, Center, Center), grid)

zscale = 300
zu = zu .* zscale
zb = zb .* zscale
zζ = zζ .* zscale

zonal_slice_displacement = 1.4

b_slices = (
      west = @lift(Array(slice_files.west["timeseries/b/"   * string($iter)][1, :, :])),
      east = @lift(Array(slice_files.east["timeseries/b/"   * string($iter)][1, :, :])),
     south = @lift(Array(slice_files.south["timeseries/b/"  * string($iter)][:, 1, :])),
     north = @lift(Array(slice_files.north["timeseries/b/"  * string($iter)][:, 1, :])),
    bottom = @lift(Array(slice_files.bottom["timeseries/b/" * string($iter)][:, :, 1])),
       top = @lift(Array(slice_files.top["timeseries/b/"    * string($iter)][:, :, 1]))
)

clims_b = @lift extrema(slice_files.top["timeseries/b/" * string($iter)][:])
kwargs_b = (colorrange=clims_b, colormap=:deep, show_axis=false)

surface!(ax_b, yb, zb, b_slices.west;   transformation = (:yz, xb[1]),   kwargs_b...)
surface!(ax_b, yb, zb, b_slices.east;   transformation = (:yz, xb[end]), kwargs_b...)
surface!(ax_b, xb, zb, b_slices.south;  transformation = (:xz, yb[1]),   kwargs_b...)
surface!(ax_b, xb, zb, b_slices.north;  transformation = (:xz, yb[end]), kwargs_b...)
surface!(ax_b, xb, yb, b_slices.bottom; transformation = (:xy, zb[1]),   kwargs_b...)
surface!(ax_b, xb, yb, b_slices.top;    transformation = (:xy, zb[end]), kwargs_b...)

b_avg = @lift zonal_file["timeseries/b/" * string($iter)][1, :, :]
u_avg = @lift zonal_file["timeseries/u/" * string($iter)][1, :, :]

clims_u = @lift extrema(zonal_file["timeseries/u/" * string($iter)][1, :, :])
clims_u = (-0.4, 0.4)

surface!(ax_b, yu, zu, u_avg; transformation = (:yz, zonal_slice_displacement * xu[end]), colorrange=clims_u, colormap=:balance)
contour!(ax_b, yb, zb, b_avg; levels = 25, color = :black, linewidth = 2, transformation = (:yz, zonal_slice_displacement * xb[end]), show_axis=false)

rotate_cam!(ax_b.scene, axis_rotation_angles)

ζ_slices = (
      west = @lift(Array(slice_files.west["timeseries/ζ/"   * string($iter)][1, :, :])),
      east = @lift(Array(slice_files.east["timeseries/ζ/"   * string($iter)][1, :, :])),
     south = @lift(Array(slice_files.south["timeseries/ζ/"  * string($iter)][:, 1, :])),
     north = @lift(Array(slice_files.north["timeseries/ζ/"  * string($iter)][:, 1, :])),
    bottom = @lift(Array(slice_files.bottom["timeseries/ζ/" * string($iter)][:, :, 1])),
       top = @lift(Array(slice_files.top["timeseries/ζ/"    * string($iter)][:, :, 1]))
)

# u_slices = (
#       west = @lift(Array(slice_files.west["timeseries/u/"   * string($iter)][1, :, :])),
#       east = @lift(Array(slice_files.east["timeseries/u/"   * string($iter)][1, :, :])),
#      south = @lift(Array(slice_files.south["timeseries/u/"  * string($iter)][:, 1, :])),
#      north = @lift(Array(slice_files.north["timeseries/u/"  * string($iter)][:, 1, :])),
#     bottom = @lift(Array(slice_files.bottom["timeseries/u/" * string($iter)][:, :, 1])),
#        top = @lift(Array(slice_files.top["timeseries/u/"    * string($iter)][:, :, 1]))
# )

# clims_ζ = @lift extrema(slice_files.west["timeseries/ζ/" * string($iter)][:])
clims_ζ = (-5e-5, 5e-5)
kwargs_ζ = (colorrange = clims_ζ, colormap=:curl, show_axis=false)

surface!(ax_ζ, yζ, zζ, ζ_slices.west;   transformation = (:yz, xζ[1]),   kwargs_ζ...)
surface!(ax_ζ, yζ, zζ, ζ_slices.east;   transformation = (:yz, xζ[end]), kwargs_ζ...)
surface!(ax_ζ, xζ, zζ, ζ_slices.south;  transformation = (:xz, yζ[1]),   kwargs_ζ...)
surface!(ax_ζ, xζ, zζ, ζ_slices.north;  transformation = (:xz, yζ[end]), kwargs_ζ...)
surface!(ax_ζ, xζ, yζ, ζ_slices.bottom; transformation = (:xy, zζ[1]),   kwargs_ζ...)
surface!(ax_ζ, xζ, yζ, ζ_slices.top;    transformation = (:xy, zζ[end]), kwargs_ζ...)

b_avg = @lift zonal_file["timeseries/b/" * string($iter)][1, :, :]
u_avg = @lift zonal_file["timeseries/u/" * string($iter)][1, :, :]

surface!(ax_ζ, yu, zu, u_avg; transformation = (:yz, zonal_slice_displacement * xu[end]), colorrange=clims_u, colormap=:balance)
contour!(ax_ζ, yb, zb, b_avg; levels = 25, color = :black, linewidth = 2, transformation = (:yz, zonal_slice_displacement * xb[end]), show_axis=false)

rotate_cam!(ax_ζ.scene, axis_rotation_angles)

title = @lift(string("Buoyancy and relative vorticity at t = ",
                     prettytime(zonal_file["timeseries/t/" * string($iter)])))

fig[0, :] = Label(fig, title, textsize=50)

iterations = parse.(Int, keys(zonal_file["timeseries/t"]))

record(fig, "eddying_channel_convadj.mp4", iterations, framerate=12) do i
    @info "Plotting iteration $i of $(iterations[end])..."
    iter[] = i
end

for file in slice_files
    close(file)
end

close(zonal_file)
=#

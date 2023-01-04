using Oceananigans
using Oceananigans.Units
using Printf
using CairoMakie

# f uz = - by
@inline ramp(y, Δy) = (1 + tanh(y / Δy)) / 2
d_ramp_dy(y, Δy) = sech(y / Δy)^2 / (2Δy)

architecture = GPU()
α = 4
Nz = 128
Ny = 2Nz
Nx = 2α * Ny
Lx = α * 1000kilometers # east-west extent [m]
Ly = 1000kilometers     # north-south extent [m]
Lz = 1000meters         # depth [m]
N² = 1e-5              # [s⁻²] buoyancy frequency / stratification
M² = 2e-7              # [s⁻²] horizontal buoyancy gradient
Δy = 100kilometers      # width of the region of the front
Δb = Δy * M²           # buoyancy jump associated with the front
ϵᵇ = 1e-2 * Δb         # noise amplitude
f = -1e-4
Cᴰ = 2e-3
Δt₀ = 2minutes
stop_time = 80days
wizard_iteration_interval = 20
progress_iteration_interval = 100
filename = "mixed_layer_baroclinic_adjustment"
mixed_layer_depth = h = 200meters
shear_layer_depth = H = 500meters
save_fields_interval = 0.5day

parameters = (; Cᴰ, N², Δb, Δy, h, H)

@inline unstratified_b(x, y, z, p) = p.Δb * ramp(y, p.Δy)
@inline linearly_stratified_b(x, y, z, p) = p.N² * z + p.Δb * ramp(y, p.Δy)

@inline function piecewise_stratified_b(x, y, z, p)
    by = p.Δb * ramp(y, p.Δy) * exp(z / p.H)
    bz = ifelse(z < -p.h, p.N² * (z + p.h), zero(z))
    return by + bz
end

#bᵢ(x, y, z) = unstratified_b(x, y, z, parameters) + ϵᵇ * randn()
bᵢ(x, y, z) = piecewise_stratified_b(x, y, z, parameters) + ϵᵇ * randn()

# f uz = - by
uᵢ(x, y, z) = - Δb / f * d_ramp_dy(y, Δy) * H * exp(z / H)

grid = RectilinearGrid(architecture;
                       size = (Nx, Ny, Nz),
                       x = (0, Lx),
                       y = (-Ly/2, Ly/2),
                       z = (-Lz, 0),
                       topology = (Periodic, Bounded, Bounded))

@inline u_drag(x, y, t, u, v, p) = - p.Cᴰ * sqrt(u^2 + v^2) * u
@inline v_drag(x, y, t, u, v, p) = - p.Cᴰ * sqrt(u^2 + v^2) * v

u_drag_bc = FluxBoundaryCondition(u_drag; field_dependencies=(:u, :v), parameters)
v_drag_bc = FluxBoundaryCondition(v_drag; field_dependencies=(:u, :v), parameters)

u_bcs = FieldBoundaryConditions(bottom=u_drag_bc)
v_bcs = FieldBoundaryConditions(bottom=v_drag_bc)

boundary_conditions=(u=u_bcs, v=v_bcs)

model = HydrostaticFreeSurfaceModel(; grid, boundary_conditions,
                                    coriolis = FPlane(; f),
                                    buoyancy = BuoyancyTracer(),
                                    tracers = (:b, :r),
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO())

rᵢ(x, y, z) = z
set!(model, b=bᵢ, u=uᵢ, r=rᵢ)

simulation = Simulation(model; Δt=Δt₀, stop_time)
wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=20minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(wizard_iteration_interval))

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

simulation.callbacks[:print_progress] =
    Callback(print_progress,
             IterationInterval(progress_iteration_interval))

b, r = model.tracers.b
u, v, w = model.velocities
ζ = ∂x(v) - ∂y(u)

B = Field(Average(b, dims=1))
R = Field(Average(r, dims=1))
U = Field(Average(u, dims=1))

u′ = u - U
k = (u′^2 + v^2) / 2
K = Field(Average(k, dims=1))

slicers = (west = (1, :, :),
           east = (grid.Nx, :, :),
           south = (:, 1, :),
           north = (:, grid.Ny, :),
           bottom = (:, :, 1),
           top = (:, :, grid.Nz))

for side in keys(slicers)
    indices = slicers[side]

    simulation.output_writers[side] = JLD2OutputWriter(model, (; b, ζ, k);
                                                       filename = filename * "_$(side)_slice",
                                                       schedule = TimeInterval(save_fields_interval),
                                                       overwrite_existing = true,
                                                       indices)
end

simulation.output_writers[:zonal] = JLD2OutputWriter(model, (b=B, u=U, k=K);
                                                     schedule = TimeInterval(save_fields_interval),
                                                     overwrite_existing = true,
                                                     filename = filename * "_zonal_average")

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)

#fig = Figure(resolution = (1200, 800))
#ax = Axis(fig[1, 1])
#heatmap!(ax, interior(b, :, :, grid.Nz))


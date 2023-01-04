using Oceananigans
using Oceananigans.Units
using Printf
using CairoMakie

@inline ramp(y, Δy) = (1 + tanh(y / Δy)) / 2 #min(max(0, y/Δy + 1/2), 1)
d_ramp_dy(y, Δy) = sech(y / Δy)^2 / (2Δy)

architecture = GPU()
Nx = 512
Ny = 128
Nz = 64
Ly = 3000kilometers     # north-south extent [m]
Lx = 4 * Ly            # east-west extent [m]
Lz = 1000meters         # depth [m]
N² = 4e-6              # [s⁻²] buoyancy frequency / stratification
M² = 8e-8              # [s⁻²] horizontal buoyancy gradient
Δy = 200kilometers      # width of the region of the front
Δb = Δy * M²           # buoyancy jump associated with the front
ϵᵇ = 1e-2 * Δb         # noise amplitude
f = -1e-4
τ = 30days
Δt₀ = 10minutes
stop_time = 10years
wizard_iteration_interval = 20
progress_iteration_interval = 100
filename = "mixed_layer_baroclinic_equilibrium"
mixed_layer_depth = h = 100meters
save_fields_interval = 10days

parameters = (; N², Δb, Δy, h, τ)

@inline unstratified_b(x, y, z, p) = p.Δb * ramp(y, p.Δy)
@inline linearly_stratified_b(x, y, z, p) = p.N² * z + p.Δb * ramp(y, p.Δy)

@inline function piecewise_stratified_b(x, y, z, p)
    by = p.Δb * ramp(y, p.Δy)
    bz = ifelse(z < -p.h, p.N² * (z + p.h), zero(z))
    return by + bz
end

bᵢ(x, y, z) = linearly_stratified_b(x, y, z, parameters) + ϵᵇ * randn()

grid = RectilinearGrid(architecture;
                       size = (Nx, Ny, Nz),
                       x = (0, Lx),
                       y = (-Ly/2, Ly/2),
                       z = (-Lz, 0),
                       topology = (Periodic, Bounded, Bounded))

@inline b_relaxation(x, y, z, t, b, p) = (linearly_stratified_b(x, y, z, p) - b) / p.τ 
b_forcing = Forcing(b_relaxation; parameters, field_dependencies=:b)

model = HydrostaticFreeSurfaceModel(; grid,
                                    coriolis = FPlane(; f),
                                    buoyancy = BuoyancyTracer(),
                                    tracers = :b,
                                    forcing = (; b=b_forcing),
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO())

set!(model, b=bᵢ)

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

b = model.tracers.b
u, v, w = model.velocities
ζ = ∂x(v) - ∂y(u)

B = Field(Average(b, dims=1))
U = Field(Average(u, dims=1))

slicers = (west = (1, :, :),
           east = (grid.Nx, :, :),
           south = (:, 1, :),
           north = (:, grid.Ny, :),
           bottom = (:, :, 1),
           top = (:, :, grid.Nz))

for side in keys(slicers)
    indices = slicers[side]

    simulation.output_writers[side] = JLD2OutputWriter(model, (; b, ζ);
                                                       filename = filename * "_$(side)_slice",
                                                       schedule = TimeInterval(save_fields_interval),
                                                       overwrite_existing = true,
                                                       indices)
end

simulation.output_writers[:zonal] = JLD2OutputWriter(model, (b=B, u=U);
                                                     schedule = TimeInterval(save_fields_interval),
                                                     overwrite_existing = true,
                                                     filename = filename * "_zonal_average")

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)

#fig = Figure(resolution = (1200, 800))
#ax = Axis(fig[1, 1])
#heatmap!(ax, interior(b, :, :, grid.Nz))


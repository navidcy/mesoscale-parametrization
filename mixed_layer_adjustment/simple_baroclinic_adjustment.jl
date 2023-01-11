using Oceananigans
using Oceananigans.Units
using Printf

architecture = GPU()
Nz = 32
Ny = 256
Nx = 1024
Lx = 8000kilometers     # east-west extent [m]
Ly = 2000kilometers     # north-south extent [m]
Lz = 1000meters         # depth [m]
N² = 1e-5               # [s⁻²] buoyancy frequency / stratification
M² = 2e-7               # [s⁻²] horizontal buoyancy gradient
Δy = 100kilometers      # width of the region of the front
Δb = Δy * M²            # buoyancy jump associated with the front
ϵᵇ = 1e-2 * Δb          # noise amplitude
Cᴰ = 2e-3
coriolis = BetaPlane(latitude=-45)
f = coriolis.f₀
stop_time = 30days
filename = "simple_baroclinic_adjustment"
save_interval = 6hours

ramp(y, Δy) = (1 + tanh(y / Δy)) / 2
d_ramp_dy(y, Δy) = sech(y / Δy)^2 / (2Δy)

bᵢ(x, y, z) =   Δb * ramp(y, Δy) + N² * z + ϵᵇ * (2rand() - 1)
uᵢ(x, y, z) = - Δb / f * d_ramp_dy(y, Δy) * (z + Lz/2)

grid = RectilinearGrid(architecture;
                       size = (Nx, Ny, Nz),
                       x = (0, Lx),
                       y = (-Ly/2, Ly/2),
                       z = (-Lz, 0),
                       topology = (Periodic, Bounded, Bounded))

@inline ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2
@inline speedᶠᶜᶜ(i, j, k, grid, u, v) = @inbounds sqrt(u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², v))
@inline speedᶜᶠᶜ(i, j, k, grid, u, v) = @inbounds sqrt(v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², u))
@inline u_drag(i, j, grid, clock, fields, Cᴰ) = @inbounds - Cᴰ * u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields.u, fields.v)
@inline v_drag(i, j, grid, clock, fields, Cᴰ) = @inbounds - Cᴰ * v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields.u, fields.v)

u_drag_bc = FluxBoundaryCondition(u_drag; discrete_form=true, parameters=Cᴰ)
v_drag_bc = FluxBoundaryCondition(v_drag; discrete_form=true, parameters=Cᴰ)

#@inline u_drag(x, y, t, u, v, Cᴰ) = - Cᴰ * sqrt(u^2 + v^2) * u
#@inline v_drag(x, y, t, u, v, Cᴰ) = - Cᴰ * sqrt(u^2 + v^2) * v
#
#u_drag_bc = FluxBoundaryCondition(u_drag; field_dependencies=(:u, :v), parameters=Cᴰ)
#v_drag_bc = FluxBoundaryCondition(v_drag; field_dependencies=(:u, :v), parameters=Cᴰ)

u_bcs = FieldBoundaryConditions(bottom=u_drag_bc)
v_bcs = FieldBoundaryConditions(bottom=v_drag_bc)

boundary_conditions=(u=u_bcs, v=v_bcs)

model = HydrostaticFreeSurfaceModel(; grid, boundary_conditions, coriolis,
                                    buoyancy = BuoyancyTracer(),
                                    tracers = (:b, :r),
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO())

rᵢ(x, y, z) = z
set!(model, b=bᵢ, u=uᵢ, r=rᵢ)

simulation = Simulation(model; Δt=10minutes, stop_time)
wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=20minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(20))

b, r = model.tracers
u, v, w = model.velocities
ζ = ∂x(v) - ∂y(u)

B = Field(Average(b, dims=1))
R = Field(Average(r, dims=1))
U = Field(Average(u, dims=1))
Uz = Field(Average(∂z(u), dims=1))
N² = ∂z(b)
N²avg = Field(Average(N², dims=1))

u′ = u - U
k = @at (Center, Center, Center) (u′^2 + v^2) / 2
K = Field(Average(k, dims=1))

slicers = (west = (1, :, :),
           east = (grid.Nx, :, :),
           south = (:, 1, :),
           north = (:, grid.Ny, :),
           bottom = (:, :, 1),
           top = (:, :, grid.Nz))

for side in keys(slicers)
    indices = slicers[side]

    simulation.output_writers[side] = JLD2OutputWriter(model, (; b, ζ, k, N²);
                                                       filename = filename * "_$(side)_slice",
                                                       schedule = TimeInterval(save_interval),
                                                       overwrite_existing = true,
                                                       indices)
end

simulation.output_writers[:zonal] = JLD2OutputWriter(model, (b=B, u=U, k=K, uz=Uz, N²=N²avg);
                                                     schedule = TimeInterval(save_interval),
                                                     overwrite_existing = true,
                                                     filename = filename * "_zonal_average")

wall_clock = Ref(time_ns())

function print_progress(sim)

    compute!(U)
    compute!(K)

    msg1 = @sprintf("i: % 4d, t: % 12s, wall time: % 12s, extrema(U): (%6.3e, %6.3e), max(K): %6.3e, ",
                    iteration(sim),
                    prettytime(sim),
                    prettytime(1e-9 * (time_ns() - wall_clock[])),
                    maximum(U),
                    minimum(U),
                    maximum(K))

    msg2 = @sprintf("max(u): (%6.3e, %6.3e, %6.3e) m/s, next Δt: %s",
                    maximum(abs, sim.model.velocities.u),
                    maximum(abs, sim.model.velocities.v),
                    maximum(abs, sim.model.velocities.w),
                    prettytime(sim.Δt))
                    
    @info msg1 * msg2

    wall_clock[] = time_ns()
    
    return nothing
end

simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(100))

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)

#fig = Figure(resolution = (1200, 800))
#ax = Axis(fig[1, 1])
#heatmap!(ax, interior(b, :, :, grid.Nz))


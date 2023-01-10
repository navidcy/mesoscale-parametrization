using Printf
using Statistics
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: xnode, ynode, znode

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    MixingLength,
    TurbulentKineticEnergyEquation,
    convective_mixing_lengthᶜᶜᶠ,
    TKE_mixing_lengthᶜᶜᶠ,
    momentum_mixing_lengthᶜᶜᶠ,
    tracer_mixing_lengthᶜᶜᶠ,
    Kuᶜᶜᶠ,
    Kcᶜᶜᶠ,
    Keᶜᶜᶠ

using Oceananigans.TurbulenceClosures: RiBasedVerticalDiffusivity
using Random

using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ

Random.seed!(123)

arch = GPU()
fileprefix = "eddying_channel_catke"
save_interval = 2days
stop_time = 20years

# Horizontal domain
Lx = 2000kilometers # zonal domain length [m]
Ly = 2000kilometers # meridional domain length [m]

# number of grid points
Nx = 128
Ny = 128

# Vertical grid generation
# Fixed spacing in the upper ocean
Δz₀ = 5.0         # surface layer grid spacing
h₀ = 100          # surface layer extent

# Generate surface layer grid
z = [-Δz₀ * (k-1) for k = 1:ceil(h₀ / Δz₀)]

# Generate stretched interior grid
γ = 1.03 # stretching parameter
Lz₀ = 2kilometers # minimum domain depth

while z[end] > - Lz₀
    Δz = (z[end-1] - z[end])^γ
    push!(z, round(z[end] - Δz, digits=1))
end

# Reverse grid to be right-side-up
z = reverse(z)

# Infer domain parameters
Lz = z[1]
Nz = length(z) - 1

grid = RectilinearGrid(arch;
                       topology = (Periodic, Bounded, Bounded),
                       size = (Nx, Ny, Nz),
                       halo = (4, 4, 4),
                       x = (0, Lx),
                       y = (0, Ly),
                       z)

@info "Built a grid: $grid."

parameters = (Ly = Ly,
              Qᵇ = + 1e-8,
              Qᵘ = - 2e-4,
              Cᴰ = 2e-3)

# Buoyancy fluxes
@inline function buoyancy_flux(i, j, grid, clock, model_fields, p)
    y = ynode(Center(), j, grid)
    return ifelse(y > 5 * p.Ly / 6, zero(grid), p.Qᵇ * cos(3π * y / p.Ly))
end

buoyancy_flux_bc = FluxBoundaryCondition(buoyancy_flux, discrete_form=true, parameters=parameters)

# Wind stress
@inline function u_stress(i, j, grid, clock, model_fields, p)
    y = ynode(Center(), j, grid)
    return p.Qᵘ * sin(π * y / p.Ly)
end

u_stress_bc = FluxBoundaryCondition(u_stress, discrete_form=true, parameters=parameters)

# Bottom drag
@inline ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2 
@inline speedᶠᶜᶜ(i, j, k, grid, u, v) = @inbounds sqrt(u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², v))
@inline speedᶜᶠᶜ(i, j, k, grid, u, v) = @inbounds sqrt(v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², u))

@inline u_drag(i, j, grid, clock, f, p) = @inbounds - p.Cᴰ * f.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, f.u, f.v)
@inline v_drag(i, j, grid, clock, f, p) = @inbounds - p.Cᴰ * f.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, f.u, f.v)

u_drag_bc = FluxBoundaryCondition(u_drag, discrete_form=true, parameters=parameters)
v_drag_bc = FluxBoundaryCondition(v_drag, discrete_form=true, parameters=parameters)

b_bcs = FieldBoundaryConditions(top = buoyancy_flux_bc)
u_bcs = FieldBoundaryConditions(top = u_stress_bc, bottom = u_drag_bc)
v_bcs = FieldBoundaryConditions(bottom = v_drag_bc)

#####
##### Coriolis
#####

#vertical_mixing_closure = RiBasedVerticalDiffusivity()
vertical_mixing_closure = CATKEVerticalDiffusivity()

#####
##### Model building
#####

@info "Building a model..."

model = HydrostaticFreeSurfaceModel(; grid,
                                    free_surface = ImplicitFreeSurface(),
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO(),
                                    buoyancy = BuoyancyTracer(),
                                    coriolis = BetaPlane(latitude=-45),
                                    closure = vertical_mixing_closure,
                                    boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs),
                                    tracers = (:b, :e))

@info "Built $model."

#####
##### Baroclinically unstable initial condition
#####

Δy = 100kilometers
N² = 1e-5              # [s⁻²] buoyancy frequency / stratification
M² = 2e-7              # [s⁻²] horizontal buoyancy gradient
Δb = Δy * M²           # buoyancy jump associated with the front
f = model.coriolis.f₀
ϵᵘ = 1e-6 * f * Δb / Δy
shear_layer_depth = H = 500meters

ramp(y, Δy) = (1 + tanh(y / Δy)) / 2
d_ramp_dy(y, Δy) = sech(y / Δy)^2 / 2Δy

function piecewise_stratified_b(x, y, z)
    by = Δb * ramp(y, Δy) * exp(z / H)
    bz = N² * z
    return by + bz
end

bᵢ(x, y, z) = piecewise_stratified_b(x, y-Ly/2, z)
uᵢ(x, y, z) = - Δb / f * d_ramp_dy(y-Ly/2, Δy) * H * exp(z / H) + ϵᵘ * (2rand() - 1)
vᵢ(x, y, z) = ϵᵘ * (2rand() - 1)

set!(model, b=bᵢ, u=uᵢ, v=vᵢ, e=1e-6)

#####
##### Simulation building
#####

simulation = Simulation(model; Δt=20minutes, stop_time)

# add timestep wizard callback
wizard = TimeStepWizard(cfl=0.1, max_change=1.1, max_Δt=20minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# add progress callback
wall_clock = Ref(time_ns())

function print_progress(sim)
    msg1 = @sprintf("i: % 8d, t: % 16s, wall time: % 16s, extrema(e): (%6.3e, %6.3e) m² s⁻², ",
                    sim.model.clock.iteration,
                    prettytime(sim.model.clock.time),
                    prettytime(1e-9 * (time_ns() - wall_clock[])),
                    maximum(sim.model.tracers.e),
                    minimum(sim.model.tracers.e))

    msg2 = @sprintf("max|u|: (%6.3e, %6.3e, %6.3e) m s⁻¹, next Δt: %s",
                    maximum(abs, sim.model.velocities.u),
                    maximum(abs, sim.model.velocities.v),
                    maximum(abs, sim.model.velocities.w),
                    prettytime(sim.Δt))

    @info msg1 * msg2

    wall_clock[] = time_ns()

    return nothing
end

simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(100))

#####
##### Diagnostics
#####

u, v, w = model.velocities
b = model.tracers.b

U = Field(Average(u, dims=1))
B = Field(Average(b, dims=1))
u′ = u - U
k = sqrt(u′^2 + v^2)
K = Average(k, dims=1)

ζ = ∂x(v) - ∂y(u)
N² = ∂z(b)

#=
closure = vertical_mixing_closure
grid = model.grid
velocities = model.velocities
tracers = model.tracers
buoyancy = model.buoyancy
clock = model.clock
top_tracer_bcs = NamedTuple(c => tracers[c].boundary_conditions.top for c in propertynames(tracers))

computed_dependencies = (closure, velocities, tracers, buoyancy, clock, top_tracer_bcs)
ℓe = KernelFunctionOperation{Center, Center, Face}(TKE_mixing_lengthᶜᶜᶠ, grid; computed_dependencies)
ℓc = KernelFunctionOperation{Center, Center, Face}(tracer_mixing_lengthᶜᶜᶠ, grid; computed_dependencies)
ℓu = KernelFunctionOperation{Center, Center, Face}(momentum_mixing_lengthᶜᶜᶠ, grid; computed_dependencies)
κe = KernelFunctionOperation{Center, Center, Face}(Keᶜᶜᶠ, grid; computed_dependencies)
κc = KernelFunctionOperation{Center, Center, Face}(Kcᶜᶜᶠ, grid; computed_dependencies)
κu = KernelFunctionOperation{Center, Center, Face}(Kuᶜᶜᶠ, grid; computed_dependencies)
outputs = merge(model.velocities, model.tracers, (; ζ, N², k, κe, ℓe, κc, ℓc, κu, ℓu))
=#

outputs = merge(model.velocities, model.tracers, (; ζ, k, N²))
averaged_outputs = (u=U, b=B, k=K)

#####
##### Build checkpointer and output writer
#####

simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                        schedule = TimeInterval(10years),
                                                        prefix = fileprefix,
                                                        overwrite_existing = true)

slicers = (west = (1, :, :),
           east = (grid.Nx, :, :),
           south = (:, 1, :),
           north = (:, grid.Ny, :),
           bottom = (:, :, 1),
           top = (:, :, grid.Nz))

for side in keys(slicers)
    indices = slicers[side]

    simulation.output_writers[side] = JLD2OutputWriter(model, outputs;
                                                       schedule = TimeInterval(save_interval),
                                                       indices,
                                                       filename = fileprefix * "_$(side)_slice",
                                                       overwrite_existing = true)
end

simulation.output_writers[:zonal] = JLD2OutputWriter(model, averaged_outputs;
                                                     schedule = TimeInterval(save_interval),
                                                     filename = fileprefix * "_zonal_average",
                                                     overwrite_existing = true)

@info "Running the simulation..."

#using Oceananigans.Simulations: NaNChecker
#simulation.callbacks[:nan_checker] = Callback(NaNChecker(; fields=(; e=model.tracers.e)))

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)


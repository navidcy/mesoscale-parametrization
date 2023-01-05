using Printf
using Statistics
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: xnode, ynode, znode
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using Oceananigans.TurbulenceClosures: RiBasedVerticalDiffusivity
using Random

Random.seed!(123)
arch = GPU()
fileprefix = "eddying_channel_catke"
saveinterval = 2days

# Domain
Lx = 2000kilometers # zonal domain length [m]
Ly = 1000kilometers # meridional domain length [m]
Lz = 2kilometers    # depth [m]

# number of grid points
Nx = 256
Ny = 128
Nz = 64

# Stretched grid implementations
σ = 1.04 # linear stretching factor
Δz_center_linear(k) = Lz * (σ - 1) * σ^(Nz - k) / (σ^Nz - 1) # k=1 is the bottom-most cell, k=Nz is the top cell
linearly_spaced_faces(k) = k==1 ? -Lz : - Lz + sum(Δz_center_linear.(1:k-1))

refinement = 2 # controls spacing near surface (higher means finer spaced)
stretching = 4  # controls rate of stretching at bottom

h(k) = (k - 1) / Nz
ζ₀(k) = 1 + (h(k) - 1) / refinement
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(arch;
                       topology = (Periodic, Bounded, Bounded),
                       size = (Nx, Ny, Nz),
                       halo = (4, 4, 4),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = (-Lz, 0)) # z_faces)

@info "Built a grid: $grid."

parameters = (Ly = Ly,
              Qᵇ = 1e-8,
              Qᵘ = - 1e-4,
              Cᴰ = 2e-3)

# Buoyancy fluxes
@inline function buoyancy_flux(i, j, grid, clock, model_fields, p)
    y = ynode(Center(), j, grid)
    return p.Qᵇ * cos(3π * y / p.Ly)
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

@inline u_drag(i, j, grid, clock, f, p) = @inbounds - p.Cᴰ * f.u[i, j, k] * speedᶠᶜᶜ(i, j, 1, grid, f.u, f.v)
@inline v_drag(i, j, grid, clock, f, p) = @inbounds - p.Cᴰ * f.v[i, j, k] * speedᶜᶠᶜ(i, j, 1, grid, f.u, f.v)

u_drag_bc = FluxBoundaryCondition(u_drag, discrete_form=true, parameters=parameters)
v_drag_bc = FluxBoundaryCondition(v_drag, discrete_form=true, parameters=parameters)

b_bcs = FieldBoundaryConditions(top = buoyancy_flux_bc)

u_bcs = FieldBoundaryConditions(top = u_stress_bc, bottom = u_drag_bc)
v_bcs = FieldBoundaryConditions(bottom = v_drag_bc)

#####
##### Coriolis
#####

#vertical_mixing_closure = CATKEVerticalDiffusivity()
vertical_mixing_closure = RiBasedVerticalDiffusivity()

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
                                    tracers = (:b, :e),
                                    boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs))

@info "Built $model."

#####
##### Baroclinically unstable initial condition
#####

Δy = 100kilometers
N² = 1e-5              # [s⁻²] buoyancy frequency / stratification
M² = 2e-7              # [s⁻²] horizontal buoyancy gradient
Δb = Δy * M²           # buoyancy jump associated with the front
ϵᵇ = 1e-2 * Δb         # noise amplitude
f = coriolis.f₀
mixed_layer_depth = h = 200meters
shear_layer_depth = H = 500meters

ramp(y, Δy) = (1 + tanh(y / Δy)) / 2
d_ramp_dy(y, Δy) = sech(y / Δy)^2 / (2Δy)

function piecewise_stratified_b(x, y, z)
    by = Δb * ramp(y, Δy) * exp(z / H)
    bz = ifelse(z < -h, N² * (z + h), zero(z))
    return by + bz
end

bᵢ(x, y, z) = piecewise_stratified_b(x, y, z) + ϵᵇ * randn()
uᵢ(x, y, z) = - Δb / f * d_ramp_dy(y, Δy) * H * exp(z / H)

set!(model, b=bᵢ, u=uᵢ)

#####
##### Simulation building
#####

simulation = Simulation(model, Δt=1minute, stop_time=10days)

# add timestep wizard callback
wizard = TimeStepWizard(cfl=0.1, max_change=1.1, max_Δt=20minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

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
b = model.tracers.b

U = Field(Average(u, dims=1))
B = Field(Average(b, dims=1))
u′ = u - U
k = sqrt(u′^2 + v^2)
K = Average(k, dims=1)

ζ = Field(∂x(v) - ∂y(u))

outputs = merge(model.velocities, model.tracers, (; ζ, k))
averaged_outputs = (u=U, b=B, k=K)

#####
##### Build checkpointer and output writer
#####

simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                        schedule = TimeInterval(10years),
                                                        prefix = fileprefix,
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
                                                       schedule = TimeInterval(saveinterval),
                                                       indices,
                                                       filename = fileprefix * "_$(side)_slice",
                                                       overwrite_existing = true)
end

simulation.output_writers[:zonal] = JLD2OutputWriter(model, averaged_outputs;
                                                     schedule = TimeInterval(saveinterval),
                                                     filename = fileprefix * "_zonal_average",
                                                     overwrite_existing = true)

@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)


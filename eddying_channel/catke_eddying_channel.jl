using Printf
using Statistics
using JLD2

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: xnode, ynode, znode
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity

using Random
Random.seed!(1234)

arch = GPU()
filename = "catke_eddying_channel"

# Domain
const Lx = 2000kilometers # zonal domain length [m]
const Ly = 2000kilometers # meridional domain length [m]
const Lz = 2kilometers    # depth [m]

# number of grid points
Nx = 128
Ny = 128
Nz = 48

save_fields_interval = 7days
stop_time = 100years
Δt₀ = 10minutes

# stretched grid

# we implement here a linearly streched grid in which the top grid cell has Δzₜₒₚ
# and every other cell is bigger by a factor σ, e.g.,
# Δzₜₒₚ, Δzₜₒₚ * σ, Δzₜₒₚ * σ², ..., Δzₜₒₚ * σᴺᶻ⁻¹,
# so that the sum of all cell heights is Lz

# Given Lz and stretching factor σ > 1 the top cell height is Δz_top = Lz * (σ - 1) / σ^(Nz - 1)

σ = 1.05 # linear stretching factor
Δz_center_linear(k) = Lz * (σ - 1) * σ^(Nz - k) / (σ^Nz - 1) # k=1 is the bottom-most cell, k=Nz is the top cell
linearly_spaced_faces(k) = k == 1 ? -Lz : - Lz + sum(Δz_center_linear.(1:k-1))

grid = RectilinearGrid(arch;
                       topology = (Periodic, Bounded, Bounded),
                       size = (Nx, Ny, Nz),
                       halo = (3, 3, 3),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = linearly_spaced_faces)

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
              λt = 7days)                          # relaxation time scale [s]

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

Fb = Forcing(buoyancy_relaxation, discrete_form=true, parameters=parameters)

# Turbulence closures

κh = 0.5e-5 # [m²/s] horizontal diffusivity
νh = 30.0   # [m²/s] horizontal viscosity

Δh = grid.Δxᶜᵃᵃ
κ₄h = 1 / 20days * Δh^4 

κz = 0.5e-5 # [m²/s] vertical diffusivity
νz = 3e-4   # [m²/s] vertical viscocity

vertical_diffusivity = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν = νz, κ = κz)
horizontal_diffusivity = HorizontalScalarDiffusivity(ν=νh, κ=κh)
#horizontal_diffusivity = HorizontalScalarBiharmonicDiffusivity(ν=κ₄h, κ=κ₄h)
catke = CATKEVerticalDiffusivity()
convective_adjustment = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = 1.0,
                                                                convective_νz = 0.0)

closures = (catke, vertical_diffusivity, horizontal_diffusivity)
#closures = (convective_adjustment, vertical_diffusivity, horizontal_diffusivity)

#####
##### Model building
#####

@info "Building a model..."

model = HydrostaticFreeSurfaceModel(grid = grid,
                                    free_surface = ImplicitFreeSurface(),
                                    momentum_advection = WENO5(; grid),
                                    tracer_advection = WENO5(; grid),
                                    buoyancy = BuoyancyTracer(),
                                    coriolis = BetaPlane(f₀ = -1e-4, β = 1e-11),
                                    closure = closures,
                                    tracers = (:b, :c, :e),
                                    boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs),
                                    forcing = (; b=Fb))

@info "Built $model."

#####
##### Initial conditions
#####

# resting initial condition
ε(σ) = σ * randn()
bᵢ(x, y, z) = parameters.ΔB * (exp(z / parameters.h) - exp(-Lz / parameters.h)) / (1 - exp(-Lz / parameters.h)) + ε(1e-8)
uᵢ(x, y, z) = ε(1e-6)
vᵢ(x, y, z) = ε(1e-6)
wᵢ(x, y, z) = ε(1e-6)

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
wizard = TimeStepWizard(cfl=0.1, max_change=1.1, max_Δt=10minutes)
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

simulation.callbacks[:print_progress] = Callback(print_progress, IterationInterval(100))

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

u′v′ = Average(uv_op, dims=1)
v′w′ = Average(vw_op, dims=1)
u′w′ = Average(uw_op, dims=1)
b′b′ = Average(b′ * b′, dims=1)
v′b′ = Average(b′ * v′, dims=1)
w′b′ = Average(b′ * w′, dims=1)
c′c′ = Average(c′ * c′, dims=1)
v′c′ = Average(c′ * v′, dims=1)
w′c′ = Average(c′ * w′, dims=1)

#=
f₀ = model.coriolis.f₀
const β = model.coriolis.β
βy(x, y, z) = β * y
N² = ∂z(B)
q = @at (Center, Center, Center) ζ / f₀ + ∂z(b / N²) + βy
Q = Field(Average(q, dims=1))
q′ = q - Q
q′q′ = Average(q′ * q′, dims=1)
=#

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

#=
simulation.output_writers[:zonal] = JLD2OutputWriter(model, zonally_averaged_outputs;
                                                     schedule = TimeInterval(save_fields_interval),
                                                     prefix = filename * "_zonal_average",
                                                     force = true)
=#

@info "Running the simulation..."

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)

using GLMakie

ζ_ts = FieldTimeSeries(filename * "_top_slice.jld2", "ζ", architecture=CPU())
Nt = length(ζ_ts.times)

fig = Figure()
ax = Axis(fig[1, 1])
heatmap!(ax, interior(ζ_ts[Nt])[:, :, 1])

display(fig)
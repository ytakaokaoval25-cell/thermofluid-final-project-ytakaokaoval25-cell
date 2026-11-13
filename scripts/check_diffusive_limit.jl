# 流れなしで Nu = 1 になることの確認（issue #5）
#
# 浮力を切った（buoyancy = nothing，advection = nothing）モデルで純拡散だけを解き，
# 純拡散解 b_diff から出発して離散的な定常状態まで進める．その状態の Nu
# （公式例と同じ演算子 / 数値の χ_diff）が 1 に近いことを確かめる．
# 拡散の定常解は κ によらないので，時間を短くするため κ = 1 で解く．
# 格子と境界条件は run_horizontal_convection.jl と同じ．
#
# 使い方：julia --project=. scripts/check_diffusive_limit.jl
#   HC_NX, HC_NZ（既定 128, 64），HC_DIFF_STOP_TIME（既定 1.5），HC_OUTPUT（既定 output）
# 出力：<HC_OUTPUT>/diffusive_limit_<Nx>x<Nz>.toml

using Oceananigans
using Printf

include(joinpath(@__DIR__, "nusselt_reference.jl"))

envnum(T, key, default) = parse(T, get(ENV, key, string(default)))
Nx = envnum(Int, "HC_NX", 128)
Nz = envnum(Int, "HC_NZ", 64)
stop_time = envnum(Float64, "HC_DIFF_STOP_TIME", 1.5)
outdir = get(ENV, "HC_OUTPUT", "output")
mkpath(outdir)

H, Lx, b★ = 1.0, 2.0, 1.0
κ = 1.0

grid = RectilinearGrid(size = (Nx, Nz), x = (-Lx/2, Lx/2), z = (-H, 0),
                       topology = (Bounded, Flat, Bounded))

b_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(bˢ, parameters=(; b★, Lx)))
model = NonhydrostaticModel(grid; advection = nothing, buoyancy = nothing, tracers = :b,
                            timestepper = :RungeKutta3,
                            closure = ScalarDiffusivity(; ν = 0, κ),
                            boundary_conditions = (; b = b_bcs))

b_diff = diffusive_buoyancy(grid; b★, Lx, H)
χ_diff = mean_dissipation(b_diff, κ; Lx, H)
set!(model, b = (x, z) -> -b★ * cos(2π * x / Lx) * cosh(2π * (H + z) / Lx) / cosh(2π * H / Lx))

Δt = 0.1 * min(Lx / Nx, H / Nz)^2 / κ
simulation = Simulation(model; Δt, stop_time)
Nu_of(m) = mean_dissipation(m.tracers.b, κ; Lx, H) / χ_diff

history = Tuple{Float64, Float64}[(0.0, Nu_of(model))]
simulation.callbacks[:nu] = Callback(s -> push!(history, (time(s), Nu_of(s.model))), TimeInterval(0.25))
wall = @elapsed run!(simulation)
push!(history, (time(simulation), Nu_of(model)))

b = interior(model.tracers.b)
max_change = maximum(abs, b .- interior(b_diff))

open(joinpath(outdir, @sprintf("diffusive_limit_%dx%d.toml", Nx, Nz)), "w") do io
    println(io, "# 流れなし（純拡散）での Nu の確認（scripts/check_diffusive_limit.jl が生成）")
    println(io, "Nx = ", Nx, "\nNz = ", Nz, "\nkappa = ", κ, "\nstop_time = ", stop_time, "\ndt = ", Δt)
    println(io, "iterations = ", iteration(simulation), "\nwall_time_seconds = ", wall)
    println(io, "chi_diff_numeric_over_analytic = ", χ_diff / analytic_chi_diff(κ; b★, Lx, H))
    println(io, "Nu_final = ", last(history)[2], "  # 離散定常状態の Nu（1 に近いべき）")
    println(io, "max_abs_b_change_from_b_diff = ", max_change)
    println(io, "t_Nu = [", join((@sprintf("[%.2f, %.8f]", h...) for h in history), ", "), "]")
end
println(read(joinpath(outdir, @sprintf("diffusive_limit_%dx%d.toml", Nx, Nz)), String))

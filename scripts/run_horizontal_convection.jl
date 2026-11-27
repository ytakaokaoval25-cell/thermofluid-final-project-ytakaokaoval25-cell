# 水平対流の計算（公式例の再現用）
#
# 出典：Oceananigans.jl v0.110.15 examples/horizontal_convection.jl（MIT License）
#       無改変のコピーは upstream/horizontal_convection.jl
#
# 公式例からの変更点（物理・数値設定は変えていない）：
#   1. 格子 (Nx, Nz)，終了時刻，Ra，出力先を環境変数で切り替えられるようにした
#      （既定値は公式例と同じ：128×64，t = 40，Ra = 1e8）
#   2. シミュレーション部分と作図部分を分けた（作図は scripts/plot_horizontal_convection.jl）
#   3. 実行条件と実行時間を output/<名前>_run.toml に記録するようにした
#   4. 検証（離散発散）のため，保存する量に速度 u, w を追加した
#   5. （issue #5）時間刻みの拡散 CFL 上限を環境変数 HC_DIFFUSIVE_CFL で指定できるようにした．
#      既定は Inf（公式例と同じく拡散 CFL で制限しない）．Ra = 1e6 では κ が大きく，
#      max_Δt = 0.1 に達する前に拡散の安定限界（RK3）を超えて NaN になるため，そのケースだけ使う．
#
# 使い方（リポジトリのルートで）：
#   julia --project=. scripts/run_horizontal_convection.jl
#   HC_NX=64 HC_NZ=32 HC_STOP_TIME=2 HC_NAME=smoke julia --project=. scripts/run_horizontal_convection.jl

using Oceananigans
using Printf
using Dates

envnum(T, key, default) = parse(T, get(ENV, key, string(default)))

Nx = envnum(Int, "HC_NX", 128)
Nz = envnum(Int, "HC_NZ", 64)
stop_time = envnum(Float64, "HC_STOP_TIME", 40.0)
Ra = envnum(Float64, "HC_RA", 1e8)
name = get(ENV, "HC_NAME", "horizontal_convection")
diffusive_cfl = envnum(Float64, "HC_DIFFUSIVE_CFL", Inf)
outdir = get(ENV, "HC_OUTPUT", "output")
mkpath(outdir)

# ### The grid（公式例と同じ）

H = 1            # vertical domain extent
Lx = 2H          # horizontal domain extent

grid = RectilinearGrid(size = (Nx, Nz),
                       x = (-Lx/2, Lx/2),
                       z = (-H, 0),
                       topology = (Bounded, Flat, Bounded))

# ### Boundary conditions（公式例と同じ）

b★ = 1
@inline bˢ(x, t, p) = - p.b★ * cos(2π * x / p.Lx)

b_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(bˢ, parameters=(; b★, Lx)))

# ### Non-dimensional control parameters（Ra だけ環境変数で変えられる）

Pr = 1      # Prandtl number

ν = sqrt(Pr * b★ * Lx^3 / Ra)  # Laplacian viscosity
κ = ν / Pr                     # Laplacian diffusivity

# ## Model instantiation（公式例と同じ）

model = NonhydrostaticModel(grid;
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            closure = ScalarDiffusivity(; ν, κ),
                            boundary_conditions = (; b=b_bcs))

# ## Simulation set-up（公式例と同じ）

simulation = Simulation(model, Δt=1e-2, stop_time=stop_time)

conjure_time_step_wizard!(simulation, IterationInterval(50), cfl=0.7, max_Δt=1e-1, diffusive_cfl=diffusive_cfl)

progress(sim) = @printf("Iter: % 6d, sim time: % 1.3f, wall time: % 10s, Δt: % 1.4f, advective CFL: %.2e, diffusive CFL: %.2e\n",
                        iteration(sim), time(sim), prettytime(sim.run_wall_time),
                        sim.Δt, AdvectiveCFL(sim.Δt)(sim.model), DiffusiveCFL(sim.Δt)(sim.model))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

# ### Output（公式例の s, b, ζ に検証用の u, w を追加．保存先は output/）

u, v, w = model.velocities
b = model.tracers.b

s = @at (Center, Center, Center) sqrt(u^2 + w^2)
ζ = ∂z(u) - ∂x(w)

saved_output_filename = joinpath(outdir, name * ".jld2")

simulation.output_writers[:fields] = JLD2Writer(model, (; s, b, ζ, u, w),
                                                schedule = TimeInterval(0.5),
                                                filename = saved_output_filename,
                                                with_halos = true,
                                                overwrite_existing = true)

Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)

wall = @elapsed run!(simulation)

# ## 実行条件と時間の記録（追加）

open(joinpath(outdir, name * "_run.toml"), "w") do io
    println(io, "date = \"", Dates.now(), "\"")
    println(io, "julia = \"", VERSION, "\"")
    println(io, "oceananigans = \"", pkgversion(Oceananigans), "\"")
    println(io, "threads = ", Threads.nthreads())
    println(io, "Nx = ", Nx)
    println(io, "Nz = ", Nz)
    println(io, "Ra = ", Ra)
    println(io, "Pr = ", Pr)
    println(io, "stop_time = ", stop_time)
    println(io, "diffusive_cfl = ", isinf(diffusive_cfl) ? "\"Inf\"" : diffusive_cfl, "  # 時間刻みの拡散 CFL 上限（公式例は Inf）")
    println(io, "iterations = ", iteration(simulation))
    println(io, "wall_time_seconds = ", wall, "  # run! 全体（初回はコンパイル時間を含む）")
    println(io, "output = \"", saved_output_filename, "\"")
end

@info "done" saved_output_filename wall

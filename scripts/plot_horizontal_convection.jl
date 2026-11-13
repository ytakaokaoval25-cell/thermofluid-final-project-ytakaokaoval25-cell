# 水平対流の作図（公式例の後半部分）
#
# 出典：Oceananigans.jl v0.110.15 examples/horizontal_convection.jl（MIT License）の
#       "Load saved output, process, visualize" 以降．
# 変更点：入力・出力のパスを変えた．最終時刻の場の図と Nu・K の時系列図を png で保存するようにした．
#
# 使い方：julia --project=. scripts/plot_horizontal_convection.jl

using CairoMakie
using Oceananigans
using Oceananigans.Fields
using Printf

name = get(ENV, "HC_NAME", "horizontal_convection")
outdir = get(ENV, "HC_OUTPUT", "output")
figdir = get(ENV, "HC_FIGURES", "results")
mkpath(figdir)
Ra = parse(Float64, get(ENV, "HC_RA", "1e8"))

H, Lx, b★, Pr = 1, 2, 1, 1
ν = sqrt(Pr * b★ * Lx^3 / Ra)
κ = ν / Pr

saved_output_filename = joinpath(outdir, name * ".jld2")

s_timeseries = FieldTimeSeries(saved_output_filename, "s")
b_timeseries = FieldTimeSeries(saved_output_filename, "b")
ζ_timeseries = FieldTimeSeries(saved_output_filename, "ζ")

times = b_timeseries.times

χ_timeseries = deepcopy(b_timeseries)

for n in 1:length(times)
    bn = b_timeseries[n]
    χ_timeseries[n] .= @at (Center, Center, Center) κ * (∂x(bn)^2 + ∂z(bn)^2)
end

n = Observable(length(times))

title = @lift @sprintf("t=%1.2f", times[$n])

sn = @lift s_timeseries[$n]
ζn = @lift ζ_timeseries[$n]
bn = @lift b_timeseries[$n]
χn = @lift χ_timeseries[$n]

slim = 0.6
blim = 0.6
ζlim = 9
χlim = 0.025

axis_kwargs = (xlabel = L"x / H",
               ylabel = L"z / H",
               limits = ((-Lx/2, Lx/2), (-H, 0)),
               aspect = Lx / H,
               titlesize = 20)

fig = Figure(size = (600, 1100))

ax_s = Axis(fig[2, 1]; title = L"speed, $(u^2+w^2)^{1/2} / (L_x b_*)^{1/2}$", axis_kwargs...)
ax_b = Axis(fig[3, 1]; title = L"buoyancy, $b / b_*$", axis_kwargs...)
ax_ζ = Axis(fig[4, 1]; axis_kwargs...,
            title = L"vorticity, $(∂u/∂z - ∂w/∂x) \, (L_x / b_*)^{1/2}$")
ax_χ = Axis(fig[5, 1]; axis_kwargs...,
            title = L"buoyancy dissipation, $κ |\mathbf{\nabla}b|^2 \, (L_x / {b_*}^5)^{1/2}$")

fig[1, :] = Label(fig, title, fontsize=24, tellwidth=false)

hm_s = heatmap!(ax_s, sn; colorrange=(0, slim), colormap=:speed)
Colorbar(fig[2, 2], hm_s)
hm_b = heatmap!(ax_b, bn; colorrange=(-blim, blim), colormap=:thermal)
Colorbar(fig[3, 2], hm_b)
hm_ζ = heatmap!(ax_ζ, ζn; colorrange=(-ζlim, ζlim), colormap=:balance)
Colorbar(fig[4, 2], hm_ζ)
hm_χ = heatmap!(ax_χ, χn; colorrange=(0, χlim), colormap=:dense)
Colorbar(fig[5, 2], hm_χ)

save(joinpath(figdir, name * "_fields_final.png"), fig)

frames = 1:length(times)
record(fig, joinpath(outdir, name * ".mp4"), frames, framerate=8) do i
    n[] = i
end

# ### Nusselt 数と運動エネルギー（公式例と同じ計算）

χ_diff = κ * b★^2 * π * tanh(2π * H / Lx) / (Lx * H)

t = b_timeseries.times
kinetic_energy, Nu = zeros(length(t)), zeros(length(t))

for n = 1:length(t)
    ke = ∫dV(1/2 * s_timeseries[n]^2 / (Lx * H))
    kinetic_energy[n] = ke[1, 1, 1]
    χ = ∫dV(χ_timeseries[n] / (Lx * H))
    Nu[n] = χ[1, 1, 1] / χ_diff
end

fig2 = Figure(size = (850, 450))
ax_KE = Axis(fig2[1, 1], xlabel = L"t \, (b_* / L_x)^{1/2}", ylabel = L"KE $ / (L_x b_*)$")
lines!(ax_KE, t, kinetic_energy; linewidth = 3)
ax_Nu = Axis(fig2[2, 1], xlabel = L"t \, (b_* / L_x)^{1/2}", ylabel = L"Nu")
lines!(ax_Nu, t, Nu; linewidth = 3)
save(joinpath(figdir, name * "_Nu_KE.png"), fig2)

open(joinpath(figdir, name * "_official_Nu.csv"), "w") do io
    println(io, "t,KE,Nu  # 公式例と同じ方法（Oceananigans の演算子）で計算")
    for n in 1:length(t)
        @printf(io, "%.3f,%.6g,%.6g\n", t[n], kinetic_energy[n], Nu[n])
    end
end
@info "figures saved" figdir

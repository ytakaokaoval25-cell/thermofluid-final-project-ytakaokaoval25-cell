# 流れ関数 ψ の図（スライド用）
#
# 自作（公式例にはない図）．入力は公式例と同じ計算の出力 output/<名前>.jld2 に保存した u（run_horizontal_convection.jl が保存）．
#
# 流れ関数の定義：u = ∂ψ/∂z，w = −∂ψ/∂x（2 次元・非圧縮なので ∇·u = 0 を自動的に満たす）．
# 計算方法：底 z = −H で ψ = 0（底は不透過）とし，各格子列で u を下から上へ足し上げる．
#   ψ(x_i, z_{k+1}) = ψ(x_i, z_k) + u(x_i, z_k) Δz        （x_i は u の点 = セルの x 面，z_k はセルの z 面）
# こうすると ψ は格子の角（x 面 × z 面）で定義され，離散発散が 0 なら −Δψ/Δx は w と一致する．
# 上面・側面でも ψ ≈ 0（不透過）になる．数値は check_horizontal_convection.jl の Psi と同じ規格化 ψ/(Lx U★)，U★ = √(Lx b★)．
#
# 読み方：ψ の等値線が流線．線の間隔が狭いほど流れが速い．
#   符号は回る向き（ψ < 0 は時計回り，ψ > 0 は反時計回り．x 右向き，z 上向きの図で）．
#   x = 0 で冷やして沈むので，左半分（x < 0）は時計回りで ψ < 0，右半分は反時計回りで ψ > 0 になるはず．
#
# 使い方（リポジトリのルートで）：
#   julia +1.10 --project=. scripts/plot_streamfunction.jl
#   HC_NAME=ra1e8 HC_T0=50 HC_T1=100 HC_FIGURES=02_ra_sweep/results/ra1e8 julia +1.10 --project=. scripts/plot_streamfunction.jl
# 出力：<HC_FIGURES>/<名前>_streamfunction.png（左：平均窓 [HC_T0, HC_T1] の時間平均，右：最終時刻の浮力と流線）
#       <HC_FIGURES>/<名前>_streamfunction.toml（ψ の最大・最小）

using CairoMakie
using Oceananigans
using Printf
using Statistics

name = get(ENV, "HC_NAME", "horizontal_convection")
outdir = get(ENV, "HC_OUTPUT", "output")
figdir = get(ENV, "HC_FIGURES", joinpath("00_reproduce", "results"))
T0 = parse(Float64, get(ENV, "HC_T0", "20"))
T1 = parse(Float64, get(ENV, "HC_T1", "40"))
mkpath(figdir)

H, Lx, b★ = 1.0, 2.0, 1.0
U★ = sqrt(Lx * b★)

path = joinpath(outdir, name * ".jld2")
uts = FieldTimeSeries(path, "u")
bts = FieldTimeSeries(path, "b")
t = uts.times

u1 = interior(uts[1])[:, 1, :]
Nxf, Nz = size(u1)            # u は x 面に乗るので Nxf = Nx + 1
Nx = Nxf - 1
Δx, Δz = Lx / Nx, H / Nz       # 等間隔格子（公式例と同じ）
xf = [-Lx/2 + (i - 1) * Δx for i in 1:Nx+1]     # x 面
zf = [-H + (k - 1) * Δz for k in 1:Nz+1]        # z 面
xc = [-Lx/2 + (i - 0.5) * Δx for i in 1:Nx]     # セル中心
zc = [-H + (k - 0.5) * Δz for k in 1:Nz]

"""u（x 面 × セル中心 z，サイズ (Nx+1, Nz)）から ψ（x 面 × z 面，サイズ (Nx+1, Nz+1)）を求める．"""
function streamfunction(u)
    ψ = zeros(size(u, 1), size(u, 2) + 1)
    for i in axes(u, 1), k in axes(u, 2)
        ψ[i, k+1] = ψ[i, k] + u[i, k] * Δz
    end
    return ψ ./ (Lx * U★)
end

# 平均窓 [T0, T1] の時間平均（保存間隔は一定なので単純平均）
window = findall(τ -> T0 - 1e-8 <= τ <= T1 + 1e-8, t)
isempty(window) && error("平均窓 [$T0, $T1] に保存時刻がない（保存時刻は $(first(t))〜$(last(t))）")
ψ_mean = mean(streamfunction(interior(uts[n])[:, 1, :]) for n in window)

nlast = length(t)
ψ_last = streamfunction(interior(uts[nlast])[:, 1, :])
b_last = interior(bts[nlast])[:, 1, :]

# 上面での ψ（不透過なら 0 に近いはず）：計算の確認用
top_residual = maximum(abs, ψ_last[:, end]) / max(maximum(abs, ψ_last), eps())

# ---- 図（16:9 のスライドに 1 枚で載る横長） ----

ψlim = maximum(abs, ψ_mean)
levels = range(-ψlim, ψlim, length = 17)      # 等値線（流線）は 16 本の間隔

set_theme!(fontsize = 22)
fig = Figure(size = (1600, 620))

axis_kwargs = (xlabel = L"x / H", ylabel = L"z / H",
               limits = ((-Lx/2, Lx/2), (-H, 0)), aspect = DataAspect())

ax1 = Axis(fig[1, 1]; title = @sprintf("time-mean streamfunction, t = %g–%g", t[first(window)], t[last(window)]), axis_kwargs...)
cf = contourf!(ax1, xf, zf, ψ_mean; levels, colormap = :balance)
contour!(ax1, xf, zf, ψ_mean; levels, color = :black, linewidth = 1)
Colorbar(fig[2, 1], cf; vertical = false, label = L"\psi / (L_x \sqrt{L_x b_*})", flipaxis = false)

ax2 = Axis(fig[1, 2]; title = @sprintf("buoyancy and streamlines, t = %g", t[nlast]), axis_kwargs...)
hm = heatmap!(ax2, xc, zc, b_last; colorrange = (-0.6, 0.6), colormap = :thermal)
ψl = maximum(abs, ψ_last)
contour!(ax2, xf, zf, ψ_last; levels = range(-ψl, ψl, length = 17), color = :white, linewidth = 1.2)
Colorbar(fig[2, 2], hm; vertical = false, label = L"b / b_*", flipaxis = false)

Label(fig[0, :], @sprintf("Horizontal convection, Ra = %s  (cooling at x = 0, heating at x = ±1)",
                          get(ENV, "HC_RA", "1e8")), fontsize = 26, tellwidth = false)

figpath = joinpath(figdir, name * "_streamfunction.png")
save(figpath, fig; px_per_unit = 2)

open(joinpath(figdir, name * "_streamfunction.toml"), "w") do io
    println(io, "# 流れ関数 ψ/(Lx U★)（scripts/plot_streamfunction.jl が生成）．ψ<0：時計回り，ψ>0：反時計回り")
    println(io, "window = [", t[first(window)], ", ", t[last(window)], "]")
    println(io, "mean_max = ", maximum(ψ_mean))
    println(io, "mean_min = ", minimum(ψ_mean))
    println(io, "final_time = ", t[nlast])
    println(io, "final_max = ", maximum(ψ_last))
    println(io, "final_min = ", minimum(ψ_last))
    println(io, "final_top_residual_relative = ", top_residual, "  # 上面での |ψ| / max|ψ|（0 に近いはず）")
end

@info "saved" figpath
println(read(joinpath(figdir, name * "_streamfunction.toml"), String))

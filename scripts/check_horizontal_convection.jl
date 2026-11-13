# 手順1の検証：離散発散，浮力分散の収支，純拡散解との比較，Nu・K の時間平均
#
# 使い方（run_horizontal_convection.jl の後に）：
#   julia --project=. scripts/check_horizontal_convection.jl            # output/horizontal_convection.jld2
#   HC_NAME=smoke julia --project=. scripts/check_horizontal_convection.jl
#
# 微分は Oceananigans の演算子を使わず，スタッガード格子上の 2 次差分を
# 配列から独立に計算する（公式例の Nu 計算との照合にもなる）．
#
# 出力：output/<名前>_checks.toml，output/<名前>_timeseries.csv
#
# 変更履歴：平均窓 [20, 40] に保存時刻が 2 点未満のとき（smoke 計算など）は
#           収支 [budget] と時間平均 [means] を省略するようにした（以前は BoundsError で停止）．

using Oceananigans
using Printf

name = get(ENV, "HC_NAME", "horizontal_convection")
outdir = get(ENV, "HC_OUTPUT", "output")
path = joinpath(outdir, name * ".jld2")
Ra = parse(Float64, get(ENV, "HC_RA", "1e8"))

H, Lx, b★, Pr = 1.0, 2.0, 1.0, 1.0
ν = sqrt(Pr * b★ * Lx^3 / Ra)
κ = ν / Pr
U★ = sqrt(Lx * b★)
χ_diff = κ * b★^2 * π * tanh(2π * H / Lx) / (Lx * H)   # 公式例と同じ解析値

bts = FieldTimeSeries(path, "b")
uts = FieldTimeSeries(path, "u")
wts = FieldTimeSeries(path, "w")
t = bts.times
Nt = length(t)

b1 = interior(bts[1])
Nx, Nz = size(b1, 1), size(b1, 3)
Δx, Δz = Lx / Nx, H / Nz
xc = [-Lx/2 + (i - 0.5) * Δx for i in 1:Nx]
zc = [-H + (k - 0.5) * Δz for k in 1:Nz]
bs = [-b★ * cos(2π * x / Lx) for x in xc]          # 上面の浮力（境界条件）

# ---- 離散演算（2 次有限体積，等間隔格子） ----

# 浮力分散散逸 D = ∫ κ|∇b|² dA（内部の面 + 上面の境界面．側面・底面は zero-flux で 0）
function dissipation(b)
    D = 0.0
    for k in 1:Nz, i in 1:Nx-1
        D += ((b[i+1, k] - b[i, k]) / Δx)^2 * Δx * Δz
    end
    for k in 1:Nz-1, i in 1:Nx
        D += ((b[i, k+1] - b[i, k]) / Δz)^2 * Δx * Δz
    end
    for i in 1:Nx   # 上面：値境界条件．面の勾配は (bs - b_N)/(Δz/2)，代表面積 Δx·Δz/2
        D += ((bs[i] - b[i, Nz]) / (Δz / 2))^2 * Δx * (Δz / 2)
    end
    return κ * D
end

# 注：上の D には上面面の寄与を含めたので，対応する B は上面の値 bs で評価する
boundary_supply_face(b) = sum(κ * bs[i] * (bs[i] - b[i, Nz]) / (Δz / 2) * Δx for i in 1:Nx)

# 散逸の鉛直重心による境界層厚さ δχ/H = ∫(-z)χ dA / (H ∫χ dA)（各面の χ をその面の z に置く）
function dissipation_depth(b)
    num, den = 0.0, 0.0
    for k in 1:Nz, i in 1:Nx-1
        c = ((b[i+1, k] - b[i, k]) / Δx)^2 * Δx * Δz
        num += c * (-zc[k]); den += c
    end
    for k in 1:Nz-1, i in 1:Nx
        c = ((b[i, k+1] - b[i, k]) / Δz)^2 * Δx * Δz
        num += c * (-(-H + k * Δz)); den += c
    end
    for i in 1:Nx
        c = ((bs[i] - b[i, Nz]) / (Δz / 2))^2 * Δx * (Δz / 2)
        num += c * (Δz / 4); den += c
    end
    return num / (H * den)
end

variance(b) = 0.5 * sum(abs2, b) * Δx * Δz

function max_divergence(u, w)
    m = 0.0
    for k in 1:Nz, i in 1:Nx
        d = (u[i+1, k] - u[i, k]) / Δx + (w[i, k+1] - w[i, k]) / Δz
        m = max(m, abs(d))
    end
    return m
end

kinetic(u, w) = sum(0.5 * (0.5 * (u[i, k] + u[i+1, k]))^2 + 0.5 * (0.5 * (w[i, k] + w[i, k+1]))^2
                    for k in 1:Nz, i in 1:Nx) * Δx * Δz / (Lx * H)

# 流れ関数 ψ（u = ∂zψ，ψ(z=-H) = 0 から上向きに積分．格子の角で定義）
function psi_max(u)
    m = 0.0
    for i in 1:Nx+1
        ψ = 0.0
        for k in 1:Nz
            ψ += u[i, k] * Δz
            m = max(m, abs(ψ))
        end
    end
    return m
end

# ---- 時系列 ----

Eb = zeros(Nt); D = zeros(Nt); B = zeros(Nt); Nu = zeros(Nt); K = zeros(Nt)
Ψ = zeros(Nt); divmax = zeros(Nt); δ = zeros(Nt)
for n in 1:Nt
    b = interior(bts[n])[:, 1, :]
    u = interior(uts[n])[:, 1, :]
    w = interior(wts[n])[:, 1, :]
    Eb[n] = variance(b)
    D[n] = dissipation(b)
    B[n] = boundary_supply_face(b)
    Nu[n] = D[n] / (Lx * H) / χ_diff
    K[n] = kinetic(u, w) / (Lx * b★)
    Ψ[n] = psi_max(u) / (Lx * U★)
    divmax[n] = max_divergence(u, w)
    δ[n] = dissipation_depth(b)
end

# ---- 収支残差：E_b(t1) - E_b(t0) - ∫(B - D)dt （台形則．保存間隔 0.5 による誤差を含む） ----

function budget_residual(n0, n1)
    integral = sum(0.5 * ((B[n] - D[n]) + (B[n+1] - D[n+1])) * (t[n+1] - t[n]) for n in n0:n1-1)
    ΔE = Eb[n1] - Eb[n0]
    R = ΔE - integral
    scale = sum(0.5 * (D[n] + D[n+1]) * (t[n+1] - t[n]) for n in n0:n1-1)
    return R, R / scale
end

# ---- 純拡散解との比較：解析解 b_diff を同じ離散勾配で評価したときの Nu（→ 1 になるべき） ----

function nu_of_diffusive_solution(Nxd, Nzd)
    Δxd, Δzd = Lx / Nxd, H / Nzd
    x = [-Lx/2 + (i - 0.5) * Δxd for i in 1:Nxd]
    z = [-H + (k - 0.5) * Δzd for k in 1:Nzd]
    s = [-b★ * cos(2π * xi / Lx) for xi in x]
    bd = [s[i] * cosh(2π * (H + z[k]) / Lx) / cosh(2π * H / Lx) for i in 1:Nxd, k in 1:Nzd]
    Dd = 0.0
    for k in 1:Nzd, i in 1:Nxd-1
        Dd += ((bd[i+1, k] - bd[i, k]) / Δxd)^2 * Δxd * Δzd
    end
    for k in 1:Nzd-1, i in 1:Nxd
        Dd += ((bd[i, k+1] - bd[i, k]) / Δzd)^2 * Δxd * Δzd
    end
    for i in 1:Nxd
        Dd += ((s[i] - bd[i, Nzd]) / (Δzd / 2))^2 * Δxd * (Δzd / 2)
    end
    return κ * Dd / (Lx * H) / χ_diff
end
diffusive_table = [(n, n ÷ 2, nu_of_diffusive_solution(n, n ÷ 2)) for n in (32, 64, 128, 256, 512)]

# ---- 平均窓（全条件で共通にする） ----

window(t0, t1) = findall(τ -> t0 <= τ <= t1, t)
mean_over(x, idx) = sum(x[idx]) / length(idx)
w_all, w_a, w_b = window(20, 40), window(20, 30), window(30, 40)

# 窓 [20, 40] に保存時刻が 2 点未満（smoke など短い計算）のときは収支・平均を省略する
has_window = length(w_all) >= 2
if has_window
    n0, n1 = first(w_all), last(w_all)
    R, Rrel = budget_residual(n0, n1)
else
    @warn "平均窓 [20, 40] に十分なデータがないため，収支と時間平均を省略します" final_time = t[end]
end

open(joinpath(outdir, name * "_checks.toml"), "w") do io
    println(io, "# 手順1の検証結果（scripts/check_horizontal_convection.jl が生成）")
    println(io, "file = \"", path, "\"")
    println(io, "Nx = ", Nx, "\nNz = ", Nz, "\nRa = ", Ra, "\nPr = ", Pr)
    println(io, "final_time = ", t[end])
    println(io, "\n[divergence]")
    println(io, "max_abs = ", maximum(divmax))
    println(io, "max_abs_normalized = ", maximum(divmax) / (U★ / Lx), "  # U★/Lx で規格化")
    if has_window
        println(io, "\n[budget]  # 窓 t = ", t[n0], "〜", t[n1])
        println(io, "delta_Eb = ", Eb[n1] - Eb[n0])
        println(io, "integral_B_minus_D = ", (Eb[n1] - Eb[n0]) - R)
        println(io, "residual = ", R)
        println(io, "residual_relative_to_integral_D = ", Rrel)
    else
        println(io, "\n# [budget] と [means] は省略：平均窓 [20, 40] にデータがない（final_time < 20）")
    end
    println(io, "\n[diffusive_solution]  # 解析解を同じ離散勾配で評価した Nu（1 に近づくべき）")
    for (a, c, v) in diffusive_table
        @printf(io, "Nu_%dx%d = %.6f\n", a, c, v)
    end
    if has_window
        println(io, "\n[means]  # 窓 [20, 40], 前半 [20, 30], 後半 [30, 40]")
        for (lab, x) in (("Nu", Nu), ("K", K), ("Psi", Ψ), ("delta_chi", δ))
            @printf(io, "%s = [%.5g, %.5g, %.5g]\n", lab, mean_over(x, w_all), mean_over(x, w_a), mean_over(x, w_b))
        end
    end
end

open(joinpath(outdir, name * "_timeseries.csv"), "w") do io
    println(io, "t,Nu,K,Psi,delta_chi,Eb,B,D,max_div")
    for n in 1:Nt
        @printf(io, "%.3f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.3e\n", t[n], Nu[n], K[n], Ψ[n], δ[n], Eb[n], B[n], D[n], divmax[n])
    end
end

println(read(joinpath(outdir, name * "_checks.toml"), String))

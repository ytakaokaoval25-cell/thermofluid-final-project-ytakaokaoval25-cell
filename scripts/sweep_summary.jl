# 手順2（issue #5）：Ra スイープのまとめ
#
# 入力：results/<ケース名>/<ケース名>_checks.toml（check スクリプト），
#       results/<ケース名>/<ケース名>_official_Nu.csv（plot スクリプト：公式例と同じ演算子・数値の χ_diff）
# 出力：results/ra_sweep/summary.csv，results/ra_sweep/Nu_Ra.png，results/ra_sweep/fit.toml
#
# 使い方：julia --project=. scripts/sweep_summary.jl
#   平均窓は HC_T0, HC_T1（既定 50, 100．check と同じ窓にすること）

using CairoMakie
using Printf
using TOML

T0 = parse(Float64, get(ENV, "HC_T0", "50"))
T1 = parse(Float64, get(ENV, "HC_T1", "100"))
Tm = (T0 + T1) / 2
resdir = get(ENV, "HC_RESULTS", "results")
outdir = joinpath(resdir, "ra_sweep")
mkpath(outdir)

sweep = ["ra1e6", "ra1e7", "ra1e8", "ra1e9"]   # 128×64 の 4 点（傾きの推定に使う）
fine = "ra1e9_256x128"
cases = [sweep; fine]

"""公式例と同じ演算子で求めた Nu の時系列（t, Nu）を読む．# で始まる行と見出しは飛ばす．"""
function read_official(c)
    t, Nu = Float64[], Float64[]
    for line in eachline(joinpath(resdir, c, c * "_official_Nu.csv"))
        (startswith(line, "#") || startswith(line, "t,")) && continue
        f = split(line, ",")
        push!(t, parse(Float64, f[1])); push!(Nu, parse(Float64, f[3]))
    end
    return t, Nu
end

function window_mean(t, x, a, b)
    idx = findall(τ -> a <= τ <= b, t)
    isempty(idx) && error("窓 [$a, $b] にデータがありません")
    return sum(x[idx]) / length(idx)
end

function summarize(c)
    chk = TOML.parsefile(joinpath(resdir, c, c * "_checks.toml"))
    run = TOML.parsefile(joinpath(resdir, c, c * "_run.toml"))
    w = chk["window"]
    (w["t0"] == T0 && w["t1"] == T1) || error("$c の checks.toml の窓 [$(w["t0"]), $(w["t1"])] が [$T0, $T1] と違います")
    t, Nu = read_official(c)
    m = chk["means"]
    return (; case = c, Ra = chk["Ra"], Nx = chk["Nx"], Nz = chk["Nz"],
            Nu = window_mean(t, Nu, T0, T1), Nu_a = window_mean(t, Nu, T0, Tm), Nu_b = window_mean(t, Nu, Tm, T1),
            Nu_check = m["Nu"][1], Psi = m["Psi"][1], Psi_a = m["Psi"][2], Psi_b = m["Psi"][3],
            K = m["K"][1], delta_chi = m["delta_chi"][1],
            budget_residual = chk["budget"]["residual_relative_to_integral_D"],
            max_div = chk["divergence"]["max_abs_normalized"],
            iterations = run["iterations"], wall = run["wall_time_seconds"])
end

rows = [summarize(c) for c in cases]
halfdiff(a, b) = (b - a) / ((a + b) / 2)   # 後半 − 前半（平均で規格化）

open(joinpath(outdir, "summary.csv"), "w") do io
    println(io, "# 窓 [$T0, $T1]（前半 [$T0, $Tm]，後半 [$Tm, $T1]）．Nu は公式例と同じ演算子・数値の χ_diff，Nu_check は照合用（面ベース・解析値の χ_diff）")
    println(io, "case,Ra,Nx,Nz,Nu,Nu_first,Nu_second,Nu_halfdiff,Nu_check,Psi,Psi_halfdiff,K,delta_chi,budget_residual_rel,max_div_norm,iterations,run_wall_s")
    for r in rows
        @printf(io, "%s,%.0e,%d,%d,%.5g,%.5g,%.5g,%.3g,%.5g,%.5g,%.3g,%.5g,%.5g,%.3g,%.3g,%d,%.1f\n",
                r.case, r.Ra, r.Nx, r.Nz, r.Nu, r.Nu_a, r.Nu_b, halfdiff(r.Nu_a, r.Nu_b), r.Nu_check,
                r.Psi, halfdiff(r.Psi_a, r.Psi_b), r.K, r.delta_chi, r.budget_residual, r.max_div, r.iterations, r.wall)
    end
    # Ra = 1e9 の解像度差：(128×64 − 256×128) / 256×128
    c, f = rows[end-1], rows[end]
    rel(a, b) = (a - b) / b
    @printf(io, "ra1e9_reldiff_128_vs_256,%.0e,,,%.3g,%.3g,%.3g,,%.3g,%.3g,,%.3g,%.3g,,,,\n",
            c.Ra, rel(c.Nu, f.Nu), rel(c.Nu_a, f.Nu_a), rel(c.Nu_b, f.Nu_b), rel(c.Nu_check, f.Nu_check),
            rel(c.Psi, f.Psi), rel(c.K, f.K), rel(c.delta_chi, f.delta_chi))
end

# ---- 両対数の最小二乗（128×64 の 4 点） ----

function loglog_slope(x, y)
    X, Y = log10.(x), log10.(y)
    Xm, Ym = sum(X) / length(X), sum(Y) / length(Y)
    s = sum((X .- Xm) .* (Y .- Ym)) / sum((X .- Xm) .^ 2)
    return s, Ym - s * Xm
end

base = rows[1:length(sweep)]
Ra = [r.Ra for r in base]
fits = Dict(lab => loglog_slope(Ra, [getproperty(r, sym) for r in base])
            for (lab, sym) in (("Nu", :Nu), ("Nu_check", :Nu_check), ("Psi", :Psi), ("delta_chi", :delta_chi)))
# 上位 3 点（Ra ≥ 1e7）だけの傾きも参考に出す
fits3 = loglog_slope(Ra[2:end], [r.Nu for r in base[2:end]])

open(joinpath(outdir, "fit.toml"), "w") do io
    println(io, "# 両対数の最小二乗の傾き（128×64，Ra = 1e6〜1e9 の 4 点，窓 [$T0, $T1]）")
    for k in ("Nu", "Nu_check", "Psi", "delta_chi")
        @printf(io, "slope_%s = %.4f\n", k, fits[k][1])
    end
    @printf(io, "slope_Nu_Ra_ge_1e7 = %.4f  # 参考：Ra = 1e7〜1e9 の 3 点\n", fits3[1])
    println(io, "reference_slope = 0.2  # Rossby (1965) Nu ∝ Ra^(1/5)")
end

# ---- 図 ----

fig = Figure(size = (650, 520))
ax = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = L"Ra", ylabel = L"\langle Nu \rangle",
          title = @sprintf("time mean over t = [%g, %g]", T0, T1))
s, a = fits["Nu"]
Rline = 10 .^ range(5.7, 9.3; length = 50)
lines!(ax, Rline, 10^a .* Rline .^ s; color = :black, label = @sprintf("fit (128×64): slope %.3f", s))
Nu8 = base[3].Nu
lines!(ax, Rline, Nu8 .* (Rline ./ 1e8) .^ 0.2; color = :gray, linestyle = :dash, label = "slope 1/5 (through Ra = 1e8)")
scatter!(ax, Ra, [r.Nu for r in base]; markersize = 14, label = "Nu (official operators, 128×64)")
scatter!(ax, Ra, [r.Nu_check for r in base]; marker = :rect, markersize = 10, color = :transparent,
         strokewidth = 1.5, strokecolor = :orange, label = "Nu_check (face-based, 128×64)")
scatter!(ax, [rows[end].Ra], [rows[end].Nu]; marker = :diamond, markersize = 16, color = :red, label = "Nu (256×128)")
axislegend(ax; position = :lt)
save(joinpath(outdir, "Nu_Ra.png"), fig)

println(read(joinpath(outdir, "summary.csv"), String))
println(read(joinpath(outdir, "fit.toml"), String))

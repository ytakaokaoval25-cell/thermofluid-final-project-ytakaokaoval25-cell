# 00_reproduce：公式例の再現（手順1，issue #3）

Oceananigans.jl v0.110.15 の公式ギャラリー例 `examples/horizontal_convection.jl`（無改変コピーは `upstream/`）を，
Mac の CPU でそのまま再現する．

| 条件 | 値 |
|---|---|
| Ra，Pr | 1e8，1 |
| 格子 | 128 × 64（x：Bounded，y：Flat） |
| 数値 | `NonhydrostaticModel`，5 次 WENO，3 段 RK，CFL 0.7（`conjure_time_step_wizard!`），max Δt = 0.1 |
| 境界 | 上面 b = −b★ cos(2πx/Lx)，他は zero-flux．速度は全境界 free-slip |
| 終了 | t = 40（0.5 ごとに JLD2 保存） |

## 実行方法（リポジトリのルートで）

再現には Oceananigans.jl の公式例そのもの（`upstream/horizontal_convection.jl`）を使う．
(B) の `scripts/run_horizontal_convection.jl` は公式例と物理・数値設定が同じで，変更点はファイル冒頭に列挙している（自作の書き直しではない）．

```sh
julia +1.10 --project=. -e 'using Pkg; Pkg.instantiate()'   # 初回のみ

# 小さな格子での動作確認（64×32, t = 2，1〜2 分）
HC_NX=64 HC_NZ=32 HC_STOP_TIME=2 HC_NAME=smoke julia +1.10 --project=. scripts/run_horizontal_convection.jl

# (A) 公式例をそのまま実行（upstream/horizontal_convection.jl は無改変．出力は output/upstream/）
sh scripts/run_upstream_example.sh
sh scripts/run_upstream_example.sh   # 2 回目（コンパイル済みでの時間）

# (B) 検証用：公式例に最小限の変更（出力先・計時・u, w の保存）を加えた版で run → check → plot をまとめて実行し，実行時間を results/timing_raw.txt に追記
sh scripts/run_reproduce.sh
sh scripts/run_reproduce.sh   # 2 回目（コンパイル済みでの時間を測る）
```

1 つずつ実行する場合は `scripts/run_horizontal_convection.jl` → `check_horizontal_convection.jl` →
`HC_FIGURES=00_reproduce/results ... plot_horizontal_convection.jl`．
大きな出力（JLD2，mp4，ログ）は `output/`（commit しない）にできる．

## 結果（2026-11-13，Apple M2 / 8 GB，Julia 1.10.13，Oceananigans 0.110.15，1 スレッド）

- 図：[`results/horizontal_convection_fields_final.png`](results/horizontal_convection_fields_final.png)（t = 40 の b，ζ，χ），
  [`results/horizontal_convection_Nu_KE.png`](results/horizontal_convection_Nu_KE.png)（Nu(t)，K(t)）
- 離散発散 max|∇·u|/(U★/Lx) = 5.4e-6，浮力分散の収支残差 2%（窓 [20, 40]）
- 窓 [20, 40] の平均 Nu ≈ 9.2（面ベースの照合用定義．公式例の定義と K の規格化は手順2で確定：[../02_ra_sweep/](../02_ra_sweep/)）
- 実行時間：[`results/timing.txt`](results/timing.txt)．本計算 run は実時間 1 分 51 秒（run! 74.6 秒，そのうち約 45 秒がコンパイル），
  初回の `Pkg` 準備に約 16 分．2 回目（コンパイル済み）の時間は `run_reproduce.sh` で追記する．

## 流れ関数の図（スライド用，issue #11）

`scripts/plot_streamfunction.jl`（自作．公式例にはない図）で，循環が目で分かる流れ関数 ψ の図を作る．

```sh
julia +1.10 --project=. scripts/plot_streamfunction.jl     # 本計算（output/horizontal_convection.jld2）の後に実行
```

出力：`results/horizontal_convection_streamfunction.png`（左：t = 20〜40 の時間平均 ψ，右：t = 40 の浮力と流線），
`results/horizontal_convection_streamfunction.toml`（ψ の最大・最小）．

作り方：
1. 公式例と同じ計算で保存した水平速度 u（0.5 ごとの JLD2）を読む．
2. 定義 u = ∂ψ/∂z，w = −∂ψ/∂x のうち 1 つ目を使い，底 z = −H で ψ = 0（不透過）として各列で u Δz を下から足し上げる．
   ψ は格子の角に乗る．離散発散が 0（ここでは 5e-6）なら，−∂ψ/∂x は w と自動的に一致し，上面でも ψ ≈ 0 になる（toml に上面の残差を書く）．
3. 規格化は ψ/(Lx U★)，U★ = √(Lx b★)（check スクリプトの Psi と同じ）．
4. 時間平均は平均窓 [20, 40] の保存時刻の単純平均．等値線（流線）は ±max|ψ| を 16 等分した値で引く．

読み方：等値線が流線で，間隔が狭いほど速い．x = 0 で冷やして沈むので，左半分は時計回り（ψ < 0），右半分は反時計回り（ψ > 0）．
（この計算方法は，既知の ψ から作った非発散の速度場で，丸め誤差の範囲で ψ が戻ることを確認した．）

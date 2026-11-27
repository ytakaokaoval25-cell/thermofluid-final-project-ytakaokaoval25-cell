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

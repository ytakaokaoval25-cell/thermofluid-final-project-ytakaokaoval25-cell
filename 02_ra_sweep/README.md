# 02_ra_sweep：Ra スイープ（手順2，issue #5）

Ra = 1e6, 1e7, 1e8, 1e9（128×64）と Ra = 1e9（256×128）を t = 100 まで計算し，窓 [50, 100] で平均して
Nu ∝ Ra^(1/5)（Rossby 1965）を確かめる．

## 実行方法（リポジトリのルートで，1 本ずつ順番に）

```sh
sh scripts/run_ra_sweep.sh                 # 全 5 ケース．一部なら sh scripts/run_ra_sweep.sh ra1e6 ra1e7
HC_T0=50 HC_T1=100 julia +1.10 --project=. scripts/sweep_summary.jl   # まとめ表・図
```

各ケースの図と小さな結果は `results/<ケース名>/`，まとめは `results/ra_sweep/`（`summary.csv`，`fit.toml`，`Nu_Ra.png`，`timing.txt`）．

## 定義

- Nu：公式例と同じ演算子 `κ(∂x b² + ∂z b²)` を `∫dV` し，分母 χ_diff は純拡散解に同じ演算子を適用した数値（`scripts/nusselt_reference.jl`）．
  check スクリプトの面ベースの Nu は照合用．
- K：運動エネルギー /(Lx b★)．
- Ra = 1e6 だけ拡散 CFL 上限 0.2（`HC_DIFFUSIVE_CFL=0.2`）．既定（公式例どおり Inf）では t ≈ 15 で NaN．

## 結果（2026-11-13，Apple M2）

| Ra | 格子 | Nu | 
|---|---|---|
| 1e6 | 128×64 | 3.42 |
| 1e7 | 128×64 | 5.08 |
| 1e8 | 128×64 | 7.09 |
| 1e9 | 128×64 | 12.0 |
| 1e9 | 256×128 | 17.2 |

傾き 0.178（照合用定義では 0.214），理論 0.2．Ra = 1e9 は 128×64 では解像度不足．

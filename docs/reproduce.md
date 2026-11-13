# 再現手順

## 手順1：公式例の再現（Ra = 1e8, 128×64）

リポジトリのルートで実行する．

```sh
# 1回だけ：依存関係の復元
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# 小さな CPU 例（64×32, t = 2）
HC_NX=64 HC_NZ=32 HC_STOP_TIME=2 HC_NAME=smoke julia --project=. scripts/run_horizontal_convection.jl
HC_NAME=smoke julia --project=. scripts/check_horizontal_convection.jl

# 本計算（公式例と同じ 128×64, t = 40, Ra = 1e8）
julia --project=. scripts/run_horizontal_convection.jl
julia --project=. scripts/check_horizontal_convection.jl
julia --project=. scripts/plot_horizontal_convection.jl
```

`output/` に計算結果（JLD2，動画，実行記録），`results/` に図と小さな結果ファイルができる．
公式コードの出典と変更点は `upstream/README.md` と各スクリプトの冒頭に記録している．

## 手順2：Ra スイープ（issue #5）

Ra = 1e6, 1e7, 1e8, 1e9（128×64）と Ra = 1e9（256×128）を t = 100 まで計算し，窓 [50, 100] で平均する．
1 本ずつ順番に実行する（並列にしない）．スレッド数は 1（`JULIA_NUM_THREADS=1`）．

```sh
# 5 ケースを順番に run → check（HC_T0=50 HC_T1=100）→ plot．結果は results/<ケース名>/
sh scripts/run_ra_sweep.sh                 # 一部だけなら：sh scripts/run_ra_sweep.sh ra1e6 ra1e7

# 1 ケースだけ手で回す場合（例：Ra = 1e7）
HC_NAME=ra1e7 HC_RA=1e7 HC_STOP_TIME=100 julia --project=. scripts/run_horizontal_convection.jl
HC_NAME=ra1e7 HC_RA=1e7 HC_T0=50 HC_T1=100 julia --project=. scripts/check_horizontal_convection.jl
HC_NAME=ra1e7 HC_RA=1e7 HC_FIGURES=results/ra1e7 julia --project=. scripts/plot_horizontal_convection.jl

# 流れなしで Nu = 1 になることの確認（純拡散，128×64）
julia --project=. scripts/check_diffusive_limit.jl

# まとめ：results/ra_sweep/summary.csv，fit.toml，Nu_Ra.png
HC_T0=50 HC_T1=100 julia --project=. scripts/sweep_summary.jl
```

定義（issue #5）：
- Nu は公式例と同じ演算子 `@at (Center, Center, Center) κ * (∂x(b)^2 + ∂z(b)^2)` と `∫dV` で求め，
  分母 χ_diff は純拡散解を同じ境界条件の Field に入れて同じ演算子で計算した数値（`scripts/nusselt_reference.jl`）．
  check スクリプトの面ベースの Nu（分母は解析値）は照合用．
- 運動エネルギーは K/(Lx b★)．

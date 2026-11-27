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

# 熱流体力学演習 2026 最終プロジェクト

Oceananigans.jl による水平対流と海底地形の計算．計画は issue #1 で作成中．

## フォルダ構成（issue #8）

手順（課題）ごとにフォルダを分ける．各フォルダの README.md に条件・実行方法・結果を書き，図と小さな結果ファイルはその下の `results/` に置く．

| フォルダ | 内容 |
|---|---|
| [`00_reproduce/`](00_reproduce/) | 手順1：公式例（Ra = 1e8，128×64，t = 40）の再現 |
| [`01_validation/`](01_validation/) | 検証：離散発散，収支，純拡散で Nu = 1，格子依存性 |
| [`02_ra_sweep/`](02_ra_sweep/) | 手順2：Ra スイープ（Nu ∝ Ra^(1/5)） |
| [`03_ridge/`](03_ridge/) | 手順3：海嶺の高さと深層循環 |
| [`04_extensions/`](04_extensions/) | 発展：加熱分布，回転 |
| [`slides/`](slides/) | 発表に使う図の一覧 |
| `scripts/` | 全手順で共通のスクリプト（run / check / plot など） |
| `upstream/` | Oceananigans 公式コードの無改変コピー |
| `output/` | 大きな計算出力（JLD2，mp4，ログ）．commit しない |

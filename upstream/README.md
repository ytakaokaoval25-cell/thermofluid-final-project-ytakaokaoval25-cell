# 上流（Oceananigans.jl）からの無改変コピー

| ファイル | 出典 | 用途 |
|---|---|---|
| `horizontal_convection.jl` | https://github.com/CliMA/Oceananigans.jl/blob/v0.110.15/examples/horizontal_convection.jl | 手順1（公式例の再現）の元 |
| `internal_tide.jl` | https://github.com/CliMA/Oceananigans.jl/blob/v0.110.15/examples/internal_tide.jl | 手順3（海嶺）の地形の作り方の参考 |

- タグ：v0.110.15（課題ページ記載の commit 5e62a72）
- ライセンス：MIT License（Copyright (c) 2018 Climate Modeling Alliance）
- 取得日：2026-11-13（raw.githubusercontent.com から取得）
- SHA-256：
  - horizontal_convection.jl `757fdd704c57ac16d0183d28e765adc3a2cbf31ff751f899d2301dbb5ceac063`
  - internal_tide.jl `1402dcc6ea13bac0d5b51543c2f95401466be520021171bb4c995e032a8c13c1`

メモ：v0.110.15 の `internal_tide.jl` は山を `PartialCellBottom(bottom)` で表している（セルを部分的に削る方式）．
`GridFittedBottom` はセル単位で塗りつぶす方式で，どちらを使うかは手順3で比べて決める．

これらのファイルは編集しない．変更は `scripts/` 側で行い，変更点を冒頭に書く．

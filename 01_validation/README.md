# 01_validation：検証

計算と Nu の定義が正しいことを確かめる．再現（00）と Ra スイープ（02）の両方に効く確認をここにまとめる．

| 確認 | 方法 | 結果 |
|---|---|---|
| 離散発散 | check スクリプトで max\|∇·u\| | 各ケース 5e-6 程度（U★/Lx 規格化）．`../00_reproduce/results/*_checks.toml`，`../02_ra_sweep/results/ra_sweep/summary.csv` |
| 浮力分散の収支 | 窓内の ΔE_b と ∫(B−D)dt の差 | 00 で 2%．Ra=1e9（128×64）で 7% |
| 純拡散で Nu = 1 | `scripts/check_diffusive_limit.jl`（流れを切って拡散だけ解く） | [`results/diffusive_limit_128x64.toml`](results/diffusive_limit_128x64.toml) |
| 格子依存性 | Ra=1e9 を 128×64 と 256×128 で比較 | Nu が 30% 違う → 128×64 は Ra=1e9 で解像度不足（02 を参照） |

```sh
julia +1.10 --project=. scripts/check_diffusive_limit.jl   # 128×64 で約 9 分．結果は 01_validation/results/
```

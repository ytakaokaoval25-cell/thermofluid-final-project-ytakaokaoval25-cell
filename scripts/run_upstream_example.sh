#!/bin/sh
# 公式例 upstream/horizontal_convection.jl（Oceananigans v0.110.15，無改変）をそのまま実行する．
# 使い方（リポジトリのルートで）：sh scripts/run_upstream_example.sh
#   公式例はカレントディレクトリに horizontal_convection.jld2 / .mp4 を書くので，output/upstream/ に移動して実行する．
#   実行時間（秒）は 00_reproduce/results/timing_raw.txt に追記する．2 回続けて実行すると初回（コンパイル込み）と 2 回目を比べられる．
set -eu
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-1}"
JULIA="${JULIA:-julia +1.10}"
ROOT=$(pwd)
mkdir -p output/upstream 00_reproduce/results
cd output/upstream
start=$(date +%s)
$JULIA --project="$ROOT" "$ROOT/upstream/horizontal_convection.jl" 2>&1 | tee upstream_run.log
end=$(date +%s)
echo "$(date '+%Y-%m-%d %H:%M') upstream_example(run+plot) $((end - start))" >> "$ROOT/00_reproduce/results/timing_raw.txt"
echo "完了：$((end - start)) 秒．出力は output/upstream/（horizontal_convection.jld2，horizontal_convection.mp4）"

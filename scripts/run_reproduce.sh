#!/bin/sh
# 手順1（00_reproduce）：公式例（Ra = 1e8, 128×64, t = 40）を run → check → plot の順に実行する．
# 使い方（リポジトリのルートで）：sh scripts/run_reproduce.sh
#   同じ処理を 2 回続けて実行すると，1 回目（JIT コンパイル込み）と 2 回目の実行時間を比べられる．
#   実行時間（秒）は 00_reproduce/results/timing_raw.txt に追記する．
#   小さな結果ファイル（*_run.toml，*_checks.toml，*_timeseries.csv）は 00_reproduce/results/ にコピーする．
set -eu
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-1}"
JULIA="${JULIA:-julia +1.10}"
R=00_reproduce/results
mkdir -p "$R" output
n=horizontal_convection
export HC_NAME="$n" HC_FIGURES="$R"

timed() {  # timed <ラベル> <コマンド...>
  label="$1"; shift
  start=$(date +%s)
  "$@"
  end=$(date +%s)
  echo "$(date '+%Y-%m-%d %H:%M') $label $((end - start))" >> "$R/timing_raw.txt"
  echo "  $label: $((end - start)) 秒"
}

timed "run"   $JULIA --project=. scripts/run_horizontal_convection.jl   > "output/${n}_run.log" 2>&1
timed "check" $JULIA --project=. scripts/check_horizontal_convection.jl > "output/${n}_check.log" 2>&1
timed "plot"  $JULIA --project=. scripts/plot_horizontal_convection.jl  > "output/${n}_plot.log" 2>&1
cp "output/${n}_run.toml" "output/${n}_checks.toml" "output/${n}_timeseries.csv" "$R/"
echo "完了：図と結果は $R/，ログは output/${n}_*.log"

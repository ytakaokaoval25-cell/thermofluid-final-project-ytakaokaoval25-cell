#!/bin/sh
# 手順2（issue #5）：Ra スイープを 1 本ずつ順番に実行する（並列にしない）．
# 使い方（リポジトリのルートで）：sh scripts/run_ra_sweep.sh [ケース名 ...]
#   引数なしなら全ケース．各ケースで run → check（窓 [50, 100]）→ plot を行い，
#   小さな結果ファイルを results/<ケース名>/ にコピーする．実行時間は results/ra_sweep/timing_raw.txt に追記する．
#   Ra = 1e6 だけ時間刻みの拡散 CFL 上限 0.2 をかける（既定の Inf では t ≈ 15 で拡散の安定限界を超えて NaN）．
set -eu
export JULIA_NUM_THREADS=1
JULIA="${JULIA:-julia +1.10}"
mkdir -p results/ra_sweep output

case_args() {
  case "$1" in
    ra1e6)          echo "1e6 128 64 0.2" ;;
    ra1e7)          echo "1e7 128 64 Inf" ;;
    ra1e8)          echo "1e8 128 64 Inf" ;;
    ra1e9)          echo "1e9 128 64 Inf" ;;
    ra1e9_256x128)  echo "1e9 256 128 Inf" ;;
    *) echo "unknown case $1" >&2; exit 1 ;;
  esac
}

timed() {  # timed <ラベル> <コマンド...>：実時間（秒）を timing_raw.txt に記録
  label="$1"; shift
  start=$(date +%s)
  "$@"
  end=$(date +%s)
  echo "$label $((end - start))" >> results/ra_sweep/timing_raw.txt
}

CASES="${*:-ra1e6 ra1e7 ra1e8 ra1e9 ra1e9_256x128}"
for c in $CASES; do
  set -- $(case_args "$c"); ra=$1; nx=$2; nz=$3; dcfl=$4
  mkdir -p "results/$c"
  export HC_NAME="$c" HC_RA="$ra" HC_NX="$nx" HC_NZ="$nz" HC_STOP_TIME=100 HC_DIFFUSIVE_CFL="$dcfl" HC_T0=50 HC_T1=100 HC_FIGURES="results/$c"
  echo "=== $c (Ra=$ra, ${nx}x${nz}) ==="
  timed "$c run"   $JULIA --project=. scripts/run_horizontal_convection.jl   > "output/${c}_run.log" 2>&1
  timed "$c check" $JULIA --project=. scripts/check_horizontal_convection.jl > "output/${c}_check.log" 2>&1
  timed "$c plot"  $JULIA --project=. scripts/plot_horizontal_convection.jl  > "output/${c}_plot.log" 2>&1
  cp "output/${c}_run.toml" "output/${c}_checks.toml" "output/${c}_timeseries.csv" "results/$c/"
  echo "=== $c done ==="
done

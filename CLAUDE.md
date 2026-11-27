# 作業ルール（人間・AI 共通）

このリポジトリは熱流体力学演習 2026 の最終プロジェクト「Oceananigans.jl による水平対流：海底地形が深層循環に与える影響」です．
問い・計画・計算条件は README.md が正本です．作業前に README.md を読んでください．

複数の作業者（本人，VS Code のローカル Claude Code，claude.ai の Claude など）が同じリポジトリを触るため，次のルールを守ります．

## Git の流れ

- 必ず **issue → branch → 実装 → ローカルテスト → PR → レビュー → merge** の順に進める（課題ページの決まり）．
- `main` へ直接 commit / push しない．
- PR は作成者以外がレビューしてから merge する．AI は merge しない．
- 作業を始める前に `git pull` で最新の `main` を取り込み，そこからブランチを切る．
- ブランチ名は `issue-<番号>-<内容>`（例：`issue-3-reproduce-official-example`）．
- 作業者ごとに別のブランチを使い，同じファイルを同時に編集しない．着手するときは issue に「作業中」とコメントするか自分を assign する．
- PR には「AI を使った箇所」「提案の採否と理由」「確認方法」を書く．

## 記録

- 決めたこと・進み具合・計算結果は，チャットではなくリポジトリに残す（README.md の該当節，または `docs/progress.md`）．
- 計算を実行したら，格子，Ra，終了時刻，スレッド数，実行時間（初回コンパイルと 2 回目を分ける），Julia と Oceananigans の版を記録する．`scripts/run_horizontal_convection.jl` は `output/<名前>_run.toml` に自動で書き出す．`scripts/run_reproduce.sh`，`scripts/run_ra_sweep.sh` は実行時間を各手順の `results/` に追記する．
- 公式コードを流用したら，出典（URL・タグ・commit）とライセンス，変更点をファイル冒頭とPRに書く．無改変のコピーは `upstream/` に置く．

## 環境

- Julia と Oceananigans の版は `Project.toml` / `Manifest.toml` で固定する（Oceananigans v0.110.15，Julia は 1.10 系）．勝手に `Pkg.update` しない．
- 計算は MacBook Air（Apple Silicon）の CPU で行う．GPU は使わない．
- 大きな出力（`*.jld2`，`*.mp4`，`output/`）は commit しない．図（png）と小さな結果ファイル（toml，csv）は各手順フォルダの `results/` に置いて commit する．

## フォルダ構成

- 手順ごとに `00_reproduce/`，`01_validation/`，`02_ra_sweep/`，`03_ridge/`，`04_extensions/` に分ける．各フォルダに README.md（条件・実行方法・結果）と `results/` を置く．
- スクリプトは全手順共通で `scripts/` に置く（場所を変えない）．環境変数 `HC_NX, HC_NZ, HC_STOP_TIME, HC_RA, HC_NAME, HC_OUTPUT, HC_FIGURES` の意味は変えない．
- 発表用の図は `slides/README.md` に一覧を書き，図そのものは各手順の `results/` から使う．

## 文体

- 返答・コメント・ドキュメントは日本語で書く．

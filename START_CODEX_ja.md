# Codex起動指示 — nc-p3 / 2026-11-15 / v16

GitHub uedakazushi/nc-p3 のmain先端を確認する。
研究方針 research/plan_ja.tex → AGENTS.md → 主論文 nc-p3.tex の順に読み、
必要な証明位置は research/STATE.json、計算・障害は research/notes_ja.tex、
一次文献は research/BASELINE.md から調べる。

最優先課題は ALGEBRAIC-FAMILY-008：
九変数の楕円Clifford–Ore族の一般点で全変形と既存族の接像を比較し、
成分としての最大性と族内同型判定を進める。不足方向は実際の代数的族で積分する。
具体的手順と優先順位は研究方針に集約している。

既に証明された境界包含・一般最小周期2を未解決へ戻さない。
接空間、族の像、既約成分の次元を区別する。
原稿・ノート・方針・索引を整合させ、READMEの検証を実行した範囲だけ記録する。

通常fast-forwardでmainへ直接保存し、遠隔を読み戻す。
PR・force push・GitHub Actions・無断のas4書込みをしない。
PDFとキャッシュはGit追跡せず、GitHub保存成功時はLibraryに重複保存しない。

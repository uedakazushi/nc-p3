# Codex起動指示 — nc-p3 / 2026-11-15 / v10

`uedakazushi/nc-p3` の `main` を確認し、研究方針 → AGENTS → STATE → 主論文 → 研究ノート → BASELINE の順に読む。主論文は `nc-p3.tex`、九変数族の節は `sections/elliptic_ore_family.tex`。

最優先は **ALGEBRAIC-FAMILY-008**。実際の代数的変形を求める。九変数の楕円Clifford–Ore族は全ファイバーAS正則、一般の最小周期2、同型類の像の次元9。七変数族は `t=h=0` に保持される。この成果の再実装を繰り返さず、楕円曲線と次数12の因子による同型の必要十分条件、像の閉包と最大性、別の三生成部分による追加方向を追う。

偶6＋交代3方向の積分から全15方向や14次元成分を宣言しない。`1+t*h != 0` の条件を保ち、一般のt,hで元の関係が二次Gröbner基底だと仮定しない。形式的障害計算は候補選別の補助手段とする。

READMEの検証を実行し、成功したものだけをPASSにする。証明・ノート・方針・STATEを現行内容へ統合してmainへ通常fast-forwardでpushし、遠隔を読み戻す。PR・force push・GitHub Actions・無断のas4書込みをしない。PDFとキャッシュはGitで追跡しない。GitHubへのpushが成功した場合はLibraryに保存しない。

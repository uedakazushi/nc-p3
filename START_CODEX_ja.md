# Codex起動指示 — nc-p3 / 2026-11-15 / v15

`uedakazushi/nc-p3` の `main` を確認し、研究方針 → AGENTS → STATE → 主論文 → 研究ノート → BASELINE の順に読む。主論文は `nc-p3.tex`、九変数族の節は `sections/elliptic_ore_family.tex`。

最優先は **ALGEBRAIC-FAMILY-008**。実際の代数的変形を求める。九変数の楕円Clifford–Ore族は全ファイバーAS正則、一般の最小周期2、同型類の像の次元9。七変数族は `t=h=0` に保持される。

**新しい有効な成果**：三生成Clifford部分Bを固定した一次Ore導分は通常6、符号捩れ9に完全分類済み（内導分・基底変更を除いて9係数）。九変数正規形はその構成内で一般の族を尽くす。楕円曲線の$j$-不変量`J(th)`の明示式は同型の必要条件となる（次数12の有理写像）。原稿`sections/ore_complete.tex`の`thm:ore-derivations`, `cor:ore-slice`, `prop:elliptic-j-obstruction`を参照。同型の十分条件、全モジュライ成分性は未証明。この成果の再実装を繰り返さず、楕円曲線と次数12の因子による同型の必要十分条件、像の閉包と最大性、別の三生成部分による追加方向を追う。

偶6＋交代3方向の積分から全15方向や14次元成分を宣言しない。`1+t*h != 0` の条件を保ち、一般のt,hで元の関係が二次Gröbner基底だと仮定しない。形式的障害計算は候補選別の補助手段とする。

`sections/isogeny_invariant.tex`（内在的2-isogeny、次数4の同型必要不変量）と`sections/nonore_family.tex`（全ファイバーAS正則な五変数族、像の次元2、別の交代接方向）を読む。非Ore族の一般周期2は`sections/nonore_period_two.tex`で証明済み。新たな厳密検算は`computations/verify_nonore_point_geometry.py`。九変数の軌道飽和閉包もClifford点で特異：次元39に対して接空間41以上、境界退化を含め偶7・交代4の個別に代数的積分可能な11方向を確保。ただし二次障害`lambda*c`によって滑らかな11方向同時積分は不可能。証明は`sections/nine_closure.tex`、計算は`computations/verify_nine_closure_tangent.py`。接空間次元と新成分の次元を混同しない。

五変数族は七変数族の閉包に入ることが証明済みで、七変数族の軌道飽和閉包はClifford点で特異。証明 `sections/ore_boundary.tex`、厳密検算 `computations/verify_ore_boundary.py`。今後は閉包外の代数的変形と七・九変数族の成分最大性を調べる。

READMEの検証を実行し、成功したものだけをPASSにする。証明・ノート・方針・STATEを現行内容へ統合してmainへ通常fast-forwardでpushし、遠隔を読み戻す。PR・force push・GitHub Actions・無断のas4書込みをしない。PDFとキャッシュはGitで追跡しない。GitHubへのpushが成功した場合はLibraryに保存しない。

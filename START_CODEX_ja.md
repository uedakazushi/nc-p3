# Codex 起動指示 — nc-p3

`uedakazushi/nc-p3` の `main` 上で研究する。まず `research/plan_ja.tex` → `AGENTS.md` → `research/STATE.json` → `nc-p3.tex` → `research/notes_ja.tex` → `research/BASELINE.md` を読む。

**現在の具体的課題は Z4-AS-004。** 構造定理は既に原稿に統合したので、四窓群oidやFrobenius標準化の同じ整理を繰り返さずにKoszul条件・全AS正則性の正確な証明へ進む。

1. 各四窓$R_i\subset V_i\otimes V_{i+1}$から二次双対$E$の長さ0--5を計算し、Frobenius完全対合を有限階数条件で書く。既に原稿命題`prop:finite-frobenius`に証明あり。これ自体はAS正則性を保証しないことを守る。
2. Koszul複体の完全性を補う必要十分条件、あるいは十分性を厳密に証明できる被覆部分（PBWチャート、既存`as4`の認定済み族）を構築する。genericやdomain/PBWを全体へ仮定しない。各特殊値を扱う。
3. Frobenius標準4周期ではtrace積行列`P_i^T=Q_i`、周期`p=1,2`は四つの可逆行列`g_i`を使う有限関係空間同型式で判定できる（原稿`thm:finite-period-test`）。これらを既知2周期族または4周期例で実際に適用して、同型・例外パラメータの未分類域を縮める。
4. 2周期族の`A(λ)(1)=A(-λ)`と最小周期2での`A(λ)≄A(-λ)`を保つ。同型は添字保持であり、平行移動の商ではない。
5. 過去の`as4`周期性論文・検証スクリプトを固定blob SHAで参照し、重複実装を作らない。必要な実装のみnc-p3で維持する。数学的な証明と有限計算の証拠を別に記す。
6. 研究成果は適切な節の本文・研究ノート・研究方針・STATEへ整理し、最新版だけから再開できる状態にして`main`へ通常fast-forwardでpushする。PR、force push、GitHub Actions、無断のas4書込みはしない。GitにはTeXソースだけ、PDFは`latexmk`で`build/`へ。成功した検証だけをPASSとする。

現行の成果：幾何学的同型類でのFrobenius標準4周期同型の選択独立性と、4×4 trace行列による標準化・周期1/2/4の有限判定を証明済み。モジュライstack全体の同値と一般AS locusの必要十分条件は未完成。外積代数traceのPython有理数計算はPASSしたが、既存as4の全周期性コードや全稿組版を今回再実行した意味ではない。

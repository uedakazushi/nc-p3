# nc-p3 — 非可換射影三空間の分類

**2026-11-15 / v11。完全分類は未完成。正本は GitHub `main`。**

再開は [研究方針](research/plan_ja.tex) → [研究指示](AGENTS.md) → [研究ノート](research/notes_ja.tex) → [主論文](nc-p3.tex) の順。

## 現在の中心成果

標準14641型 quadratic Koszul AS-regular indexed Z-algebra を添字保持同型で分類する。graded同型、周期構造付き同型、qgr同値とは区別する。

**楕円Clifford部分代数を動かす九変数の代数的族**：主論文の `thm:algebraic-nine`（ソース `sections/elliptic_ore_family.tex`）。基底は `Spec k[u,v,r,s,a,b,c,t,h,(1+th)^(-1)]`。全ファイバーが標準AS正則であり、順序付き単項式の自由基底を持つ。同型による重複を除いた像の次元は9。非空なZariski開集合で最小周期2となり、どのgraded代数にも由来しない。形式級数の有限次数打切りではない。

三生成部分は `YW+WY=t Z^2`, `ZW+WZ=h Y^2`, `ZY+YZ=W^2`。その点曲線は一般に滑らかな楕円曲線となる。隣接する点配置を、楕円曲線上の12点の因子と非零2等分点による平行移動で比較して非1周期性を証明する。唯一の非被約点は長さ8である。

**七変数族はそのまま保持**：`t=h=0`への制限は `eq:seven-family` と表示そのものが一致する。旧七変数族の全パラメータPBW証明、像の次元7、一般の最小周期2の具体的な点配置証明も現行論文に残す。九変数族の一般点で二次Gröbner基底は主張しない。中心平方上のClifford基底とOre拡大を使う。

**他の保持する結果**：Frobenius標準4周期、四窓による幾何学的同型類、有限Frobenius検査、周期と同型の212変数388式。一般量子成分ではgraded同型類の次元6、indexed同型類の次元3。Clifford点では接空間が偶9＋交代6、二次障害は `Q=o0*o5+o1*o2-o2*o3+o3*o4`。

**Ore導分の完全分類と同型不変量**：固定した三生成Clifford部分上では、通常の次数一導分は内3＋外3の階数6、符号捩れ導分は中心二次式値の階数9。位置依存の基底変更を除き、一般Ore表示は九変数正規形へ移る（`sections/ore_complete.tex`）。楕円曲線の内在的な不変量は `J(x)=(x-2)^3*(x^3-6*x^2+228*x-8)^3/(x*(x-8)^3*(x+1)^6)`、`x=th`。一般のindexed同型なら `J(th)=J(t'h')` が必要で、同型類を固定すると `th` は高々12値。逆は未証明。

九変数族は基点で偶6・交代3方向を積分する。**9次元の像が全モジュライの既約成分であるとは未主張**。14次元のClifford局所模型、残る方向の代数的積分、全1周期構造、全成分分類は未解決。具体的な九変数族と先行研究との最終的な新規性比較は未了。

## 次の課題

**ALGEBRAIC-FAMILY-008**：Ore構成内部の次数一導分分類は完了。次は楕円曲線の12点因子と$j$障害を使って同型の必要十分条件を求め、別のClifford部分やOre以外の代数的構成で残る方向を積分する。成分としての最大性は未解決。別の三生成部分代数・Oreチャートで残る偶・交代方向を代数的に積分する。形式的高次障害の計算は候補選別の補助に使う。

## 再現

```sh
python3 -m pip install -r requirements.txt
python3 scripts/check_project.py
python3 computations/verify_exterior_trace.py
python3 computations/verify_period_witness.py
python3 computations/verify_moduli_jets.py --output research/moduli_certificate.json
python3 computations/verify_algebraic_family.py --output research/algebraic_family_certificate.json
python3 computations/verify_elliptic_family.py --output research/elliptic_family_certificate.json
python3 computations/verify_ore_derivations.py --output research/ore_derivation_certificate.json
python3 computations/verify_third_order_samples.py
latexmk nc-p3.tex
latexmk research/plan_ja.tex
latexmk research/notes_ja.tex
```

Python＋SymPyは厳密な有理数・多項式演算。三次サンプルは観察であり普遍的証明ではない。LuaLaTeX・latexmk・luatexjaを使う。主論文の組版には `sections/` 内の各TeXソースも必要。生成PDFはGitで追跡しない。証明・一次文献の所在と未確認点は [BASELINE](research/BASELINE.md)、実際の実行・保存状態は [validation](research/validation.json) を参照。

## 保存

`uedakazushi/nc-p3` の `main` に直接通常fast-forwardで保存する。PR・force push・GitHub Actionsは使わない。**GitHubへのpushが成功した場合はLibraryに保存しない**。失敗時のみ緊急退避を扱う。別リポジトリas4や既存Library資産を無断で変更・削除しない。

# nc-p3 — 非可換射影三空間の分類

**2026-11-15 / v15。完全分類は未完成。正本は GitHub `main`。**

再開は [研究方針](research/plan_ja.tex) → [研究指示](AGENTS.md) → [研究ノート](research/notes_ja.tex) → [主論文](nc-p3.tex) の順。

## 現在の中心成果

標準14641型 quadratic Koszul AS-regular indexed Z-algebra を添字保持同型で分類する。graded同型、周期構造付き同型、qgr同値とは区別する。

**楕円Clifford部分代数を動かす九変数の代数的族**：主論文の `thm:algebraic-nine`（ソース `sections/elliptic_ore_family.tex`）。基底は `Spec k[u,v,r,s,a,b,c,t,h,(1+th)^(-1)]`。全ファイバーが標準AS正則であり、順序付き単項式の自由基底を持つ。同型による重複を除いた像の次元は9。非空なZariski開集合で最小周期2となり、どのgraded代数にも由来しない。形式級数の有限次数打切りではない。

三生成部分は `YW+WY=t Z^2`, `ZW+WZ=h Y^2`, `ZY+YZ=W^2`。その点曲線は一般に滑らかな楕円曲線となる。隣接する点配置を、楕円曲線上の12点の因子と非零2等分点による平行移動で比較して非1周期性を証明する。唯一の非被約点は長さ8である。

**七変数族はそのまま保持**：`t=h=0`への制限は `eq:seven-family` と表示そのものが一致する。旧七変数族の全パラメータPBW証明、像の次元7、一般の最小周期2の具体的な点配置証明も現行論文に残す。九変数族の一般点で二次Gröbner基底は主張しない。中心平方上のClifford基底とOre拡大を使う。

**他の保持する結果**：Frobenius標準4周期、四窓による幾何学的同型類、有限Frobenius検査、周期と同型の212変数388式。一般量子成分ではgraded同型類の次元6、indexed同型類の次元3。Clifford点では接空間が偶9＋交代6、二次障害は `Q=o0*o5+o1*o2-o2*o3+o3*o4`。

**Ore型の完全性と同型不変量**：固定された三生成Clifford代数上の次数1の通常導分は自由階数6（内導分3＋外導分3）、符号捩れ導分は自由階数9（すべて中心二次式値）。位置に依存するXの基底変更を除き、一般Ore表示の係数は9個となり、九変数正規形はこの構成内で一般の表示を尽くす（`thm:ore-derivations`, `cor:ore-slice`（`sections/ore_complete.tex`））。

楕円曲線の内在的な$j$-不変量は$x=th$の次数12の有理関数
`J(x)=(x-2)^3*(x^3-6*x^2+228*x-8)^3/(x*(x-8)^3*(x+1)^6)`。
一般点でindexed同型なら`J(th)=J(t'h')`が必要で、固定した同型類に現れる$x$は高々12個（`prop:elliptic-j-obstruction`）。逆は未証明。

**楕円2-isogenyの内在的な同型障害**：一意の非被約点から三生成部分代数を復元し、任意のindexed同型は偏極付き楕円曲線の非零2等分点$\tau$を保存する。商曲線は $th A^3+t B^3+h C^3-(th+4)ABC=0$。$x=th$に対して$U(x)=4096(x+1)^3/(x(x-8)^3)$は次数4の有理関数で、一般の同型なら$U(th)=U(t'h')$が必要。以前の$j(E)$だけによる高々12候補を高々4候補へ改善した。$x\mapsto U$のGalois閉包の群は$A_4$。証明は`sections/isogeny_invariant.tex`。

**固定した三生成部分に関して非Ore型の別チャート**：`sections/nonore_family.tex`では$YW+WY=(-1)^i a_0 ZX$を含む五変数の二周期族を構成。全ファイバー標準14641型AS正則、indexed同型類の像は2次元。Clifford点での交代接方向$o_0$は従来の九変数族の接像に含まれない。さらに点スキームの支持を直線と五点に明示し、一般点の最小周期が正確に2であることを証明（`thm:nonore-generic-two`、`sections/nonore_period_two.tex`）。

**閉包包含と特異性（今回）**：五変数族は三生成部分 `Z,W,X` を選べばOre型になる。`YW+WY=(-1)^i a ZX+d(Z^2+X^2/v^2)` の平坦な一変数変形は、`d=0`で五変数族、`d≠0`で七変数族に同型。従って二次元の同型類像は七変数族の軌道飽和閉包に含まれる。七変数族の位置0の関係空間の37次元軌道飽和閉包はClifford点で接次元39以上となり特異。証明 `thm:nonore-in-seven-closure`, `thm:seven-closure-singular`（`sections/ore_boundary.tex`）、検算 `computations/verify_ore_boundary.py`。全モジュライの成分性はなお未解決。

**九変数族の軌道閉包の特異性（v15）**：位置0の関係空間の軌道飽和閉包 `Z9` は既約39次元で、`Z7 ⊂ Z9` は余次元2。Clifford点の接空間は41次元以上で、従って特異。九変数族自身の偶6＋交代3方向に、七変数族の境界退化由来の偶1・交代1方向を足した11方向が、個別の代数的一変数族の接方向として独立。ただし既知の二次障害はこの線型包上 `lambda*c` となるため、滑らかな11次元の同時変形族は構成できない。証明 `thm:nine-closure-singular`, `prop:nine-closure-obstruction`（`sections/nine_closure.tex`）、検証 `computations/verify_nine_closure_tangent.py`。

九変数族は基点で偶6・交代3方向を積分する。**9次元の像が全モジュライの既約成分であるとは未主張**。14次元のClifford局所模型、残る方向の代数的積分、全1周期構造、全成分分類は未解決。具体的な九変数族と先行研究との最終的な新規性比較は未了。

## 次の課題

**ALGEBRAIC-FAMILY-008**：九変数族の2-isogenyと12点因子から完全な同型判定を求める。非Ore五変数族の一般最小周期2は解決。七変数族への閉包包含は解決済み。九変数族の軌道飽和閉包も特異であると分かったので、二次障害の許す未実現方向に対し既存閉包外の代数的変形を探し、七・九変数族の成分最大性と全モジュライの既約成分を判定する。全モジュライの既約成分と14次元Clifford局所模型は未決定。

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
python3 computations/verify_isogeny_invariant.py --output research/isogeny_certificate.json
python3 computations/verify_nonore_family.py --output research/nonore_certificate.json
python3 computations/verify_nonore_point_geometry.py --output research/nonore_point_certificate.json
python3 computations/verify_ore_boundary.py --output research/ore_boundary_certificate.json
python3 computations/verify_nine_closure_tangent.py --output research/nine_closure_certificate.json
python3 computations/verify_third_order_samples.py
latexmk nc-p3.tex
latexmk research/plan_ja.tex
latexmk research/notes_ja.tex
```

Python＋SymPyは厳密な有理数・多項式演算。三次サンプルは観察であり普遍的証明ではない。LuaLaTeX・latexmk・luatexjaを使う。主論文の組版には `sections/` 内の各TeXソースも必要。生成PDFはGitで追跡しない。証明・一次文献の所在と未確認点は [BASELINE](research/BASELINE.md)、実際の実行・保存状態は [validation](research/validation.json) を参照。

## 保存

`uedakazushi/nc-p3` の `main` に直接通常fast-forwardで保存する。PR・force push・GitHub Actionsは使わない。**GitHubへのpushが成功した場合はLibraryに保存しない**。失敗時のみ緊急退避を扱う。別リポジトリas4や既存Library資産を無断で変更・削除しない。

# nc-p3 — 非可換射影三空間の周期構造と代数的変形

**2026-11-22 / v40。標準14641型の完全分類は未完成。正本は GitHub main。**

## 最初に読む

[研究方針](research/plan_ja.tex) → [主論文](nc-p3.tex) → [必要な研究ノート](research/notes_ja.tex)。
Cloud Codexの環境準備・権限・止まりにくい研究ループ・再開の手順は [Cloud Codexセットアップ](docs/cloud_codex_ja.md) と
`bash scripts/cloud_preflight.sh` を参照する。
実行規則は [AGENTS.md](AGENTS.md)、終盤を見据えた分類の方法・nc-p3自体の巨大化に対する三層アクセスと差分検証は [分類アトラス設計](research/classification_architecture.md)、
一次文献と外部依存は [BASELINE.md](research/BASELINE.md)（as4は必要箇所だけ取得）、
主張と課題の索引は [STATE.json](research/STATE.json)、今回の検証範囲は [validation.json](research/validation.json) を参照する。

対象は、標数零の代数閉体上の標準14641型 quadratic Koszul AS-regular indexed Z-algebra。
同値関係は**添字保持同型**であり、graded同型、周期構造を指定した同型、qgr同値とは区別する。

## 成果と限界

| 主題 | 原稿に証明がある結果 | まだ結論していないこと |
| --- | --- | --- |
| 有限四窓 | Frobenius標準4周期、幾何学的同型類の四窓への還元、認証済みAS対象の同型・周期の有限方程式 | 四次の必要階数条件と五次までの正確性の二階数判定は証明済み。一般四窓のKoszul/AS条件の有限な必要十分判定、明示的な正則族による全被覆は未解決 |
| 量子多項式 | 一般の局所成分でgraded同型類6次元、indexed同型類3次元 | 非1周期的対象が増える現象ではない。三つのtwist方向を忘れる比較 |
| Clifford点 | 選んだgraded表示の接空間9、indexed変形の接空間15。全15変数の二次障害は、六つの交代変数上で非退化な二次式 Q 一つ | 成分次元14、高次障害の消滅、全1周期構造の分類 |
| 楕円Clifford–Ore族 | 九変数の実際の代数的族。全ファイバーAS正則、同型類の像9次元、一般の最小周期2 | 九次元像が既約成分であること、全特殊値の周期、完全な族内同型判定 |
| 同型不変量 | 内在的な偏極付き楕円曲線と非零2等分点。次数4の必要不変量 U(th) | Uの一致による十分な同型判定 |
| 境界と特異性 | 五変数の別のOreチャートは同型類像2次元、一般最小周期2。七変数族の軌道閉包に含まれる。閉包 Z7, Z9 はClifford点で特異。指定11方向の滑らかな同時積分は不可能 | 独立した新成分、全モジュライの特異性 |

「原稿に証明があること」「今回の論理監査」「厳密計算の再実行」「形式化・査読・新規性確認」は別の状態である。
新規性と第三者査読は未確定。具体的な検証範囲は validation を参照する。

### 中心となる代数的族

九変数族の基底は
\(\operatorname{Spec}k[u,v,r,s,a,b,c,t,h,(1+th)^{-1}]\)。
三生成部分は
\(YW+WY=tZ^2,\ ZW+WZ=hY^2,\ ZY+YZ=W^2\)。
全関係、自由基底、正則性、次元、非1周期性の証明を
[elliptic_ore_family.tex](sections/elliptic_ore_family.tex) に置く。

t=h=0 は [七変数族](sections/algebraic_family.tex) そのものであり、その全パラメータPBW証明と点配置の証明を保持する。
一般の t,h では二次Gröbner基底を仮定せず、Clifford基底とOre拡大を使う。1+th != 0 は不可欠である。

九変数族の全対象は、交代係数を (a,b,c) → (λa,λb,λc) と動かすことで、
全ファイバーAS正則な A1 族の中で1周期的対象へ特殊化する（cor:nine-to-period-one）。
一般の最小周期2は、1周期的対象との変形連結性を否定しない。
標準14641型で1周期的 locus と変形連結に分離した例の存在は未解決である。

nonore_* という既存ファイル名・ラベルは参照の安定性のため保持するが、その意味は**指定した三生成部分に関して非Ore**というものに限る。
[別の三生成部分によるOre表示と閉包包含](sections/ore_boundary.tex) は既に証明されている。

最上次数テンソルによる [Frobenius候補の復元](sections/frobenius_tensor.tex) と、
第三の行列階数による全対象の二分岐を証明した。
行列束の全場合分けにより階数2を除外し、全対象の階数範囲を3以上16以下へ絞った。
[枝台帳](research/classification_branches.json) の両枝の一般正規形・正則性判定は未完成。
量子部分では階数6・8・10・12の全条件を求めた。さらに [量子部分の全indexed分類](sections/quantum_indexed.tex) として、
全ての非零パラメータを三次元トーラスの明示的なS4作用へ還元した。特殊値も含み、階数6は二種類となる。
対称テンソル部分は四二次式の完全交叉によるチャートで全ファイバー正則性を認証し、
その全特殊値を含むindexed同型を33未知数41方程式へ還元した。
このチャートへの所属は三つの線型解空間内の可逆行列の存在で必要十分に判定できる。
符号を変えた線型判定は可換多項式のindexed類を認識する。有限Frobenius候補では一つの可逆な負符号証明書で十分である。残る低階数領域の非空性・排除は未解決。

可換三生成部分上の [18係数Oreチャート](sections/polynomial_ore.tex) では、
全ファイバーの正則性と `rho = 6 + rank Q` を証明した。`Q` は6×6交代行列で、
低階数部分は可換多項式のindexed類だけとなる。18は係数空間の次元である。

位置依存の可換・反可換Ore四窓でも、三次の線型条件が全AS/Koszul正則性に必要十分となる
（[可換側](sections/polynomial_ore.tex)、[反可換側](sections/anticommutative_ore.tex)）。
反可換側の三つの消せない交代方向は階数8となり、両チャートの低階数部分は既知の二枝へ入る。

さらに [対角量子Oreチャート全体](sections/quantum_ore.tex) でも30本の線型導分条件が正則性に必要十分。
量子背景への特殊化と上の正規化から、低階数なら可換型または対称型に限られ、
未分類低階数領域はこの表示を持たないことを証明した。

[Clifford点の低階数近傍](sections/clifford_rank_slice.tex) も対称チャートで被覆した。
階数6以下の有限Frobenius候補表示空間はこの点で滑らかな69次元となり、
既知の二次障害Qは交代一次変形の4×4行列のPfaffianとして現れる。
69は表示空間の次元で、全Clifford成分や粗モジュライの次元ではない。

全反可換点でも対称チャートの局所被覆を証明し、可換多項式の軌道は低階数表示空間で開となる。
三つの行列階数で切った**基底付きテンソル表示空間**には、69次元と45次元の二つの既約成分を同定した。
他の成分や境界は未分類であり、全ASモジュライの成分とは主張しない。

### 次の研究

**最優先は全四窓のAS正則条件・構造的不変量から網羅的な分類木を構築すること。**
正規形の枝、特殊値、枝間同型の条件を証明し、未分類の残余を明示する。
既存の九変数族・量子成分・Clifford点・境界族は、全分類を目指す方法の回帰試験として使う。
九変数族の最大性・完全族内同型やClifford点の高次障害は有用な補助課題だが、
全分類の前提条件にはしない。詳細は研究方針と分類アトラス設計を参照する。
機械可読な巨大分類台帳と自動生成器は、最初の網羅的分岐の証明に合わせ最小構成で導入する。

## 再現

Python 3 と SymPy 1.14.0。LuaLaTeX・latexmk・luatexja・Harano Ajiフォントを用いる。
すべてリポジトリのルートから実行する。

~~~sh
python3 -m pip install -r requirements.txt
python3 scripts/check_project.py
python3 computations/verify_exterior_trace.py
python3 computations/verify_period_witness.py
python3 computations/verify_quartic_obstruction.py
python3 computations/verify_frobenius_tensor.py --output research/frobenius_tensor_certificate.json
python3 computations/verify_quantum_tensor_rank.py
python3 computations/verify_symmetric_congruence.py
python3 computations/verify_rank_two_exclusion.py
python3 computations/verify_quantum_indexed.py
python3 computations/verify_polynomial_ore.py --output research/polynomial_ore_certificate.json
python3 computations/verify_indexed_polynomial_ore.py --output research/indexed_polynomial_ore_certificate.json
python3 computations/verify_indexed_anticommutative_ore.py --output research/indexed_anticommutative_ore_certificate.json
python3 computations/verify_quantum_ore.py --output research/quantum_ore_certificate.json
python3 computations/verify_clifford_rank_slice.py --output research/clifford_rank_slice_certificate.json
python3 computations/verify_pure_anticommutative_slice.py --output research/pure_anticommutative_slice_certificate.json
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
~~~

check_project.py は参照・入力・JSONの構造検査であり数学的証明ではない。
最後の三次サンプル検算は三方向の観察に限る。PDFは build/ に生成しGit追跡しない。
文書・コード・入力・証明書をGitHub mainから再生成する。通常のfast-forwardで直接保存し、PR・force push・GitHub Actionsは使用しない。
GitHub保存成功時はLibraryへ重複保存しない。

対称チャート全体の表示像は `PGL4^3 × U` と同型な滑らかな69次元局所閉部分である。
四二次式完全交叉のテンソルcentroidが常にスカラーであることから、一意な三つの相対フレームを回収できる。
個々のindexed同型類は正則領域内でも非閉になり得る明示例を得た。
対称表示像の和集合は有限Frobenius領域内でも相対的に閉ではない。
[明示的な非Koszul境界](sections/nonkoszul_boundary.tex)は全ての四次検査を通るが、
ある始点の五次元数が56ではなく63となる。五次の7次元Koszul障害は、平方関係による6次元と明示的な7項の混合サイクルによる1次元で完全に生成され、3座標の線型汎関数で後者の独立性を検出できる（`prop:boundary-complete-quintic-cycles`）。一般の基点三因子対称境界も、左右反転対合によるTor下界5から非Koszulと証明した。対称閉包の全境界被覆は未解決である。

三因子対称な有限Frobenius候補を、四層の切断代数のtrivial extensionとして記述した。
また一般点の有限次数障害を族全体の排除へ拡張する階数補題を得た。
二変数係数族では一般点の五次階数が77・36であることを有理関数体上で厳密計算した。
半連続性により、非零係数の特殊値も含めて標準型Koszul AS点を一様に排除した。
これは全ての一基点四二次式系の排除ではない。

二変数境界の厳密計算には追加の多項式演算バックエンドを使う：

```bash
python3 -m pip install -r requirements-symbolic.txt
python3 computations/verify_generic_one_point_boundary.py --output research/generic_one_point_boundary_certificate.json
```

係数は有理関数体 `QQ(u,v)` の元として扱い、数値標本や補間で一般点を認定しない。

交差項を含む、固定Boolean完全交叉上の全一基点二次式には、五係数域の全体で正則な評価フレームを構成した。
可換商のHilbert列は全域で `(1,4,6,4,1,1,...)`。Frobenius条件を通る点は対称チャート外にある。
この五係数域の一般五次$d_2$階数も77と厳密計算し、全域で標準型Koszul AS点を排除した。

```bash
python3 computations/verify_boolean_quadric_frame.py --output research/boolean_quadric_frame_certificate.json
python3 computations/verify_boolean_quadric_homology.py --output research/boolean_quadric_homology_certificate.json
```

後者は五変数の厳密消去で、上の二変数検算より重い。

相異なる八点の任意の三二次式完全交叉にも、一基点境界の評価フレームを拡張した。
点配置は列スカラーと射影変換の後に `[I | M]`, `M M^T = I` と表示できる（古典的な自己付随性に対応）。
一般位置の有理数標本でも、有限Frobenius・四次正確・五次障害7を検算した。
一般の八点配置全体も、左右反転対合による普遍的Tor下界によりAS対象から排除された。

```bash
python3 computations/verify_reduced_octad_frame.py --output research/reduced_octad_frame_certificate.json
```

基点を固定した四二次式系 `Gr(4,9)` の一般点を、積型の三二次式と第四式からなる8係数族と基底変換で被覆した。
接写像の20次小行列式 `-696/78125` により支配性を証明した。これはAS正則性や同型類次元の判定ではない。
八変数の五次階数が84未満なら、フレーム束上の半連続性で基点域全体を排除できる移送命題も得た。
八変数の正確な階数は未確定だが、一般の基点域全体のAS排除は独立の構造的証明で完成した。
一方、長さ2の基点部分スキームは構造的に排除した。相異なる二点も二重点も四次元数42以上を強制する。
四次検査を通る三因子対称な基点候補は一縮約点だけを持ち、縮約八点netによる表示と可換最小分解が適用できる。

左右反転対合により、戻り辺一回の商複体の負固有部分は次元 `60,200,135` を持つ。
Euler標数 `-5` が非線型 `Tor_(4,5)` を強制するため、基点を持つ三因子対称候補はすべて非Koszulである。
したがって、認証済みAS対象の対称枝所属は二つの可逆な正符号証明書で判定できる。

```bash
python3 computations/verify_product_net_dominance.py --output research/product_net_dominance_certificate.json
python3 computations/verify_product_net_homology.py --output research/product_net_homology_certificate.json
# One-return bimodule Tor: checked positive control and exact two-variable case
python3 computations/verify_trivial_extension_bitor.py --field ci --output research/trivial_extension_bitor_ci_certificate.json
python3 computations/verify_trivial_extension_bitor.py --field two --output research/trivial_extension_bitor_two_certificate.json
# Point/jet quartic obstruction and exact point-quotient regression
python3 computations/verify_basepoint_length_two.py --output research/basepoint_length_two_certificate.json
python3 computations/verify_trivial_extension_bitor.py --field two --frame ternary --point-quotient --last-only --output research/trivial_extension_bitor_quotient_certificate.json
# Signed-reversal theorem regression, also at a non-Boolean product-net point
python3 computations/verify_one_return_involution.py --output research/one_return_involution_certificate.json
python3 computations/verify_one_return_involution.py --fixture generic --output research/one_return_involution_generic_certificate.json
# All cuts through degree seven, exact differentials and two controls
python3 computations/verify_cut_tor_dimensions.py
python3 computations/verify_general_cut_tor.py --output research/general_cut_tor_certificate.json
# Control for the residual single-adjacent-negative linear matrix slice
python3 computations/verify_single_negative_slice.py --output research/single_negative_slice_certificate.json
# Exhaustive partial-alternation models and their Frobenius exclusions
python3 computations/verify_negative_witness_rigidity.py --output research/negative_witness_rigidity_certificate.json
# Optional exact eight-variable rank; exclusion theorem no longer depends on it
python3 computations/verify_trivial_extension_bitor.py --field eight --frame ternary --point-quotient --last-only --composition-check structural --output research/trivial_extension_bitor_eight_certificate.json
```

任意の有限Frobenius標準四窓も、一辺切断した32次元四層代数のtrivial extensionとして表せる。
四切断の戻り辺一回Torは内部次数7以下を覆い、最大の次元列は120,1040,2400,1600となる。
符号付きOre対照点とBoolean境界点で、七次までの全微分の階数とホモロジーを有理数上で検算した。
連続二窓の三次交叉が四次元のとき、そこから作る56次元代数が3-representation infiniteであることは標準AS延長の必要十分条件である。延長の同型類はこの二窓で決まる。この全次数条件を有限個の反復で判定する方法と、二窓の明示的な軌道分類は未解決である。

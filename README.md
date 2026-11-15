# nc-p3 — 非可換射影三空間の周期構造と代数的変形

**2026-11-16 / v22。標準14641型の完全分類は未完成。正本は GitHub main。**

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

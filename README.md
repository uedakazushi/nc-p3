# nc-p3 — 非可換射影三空間の周期構造と代数的変形

**2026-11-15 / v17。標準14641型の完全分類は未完成。正本は GitHub main。**

## 最初に読む

[研究方針](research/plan_ja.tex) → [主論文](nc-p3.tex) → [必要な研究ノート](research/notes_ja.tex)。
実行規則は [AGENTS.md](AGENTS.md)、一次文献と外部依存は [BASELINE.md](research/BASELINE.md)、
主張と課題の索引は [STATE.json](research/STATE.json)、今回の検証範囲は [validation.json](research/validation.json) を参照する。

対象は、標数零の代数閉体上の標準14641型 quadratic Koszul AS-regular indexed Z-algebra。
同値関係は**添字保持同型**であり、graded同型、周期構造を指定した同型、qgr同値とは区別する。

## 成果と限界

| 主題 | 原稿に証明がある結果 | まだ結論していないこと |
| --- | --- | --- |
| 有限四窓 | Frobenius標準4周期、幾何学的同型類の四窓への還元、認証済みAS対象の同型・周期の有限方程式 | 一般四窓のKoszul/AS条件の有限な必要十分判定、明示的な正則族による全被覆 |
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

### 次の研究

最優先は九変数族の**一般点**の変形空間を調べ、その像が成分を尽くすかを決めること。
同時に、楕円曲線・2等分点・12点因子から族内同型の必要十分条件を求める。
不足方向があれば実際の代数的族で積分する。Clifford特殊点の接空間だけから成分次元を推測しない。
詳細な着手手順と、全1周期構造・全四窓正則locusへの接続は研究方針に集約する。

## 再現

Python 3 と SymPy 1.14.0。LuaLaTeX・latexmk・luatexja・Harano Ajiフォントを用いる。
すべてリポジトリのルートから実行する。

~~~sh
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
~~~

check_project.py は参照・入力・JSONの構造検査であり数学的証明ではない。
最後の三次サンプル検算は三方向の観察に限る。PDFは build/ に生成しGit追跡しない。
文書・コード・入力・証明書をGitHub mainから再生成する。通常のfast-forwardで直接保存し、PR・force push・GitHub Actionsは使用しない。
GitHub保存成功時はLibraryへ重複保存しない。

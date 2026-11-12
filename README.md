# ASZ4 — 4次元 quadratic AS-regular Z-algebra の完全分類

**2026-11-12 / v3・研究継続。分類未完。GitHub `main` を正本とする。**

再開時は [研究方針](research/plan_ja.tex) → [Codex指示](AGENTS.md) → [研究ノート](research/notes_ja.tex) → [分類原稿](papers/asz4_classification.tex) の順に読む。次の最優先課題は **Z4-AS-004**（Koszul/AS正則な四窓locusの分類）。

## 目標

標数0の代数閉体上、標準14641型の局所有限・正連結・quadratic・Koszul、左右AS正則なZ-algebra `C` を**添字を保つZ-algebra同型**で完全分類する。全対象被覆・各閉パラメータの正則性・同型の必要十分条件・周期1/2/4・変形/既約成分の適切な記述を達成する。domain / noetherian / PBW は全体の追加仮定にしない。一般の「4次元 quadratic AS-regular」が標準14641型から外れる可能性は研究範囲の境界として明記する。

既知の一般4周期性と最小周期1/2/4の実現は **[as4 の周期性論文](https://github.com/uedakazushi/as4/blob/main/research/z_periodicity/papers/as4_z_periodicity_ja.tex)** に証明・検証資産がある。未査読・新規性未確認と区別し、同じ計算を一から組み直さない。

## 既存研究と新研究の分担

| 分野 | 基本の再利用元 | asz4で解く問題 |
|---|---|---|
| 1周期的な部分 | [as4](https://github.com/uedakazushi/as4) graded分類（未完） | Zhang twist同値、周期構造のファイバー、Z側での像と特殊化 |
| 全14641型の4周期性 | [as4周期性論文](https://github.com/uedakazushi/as4/tree/main/research/z_periodicity) | 4周期の同型を**選んだ**有限atlasと、選択を忘れた分類の差 |
| 最小周期2/4の例 | 上記のPBW反例・厳密検証コード | 開/閉境界、族内外同型、AS条件、成分/変形 |
| 一般の4周期的対象 | Bondal–Polishchuk、Mori–Nyman | 四つの生成空間と関係空間、Koszul双対・量子行列式の有限条件、全被覆 |

既存 `as4` の全成果や数GB規模のリポジトリをコピーしない。必要な定理は固定Git SHA・仮定・証明箇所を [再利用台帳](research/BASELINE.md) に記録し、asz4固有の数学を現行の原稿・ノートに置く。

## 今回の確定した数学的進展

**Frobenius標準4周期同型による分類の簡約**：標準14641型のKoszul双対の局所Frobenius traceから得る4周期同型は、traceの取り方を変えても頂点スカラー自己同型で共役になる。従って**幾何学的点の同型類**では、裸のZ代数はFrobenius標準4周期の四窓の \(\mathrm{GL}_4^4\)-軌道と一対一に対応する（定理 `thm:canonical-four`、系 `cor:canonical-slice`）。モジュライstackの同値は未主張。

標準化は、二次双対の1×3と3×1積から作る四つの4×4行列について \(P_i^\mathsf{T}=Q_i\) を満たすtraceが存在するか、という**有限条件**になる。二次双対の局所Frobenius性も次数0–5の行列階数で検査できる。ただしKoszul性は別条件。

さらにFrobenius標準四窓では、裸のCの最小周期1/2/4を、四つの可逆写像 \(g_i:V_i\to V_{i+p}\) が \((g_i\otimes g_{i+1})R_i=R_{i+p}\) を満たすかという有限可解性で完全に判定できる（定理 `thm:finite-period-test`）。

周期の存在はさらに、4個の\(16\times6\)関係行列と未知\(4\times4\)行列・補助\(6\times6\)行列による **212変数・388式の有限多項式系**に書けます。候補を厳密な有理数階数で検証する[証明書検証器](computations/verify_period_witness.py)も追加・実行済みです。

既存as4の2周期族 \(A(\lambda)\) では \(A(\lambda)(1)=A(-\lambda)\) だが、最小周期2では \(A(\lambda)\not\simeq A(-\lambda)\)（添字保持同型）。添字平行移動を誤ってmoduli同型に読み替えない。標準例の外積代数のtrace符号規格化は [正確な検証コード](computations/verify_exterior_trace.py) で確認済み。

**証明全文**：[分類論文](papers/asz4_classification.tex)（v5）。既存の最小周期2/4と非1周期的平坦族は引き続き [as4の周期性研究](https://github.com/uedakazushi/as4/tree/main/research/z_periodicity) を参照し、重複実装を避ける。

## 残る核心：Z4-AS-004

四窓の二次関係・有限Frobenius性・Frobenius標準化・裸の周期判定の構造は証明済み。次は**その四窓がKoszul AS正則であるための必要十分条件**を、PBWを全域に仮定せず全特殊値で解く。まず既存の正則性認定済みの族を標準四窓へ移し、exactnessの証明と有限の関係テンソル条件の接続を完成させる。point-schemeなどの有限長不変量は反例検出に用いるが、同型の十分条件とはしない。

**完了判定**：四窓上のAS正則条件の必要十分性、全特殊値の処理、2/4周期族での具体的な判定、そこから全対象被覆へ繋がる証明。全分類はなお未完成。

## 実行・保存

```sh
python3 scripts/check_project.py
latexmk papers/asz4_classification.tex
latexmk research/plan_ja.tex
latexmk research/notes_ja.tex
python3 computations/verify_exterior_trace.py
python3 computations/verify_period_witness.py
```

TeX Live / LuaLaTeX / latexmk / luatexja を使用し、PDF・auxは追跡せず `build/` に出力。GitHub Actionsは使わない。個々の検算は内容を説明したコードで実行する。GitHub `main` に直接fast-forwardで保存し、push成功時にLibraryへ重複保存しない。GitHubへの書込みが失敗した場合のみLibraryへ緊急退避する。

論文原稿は研究の経緯ではなく、現在有効な主張と証明だけから構成する。公表時は重要な文献を一次資料で検証する。

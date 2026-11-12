# ASZ4 — 4次元 quadratic AS-regular Z-algebra の完全分類

**2026-11-12 / v1・研究開始。分類未完。GitHub `main` を正本とする。**

再開時は [研究方針](research/plan_ja.tex) → [Codex指示](AGENTS.md) → [研究ノート](research/notes_ja.tex) → [分類原稿](papers/asz4_classification.tex) の順に読む。初動の現行課題は **Z4-WINDOW-001**。

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

## 最初の分類課題：Z4-WINDOW-001

正規化された**4周期表示を選んだ**二次Z代数を、4頂点・各4本の隣接生成元、各6次元の二次関係空間から作る。対応する `GL_4^4` の作用と有限構成を証明し、AS正則・Koszulの全条件を有限データ上の正確な証明義務として定式化する。さらに、選んだ4周期同型を固定する同型と、裸のZ代数の同型を厳密に区別する。point-schemeなどの有限長不変量は反例検出に用いるが、同型の十分条件とはしない。

**完了判定**：有限表示・ゲージ同値の証明、periodicity判定の正確な義務、graded側への忘却射の図式、既知周期2/4例への適用を揃える。これだけでは全分類完了ではない。

## 実行・保存

```sh
python3 scripts/check_project.py
latexmk papers/asz4_classification.tex
latexmk research/plan_ja.tex
latexmk research/notes_ja.tex
```

TeX Live / LuaLaTeX / latexmk / luatexja を使用し、PDF・auxは追跡せず `build/` に出力。GitHub Actionsは使わない。個々の検算は内容を説明したコードで実行する。GitHub `main` に直接fast-forwardで保存し、push成功時にLibraryへ重複保存しない。GitHubへの書込みが失敗した場合のみLibraryへ緊急退避する。

論文原稿は研究の経緯ではなく、現在有効な主張と証明だけから構成する。公表時は重要な文献を一次資料で検証する。

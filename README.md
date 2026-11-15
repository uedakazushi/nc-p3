# nc-p3 — 非可換射影三空間の分類

**2026-11-15 / v9。完全分類は未完成。正本は GitHub `main`。**

再開は [研究方針](research/plan_ja.tex) → [研究指示](AGENTS.md) → [研究ノート](research/notes_ja.tex) → [主論文](nc-p3.tex) の順。

## 現在の中心成果

標準14641型の quadratic Koszul AS-regular indexed Z-algebra を添字保持同型で分類する。graded同型、周期構造付き同型、qgr同値とは区別する。

**実際の代数的七変数族**：`nc-p3.tex` の `eq:seven-family` に多項式関係を記載した。全 affine 7-space 上で自由なPBW族であり、全ファイバーが標準AS正則。非空なZariski開集合で最小周期2となり、どのgraded代数にも由来しない。同型による重複を除いた像の次元は7。三つの独立な交代方向を同時に代数的に積分する。

**保持する結果**：Frobenius標準4周期、四窓による幾何学的同型類、有限Frobenius検査、周期と同型の212変数388式。一般の量子多項式成分ではgraded同型類の次元6、indexed同型類の次元3。Clifford点では接空間が偶9＋交代6、二次障害は `Q=o0*o5+o1*o2-o2*o3+o3*o4`。

七変数族は交代空間の最大等方3平面 `(0,-b,-a,-b,0,c)` を積分する。**この7次元の像が全モジュライの既約成分であるとは未主張**。14次元のClifford局所模型、残る方向の代数的積分、全周期1構造、全成分分類は未解決。

## 次の課題

**ALGEBRAIC-FAMILY-008**：形式的高次障害の計算より、明示的な代数的族の拡張を優先する。三生成部分代数Bを動かす、別のBや重み順序を使う、既存七変数族との同型・包含・閉包を決める。全パラメータの正則性と非1周期性を別々に証明する。

## 再現

```sh
python3 -m pip install -r requirements.txt
python3 scripts/check_project.py
python3 computations/verify_exterior_trace.py
python3 computations/verify_period_witness.py
python3 computations/verify_moduli_jets.py --output research/moduli_certificate.json
python3 computations/verify_algebraic_family.py --output research/algebraic_family_certificate.json
python3 computations/verify_third_order_samples.py
latexmk nc-p3.tex
latexmk research/plan_ja.tex
latexmk research/notes_ja.tex
```

Python＋SymPyの計算は厳密な有理数・多項式演算。三次サンプルは観察であり普遍的証明ではない。LuaLaTeX・latexmk・luatexjaを使用する。生成PDFはGitで追跡しない。証明と計算の所在、未検証の新規性は [BASELINE](research/BASELINE.md)、現在の実行・保存状態は [validation](research/validation.json) を参照。

## 保存

`uedakazushi/nc-p3` の `main` に直接通常fast-forwardで保存する。PR・force push・GitHub Actionsは使わない。**GitHubへのpushが成功した場合はLibraryに保存しない**。失敗時のみ緊急退避を扱う。別リポジトリas4や既存Library資産を無断で変更・削除しない。

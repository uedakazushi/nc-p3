# nc-p3 再開指示 — 2026-11-15 / v9

まず research/plan_ja.tex、AGENTS.md、research/notes_ja.tex、nc-p3.tex、research/STATE.json を読む。

最優先は **ALGEBRAIC-FAMILY-008**。形式的変形だけでなく、係数を代数的パラメータで書ける実際の族を求める。既に七変数族が全ファイバーでAS正則、一般点で最小周期2、同型類の像で7次元と証明された。computations/verify_algebraic_family.py を再現し、三生成部分代数Bの変形や別チャートから追加の偶・交代方向を構成する。

同型類の族の次元7と既約成分の次元を区別する。Clifford点の接次元15、二次障害Q、14次元の局所模型予想はそのまま保持するが、高次の証明済み結果と扱わない。原稿にある有限四窓・周期・同型判定・量子成分比較の証明を削らない。

成果は現行の原稿・方針・ノート・再現コードへ統合して main に通常のfast-forwardで直接保存する。PR、force push、GitHub Actionsは使わない。PDFは追跡しない。GitHub push成功時にはLibraryへ保存しない。遠隔保存に失敗した場合だけ緊急退避を検討する。

# Codex 起動指示 — ASZ4

`uedakazushi/asz4` の `main` 上で研究する。まず `research/plan_ja.tex` を読み、`AGENTS.md`、`research/STATE.json`、`research/notes_ja.tex`、`research/BASELINE.md`、`papers/asz4_classification.tex` を読む。

次の具体的課題は **Z4-WINDOW-001**。既に証明を書いた「4周期同型を選んだ二次Z代数の四窓とGL4^4作用群oid」の先へ進む。

1. 四つの関係空間 `R_i ⊂ V_i⊗V_{i+1}` の上で、標準14641型AS正則性に必要かつ十分なKoszul解消・双対完全対合の条件を整理する。PBWやgeneric仮定は全対象に課さず、全閉値の正確な適用域を守る。
2. `as4/research/z_periodicity/` の既存2周期・4周期例を、この有限四窓へ表示する。実証明・実検算コードをGit SHAで参照し、不要な再実装をしない。
3. `Csimeq C(4)` を選ぶことと、4周期同型を忘れることの差を明確化し、`(C,τ_4)` の分類結果を裸の `C` へ移す際の残る義務を1つでも数学的に減らす。
4. 新たな証明・反例・数式・検証を現在の論文・ノートへ直接統合する。単なるタスク更新や記録整理だけで作業を終えない。証明済み/未監査/計算だけの状態を区別する。
5. `python3 scripts/check_project.py` と対象の計算、利用可能なら `latexmk papers/asz4_classification.tex` を実行。成功/未実行/失敗を明記し `main` へ通常fast-forwardで保存する。PR・GitHub Actions・force pushは使わない。

新規な結果が得られた場合だけ `research/STATE.json` に数学的進捗のrevisionを反映する。現在の正本・現在の証明だけを残し、バックアップや作業履歴を積み上げない。外部`as4`に勝手にcommitしない。

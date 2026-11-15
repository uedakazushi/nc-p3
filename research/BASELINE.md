# nc-p3の既存成果・先行研究への参照（2026-11-12 / v1）

目的：必要な定理・証明位置・根拠・未検証を把握し、as4を丸ごと複製せずに現行nc-p3を再開する。as4の全分類は完成していない。以下は固定されたGit **blob SHA** を持つ資料への参照であり、as4側の内容・現行revisionをASZ4の成果として再宣言しない。

## as4の一次的研究成果（私的研究原稿）

| 対象 | パスと内容 | 2026-11-12確認blob SHA | nc-p3での使用 |
|---|---|---|---|
| 現状と研究方針 | [as4/README.md](https://github.com/uedakazushi/as4/blob/main/README.md) | `b19baa4a10997d53de8261ff7a502a018f835af2` | graded側の完成範囲・課題 |
| graded 4D分類の研究方針 | [as4/research/plan_ja.tex](https://github.com/uedakazushi/as4/blob/main/research/plan_ja.tex) | `b6192a518477fda5963af6d878c18d8bb065a584` | graded同型とZhang twistの分離 |
| 被覆表 | [as4/research/COVERAGE.md](https://github.com/uedakazushi/as4/blob/main/research/COVERAGE.md) | `f4c9a698fb5503838bd71759328829f599c636ad` | Q/V、CY、量子、多重根、共鳴の局所成果 |
| 単冪Zhang untwist | [as4/research/proofs/UNIPOTENT_UNTWIST.tex](https://github.com/uedakazushi/as4/blob/main/research/proofs/UNIPOTENT_UNTWIST.tex) | `eaa703cb8c6c68c3b0f8d1105cf240ee52e47bba` | `A=B^H` でZ同型が同一となる点の検証 |
| 周期1/2/4の論文 | [as4/research/z_periodicity/papers/as4_z_periodicity_ja.tex](https://github.com/uedakazushi/as4/blob/main/research/z_periodicity/papers/as4_z_periodicity_ja.tex) | `cc475fae425a23a02fd8559e403ad2ec375302f3` | 標準14641での4周期、2/4例、平坦族 |
| 周期性検証の再現 | [as4/research/z_periodicity/README.md](https://github.com/uedakazushi/as4/blob/main/research/z_periodicity/README.md) | `f353f3868c1b2aa86de577480104b7ccb5901fbf` | SymPyスクリプトの実行手順 |
| 2周期例の有限検証 | [as4/research/z_periodicity/code/verify_periodicity.py](https://github.com/uedakazushi/as4/blob/main/research/z_periodicity/code/verify_periodicity.py) | `e4b1a5bb8306260aee67a5bc7d556d388d5e7b77` | 原実装を利用し重複回避 |

「コードがある」「そのコードをASZ4で再実行した」は異なる。非1周期性を判定する際の左射影像の非同型は必要条件に基づく数学的証明であり、計算PASSのみからは従わない。全14641の4周期性は左右AS条件とKoszul双対による。数値型が一般quadratic AS4に自動適用されるかは別問題。

## 外部一次資料

- A. I. Bondal, A. E. Polishchuk, *Homological properties of associative algebras: the method of helices*, Izv. RAN Ser. Mat. 57 (1993), 3–50; English translation Russian Academy of Sciences Izvestiya Mathematics 42 (1994), 219–260, DOI https://doi.org/10.1070/IM1994v042n02ABEH001536 。Koszul/Frobeniusとhelix、quantum determinant。使用箇所 Proposition 4.1, §5。原文適用の仮定は論文ごとに確認。
- I. Mori, A. Nyman, *A categorical characterization of quantum projective Z-spaces*, Applied Categorical Structures 33 (2025), DOI https://doi.org/10.1007/s10485-025-09806-2, arXiv:2307.15253。周期の定義、AS-regular Zとqgr、特にLemma 2.20（版による番号差は原稿照合）。
- M. Van den Bergh, *Non-commutative quadrics*, IMRN (2011), arXiv:0807.3753。3次元quadraticの1周期性とcubicの違いの比較。
- J. J. Zhang, *Twisted graded algebras and equivalences of graded categories*, Proc. LMS 72 (1996), 281–311。Zhang twistの一般概念（詳細な等価条件の適用は一次資料で再確認）。
- S. J. Sierra, *G-algebras, twistings, and equivalences of graded categories*, Algebras and Representation Theory **14** (2011), no. 2, 377–390, DOI https://doi.org/10.1007/s10468-009-9193-y, arXiv:math/0608791。graded categoryの同値・twisting system。書誌情報はエディンバラ大学の公開出版記録で照合。
- P. Belmansの[4d-AS-regularデータベース](https://4d-as.ncag.info/)（https://github.com/pbelmans/4d-as, 2026年閲覧）。graded側の既知族56件の検索・比較に有用。ただし一般ファイバー中心であり、完成した4次元Z代数の分類でも完全な全閉値網羅でもない。サイトの引用 DOI 10.5281/zenodo.21296619（サイトで表示、原リリースとの照合は未実施）。

上記の外部結果の新規性・仮定を推測しない。ASZ4の定理の証明で依存する結果は、実際の引用版の箇所を再確認する。

## 新旧同型の区別

- 1周期の`(C,\tau_1)`とgraded algebraの同型分類は同じ内容。
- `C`自身を分類するにはその周期同型を忘れ、Zhang twisting systemによる同一視を行う。
- 4周期表示 `(C,\tau_4)` と裸の `C` も同じ対象ではない。
- 2・4周期例から既約成分数や圏同値の非成立は即座に従わない。

## 再利用原則

AS4のgraded部門を止めたり、無断でAS4に書き込んだりしない。新規のZ側主張の証明はASZ4原稿に含める。AS4の既知定理の全文コピーは避け、必要な仮定と結果・固定参照で利用する。証明書が独立成果の再現に不可欠になった場合のみ、選んだものを監査・最小化して移す。

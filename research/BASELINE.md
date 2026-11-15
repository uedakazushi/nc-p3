# 依存・一次文献・証明位置 — 2026-11-15 / v11

## 現行プロジェクト

今回の出発点は `uedakazushi/nc-p3` のmain commit `98c0401f931123a796ee6df032b4a6430d668785`、tree `03ee44eec79aab9c77cbfbfc1ecbaa6bc77e7261`。主論文blob `d171321028573f757a0101d961e43762eaa841ac` を完全照合した上で、既存本文を保持して九変数族を統合した。遠隔保存状態は `validation.json` に記録する。

## as4から再利用する結果

`uedakazushi/as4/research/z_periodicity/papers/as4_z_periodicity_ja.tex` の本文blob SHA：
`cc475fae425a23a02fd8559e403ad2ec375302f3`（2026-11-15再確認）。

参照： https://github.com/uedakazushi/as4/blob/main/research/z_periodicity/papers/as4_z_periodicity_ja.tex

- 原稿前半の局所Frobenius・4周期性の議論：有限Koszul解消・左右AS双対性からのFrobenius性と4周期性。現行論文でもtrace構成を説明する。
- 2周期族の関係式・PBW表示・二次双対・対合行列を与える箇所：2周期PBW族、全λの基底、二次双対、完全対合。今回確認した本文に所在。
- 2周期例の非1周期性を扱う箇所：五根配置を使う非1周期性。λ=1で最小周期2。計算で恒等候補が失敗することだけとは異なる。
- 最小周期4の例と平坦族も同原稿を参照するが、今回その全証明・全計算を再実行したわけではない。

基礎原稿は私的・未査読、新規性未確定。現在の本研究の新規成果として再宣言しない。基準のas4全分類は未完成。

必要に応じて参照する既存資産のblob：
`research/z_periodicity/code/verify_periodicity.py`: `e4b1a5bb8306260aee67a5bc7d556d388d5e7b77`。
`research/z_periodicity/README.md`: `f353f3868c1b2aa86de577480104b7ccb5901fbf`。
これらは以前確認された参照識別子であり、今回の全コード実行を意味しない。巨大なリポジトリ全体を複製しない。

## 外部一次文献

1. Bondal–Polishchuk, *Homological properties of associative algebras: the method of helices*, Russian Acad. Sci. Izv. Math. 42 (1994), no.2, 219–260. DOI https://doi.org/10.1070/IM1994v042n02ABEH001536 。Koszul/Frobeniusとhelix・quantum determinantの基盤。既存の参照位置はProposition 4.1, §5。今回この論文全体の独立監査はしていない。
2. Mori–Nyman, *A categorical characterization of quantum projective Z-spaces*, Applied Categorical Structures 33 (2025), https://doi.org/10.1007/s10485-025-09806-2 ; https://arxiv.org/abs/2307.15253 。周期性・qgr・AS正則Z代数の枠組み。周期1とgradedの対応は既存版のLemma 2.20、版の番号差に注意。
3. Zhang, *Twisted graded algebras and equivalences of graded categories*, Proc. London Math. Soc. (3)72 (1996), no.2, 281–311. DOI https://doi.org/10.1112/plms/s3-72.2.281 。twisting systemの基盤。
4. Sierra, *G-algebras, twistings, and equivalences of graded categories*, Algebras and Representation Theory 14 (2011), no.2, 377–390. DOI https://doi.org/10.1007/s10468-009-9193-y ; https://arxiv.org/abs/math/0608791 。Zhang twistと付随Z代数同型の対応。大学の著者公開記録で書誌も確認。
5. Van den Bergh, *Non-commutative quadrics*, https://arxiv.org/abs/0807.3753v3 (2010); 出版版 IMRN 2011, no.17, 3983–4026（著者publication listで確認）。arXiv PDF pp.3,21、Proposition 3.3、Theorem 5.5.10、§5.6式(5.10)、Proposition 5.6.1、Proposition 6.6を比較に用いた。滑らかな楕円曲線の3パラメータ、正確な1周期条件、2周期性、graded4D ambientを混同しない。等差条件だけの記述は一般jに限定する。
6. Bhatoy–Ingalls–LaRoche–Nookala, *Some components of the moduli space of Koszul Artin-Schelter regular algebras of dimension four*, https://arxiv.org/abs/2511.08390v1 (2025-11-11)。§3.1の量子多項式、Table 2 (p.18) のCliffordのgraded HH²=9と量子多項式の6を確認。Table 2およびProposition 5.1 (p.20) はPDF画像でも確認した。Proposition 5.1を無条件の接空間判定として採用しない。尖点曲線の反例と滑らかな源を使う代替補題は現行論文・ノートに記す。

## 今回の結果の位置づけ

量子6→3は既知の量子族・twist理論を結ぶ具体的比較で、新規性を主張しない。4周期変形複体、Clifford点の全15変数の二次障害、相対標準化の局所適用は現行原稿に証明・計算を置いた。先行研究との最終的な優先権比較は未完。外部PDFや検索ログをbundleへ大量収録しない。


## 代数的七変数族の依存と位置づけ

- G. M. Bergman, *The diamond lemma for ring theory*, Advances in Mathematics 29 (1978), 178–218, Theorem 1.2。著者の訂正・更新ページ https://math.berkeley.edu/~gbergman/papers/updates/diamond.html 。停止するモニック書換えと全重なり解消から自由基底を得る一般論。参照確認は七変数族作成時に実施。
- S. B. Priddy, *Koszul resolutions*, Transactions of the AMS 152 (1970), 39–60, Theorem 5.3 (p.51)。一次資料 https://math.mit.edu/~hrm/palestine/priddy-koszul-resolutions.pdf 。七変数族の作成時に該当ページを画像でも確認。PBWからKoszulへの含意を使い、始点・終点を固定したbarフィルトレーションでZ版に適用する。
- 七変数族の全AS正則性は補助重みと双対の随伴次数Frobenius性で証明し、有限次数のHilbert一致だけからは結論しない。一般の非1周期性は二つの点射影の配置から直接証明する。
- 像の次元7は一つの関係空間へのPGL4×PGL4の階数30と全階数37による。全域モジュライstackの存在を仮定しない。具体的族の新規性は未確認。

## 九変数の楕円Clifford–Ore族

主論文 `thm:algebraic-nine`、証明ソース `sections/elliptic_ore_family.tex`、厳密検査 `computations/verify_elliptic_family.py`、証明書 `research/elliptic_family_certificate.json`。

- 三生成Clifford代数・Ore拡大・二次完全交叉のKoszul性は既存の一般手法。新しい一般定理として主張しない。本稿では特殊な三二次式の正則列条件、Clifford基底、導分の具体式、フィルター付きAS証明を記述する。
- Tate, *Homology of Noetherian rings and local rings*, Illinois J. Math. 1 (1957), 14–27, DOI https://doi.org/10.1215/ijm/1255378502 は二次完全交叉の線型Tate解消の背景。書誌は確認したが今回一次本文全体は取得していない。未確認の定理番号を補わず、本稿の特殊な場合の複体と完全性の説明を本文に置く。
- Christopher Phan, *The Yoneda algebra of a graded Ore extension*, https://arxiv.org/abs/1002.2318 は関連文献候補。著者公開書誌は確認したが本文の仮定・適用範囲の監査は未了で、本稿のAS証明の未確認依存にはしない。
- 楕円曲線上の次数12の因子、反不変6＋不変3の切断空間、有限射影自己同型の回避による非1周期性は本文で証明する。単一の数値点や七変数族の特殊化だけから一般点へ外挿しない。
- Bhatoyほかのgraded Cliffordモジュライとの比較は、通常gradedの変形と添字依存の交代変形を区別する。文献照合だけで九変数族が新しい既約成分だと結論しない。
- 九変数族の既存のOre型・Clifford型・非可換P3の族との完全な同定、新規性・優先権の最終確認は未了。新規性の不存在・存在を検索結果の有無だけで主張しない。

## Hesse三次曲線の不変量と指定Ore構成の導分分類

- V. Knibbeler, S. Lombardo, A. P. Veselov, *Automorphic Lie algebras and modular forms*, International Mathematics Research Notices 2023, no. 6, 5209–5262, DOI https://doi.org/10.1093/imrn/rnab376 。出版社HTML §7.3 式(42)にHesseの古典的$j$公式を確認。原稿 `prop:elliptic-j-obstruction` の具体式へ代入する。
- `computations/verify_ore_derivations.py`: 三生成Clifford部分の30×18制約行列、通常導分階数12・捩れ導分階数9を全パラメータで証明補助する小行列式を生成し、一般のindexed同型に関する楕円$j$必要条件の有理関数を検算する。
- 導分の完全性は固定Clifford部分と指定された符号Ore形式に限定。全四窓分類、楕円$j$一致による十分同型判定、全モジュライの成分最大性は未証明。

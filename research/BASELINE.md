# 依存・一次文献・証明位置 — 2026-11-15 / v9

## 現行プロジェクト

参照した `uedakazushi/nc-p3` のmainは commit `c387e8ceb8f1a3ca48e05d28d01c00aad69d4038`、tree `1e47763961d40eafbb0751f2bcfc431933a617bc`。主論文の確認blobは `1c8bc4047e3bda3813b88304ef3bc3e5539153e5`。これは統合元の固定参照であり、現行版の遠隔保存状態は `validation.json` に記録する。

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

- G. M. Bergman, *The diamond lemma for ring theory*, Advances in Mathematics 29 (1978), 178–218, Theorem 1.2。著者の訂正・更新ページ https://math.berkeley.edu/~gbergman/papers/updates/diamond.html を確認。停止するモニック書換えと全重なり解消から自由基底を得る一般論。
- S. B. Priddy, *Koszul resolutions*, Transactions of the AMS 152 (1970), 39–60, Theorem 5.3 (p.51)。一次資料 https://math.mit.edu/~hrm/palestine/priddy-koszul-resolutions.pdf の該当ページを画像でも確認。PBWからKoszulへの含意を使い、始点・終点を固定したbarフィルトレーションでZ版に適用する。
- 全ファイバーのAS正則性は、正の補助重みと双対の随伴次数Frobenius代数によって本文で証明する。有限次数のHilbert一致だけからは結論しない。
- 一般点の非1周期性は二つの点射影の配置から直接証明するため、既存as4の非1周期性定理や全1周期構造の分類には依存しない。
- 七変数族の一般同型類の次元は、一つの関係空間へのPGL4×PGL4作用の階数30とパラメータを含めた階数37から示す。全域のモジュライstackの存在を仮定しない。
- この具体的族の文献上の新規性、既存のOre型・非可換P3の族との完全比較は未確認。既存の一般手法を新しい一般理論として扱わない。

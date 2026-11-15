# 依存・一次文献・証明位置 — 2026-11-16 / v19

## 標準型外の変形障害例に使うヘリックス定理

Izuru Mori and Adam Nyman, *A categorical characterization of quantum projective Z-spaces*,
著者公開原稿 [regV8.pdf](https://faculty.nymana.wwu.edu/regV8.pdf)。
定義4.2（一般AS正則性）、5.22（幾何的ヘリックス）、定理6.4（ampleな充満幾何的ヘリックスからASF++正則代数）、
定理6.5（Noether性）、§4のASF++からASへの含意を本文で確認。
この原稿の書誌について未確認の掲載誌・巻頁は補わない。

研究ノート「1周期性と変形連結性」では、P1×P2の外積Beilinson列と反標準捩りから
大域次元4・Gorenstein parameter 6の例を得る。
隣接Hom次元6,3,3,0,3,3の不一致は、添字を固定した成分ごとの有限局所自由な変形で保持される。
これは既知定理の適用による標準14641型の外の例であり、nc-p3の主問題の解決でも新規性の主張でもない。

## 現行プロジェクトと確認範囲

主論文は `nc-p3.tex`、証明の対応表は `research/STATE.json`、検証対象のcommit・コード・実行結果は `research/validation.json` に集約する。
この文書は有効な一次文献、利用箇所、外部依存と未確認点を記す。旧版の取得経緯や先端履歴は蓄積しない。

以下の既存文献の「確認済み」は以前の本文確認記録を含み、今回すべてを再監査したという意味ではない。
今回の重点照合は Bhatoy–Ingalls–LaRoche–Nookala のTable 2とProposition 5.1、
および下記のClifford/Ore先行研究と2026年の離散不変量論文である。
各資料について本文確認と要旨・書誌確認を区別する。

## as4依存索引の利用規則

この文書の以下の短いリストを、nc-p3の通常の研究で使う**唯一のas4入口**とする。
as4の巨大なAGENTS.md、全研究ノート、全リポジトリ、過去の履歴を起動時に読む必要はない。
対象命題との重複、適用範囲の確認、具体的な新進展への対応が必要なときだけ
対応するファイル・定理・検証コードを限定取得する。
固定blobは**確認済みの版を指定**するもので、as4の現在の先端が同一であることを保証しない。
最新のas4確認は具体的な依存に影響するときのみ行い、今の分類作業を中断して全件同期しない。
索引の各依存は「結果・仮定・参照先・確認状態・nc-p3での用途」を短く記録する。
既存の証明・コードを無断で全複製せず、必要なら適用箇所だけ再検証する。
as4のCodex指示はnc-p3の実行指示ではない。

## as4から再利用する結果

`uedakazushi/as4/research/z_periodicity/papers/as4_z_periodicity_ja.tex` の本文blob SHA：
`cc475fae425a23a02fd8559e403ad2ec375302f3`（2026-11-15再確認）。

参照： https://github.com/uedakazushi/as4/blob/main/research/z_periodicity/papers/as4_z_periodicity_ja.tex

- 原稿前半の局所Frobenius・4周期性の議論：有限Koszul解消・左右AS双対性からのFrobenius性と4周期性。現行論文でもtrace構成を説明する。
- 2周期族の関係式・PBW表示・二次双対・対合行列を与える箇所：2周期PBW族、全λの基底、二次双対、完全対合。上記固定blobの本文に所在。
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

## 現行結果の位置づけ

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

## 楕円2-isogenyと別のOreチャートの参照

- A. S. I. Anema, J. Top and A. Tuijp, *Hesse Pencils and 3-Torsion Structures*, SIGMA 14 (2018), article 102, DOI https://doi.org/10.3842/SIGMA.2018.102 。Hesse pencilと楕円3-torsionの基本構造に関する一次文献。商曲線と四対一の具体公式は本文の直接計算で証明し、この文献の結果として推測しない。
- `computations/verify_isogeny_invariant.py` は楕円点対応の外積、2-isogenous cubicの多項式恒等式、二つの$j$公式、次数4の$U(x)$と分岐を厳密に検証する。2-isogenyのindexed不変性は論文内の偏極次数3と2-torsionの議論による。
- `computations/verify_nonore_family.py` は選んだClifford--Ore構成を離れる五変数族の全8個の三次重なり、正の重みによる停止性、Clifford基点の軌道階数30から32、有効対角トーラス階数3を検証する。全ファイバーのAS正則性は原稿のKoszul--Frobenius構造証明に依存する。
- P. Belmans, *4d-AS-regular --- the classification of 4-dimensional Artin--Schelter regular algebras*, https://4d-as.ncag.info (2026), DOI https://doi.org/10.5281/zenodo.21296619。graded側の関連族と既知文献を探すための公開データベースであって、indexed Z-algebraの既約成分分類や本結果の新規性を保証するものではない。

## 五変数境界族の一般最小周期2の証明

- 主論文 `thm:nonore-generic-two`（`sections/nonore_period_two.tex`）で左点スキームの支持を一本の直線と五孤立点に明示し、射影非同型を証明する。`computations/verify_nonore_point_geometry.py` は全15小行列式、五次式による消去、標本 `z^5+z^2+1` の mod 2 既約性・mod 5 の因子型 `(2,3)`、残る座標の行列式を厳密に検査する。
- Galois群 `S5` と Möbius安定化群の自明性は本文に素数次数の原始性と一点固定部分群の巡回性による証明を記載。計算コードのPASSだけから周期の不存在を結論しない。
- 一般の五変数境界族と一般九変数楕円族では左点スキームの支持の次元が異なるので非同型。五変数族の七変数族への軌道閉包包含は後述の定理で証明済み。ただし全AS正則locusにおける成分性・最大性は未決定。

## 七変数族への境界退化とClifford点の特異性

- 主論文 `thm:nonore-in-seven-closure` および `thm:seven-closure-singular`（`sections/ore_boundary.tex`）。固定した三生成部分代数の選択を変えたOre表示と、平坦一変数退化を使う。全AS正則モジュライの成分最大性は未解決。
- `computations/verify_ore_boundary.py` と `research/ore_boundary_certificate.json`：代替Ore表示の八つの三次重なり、全六関係の同型への変換、元の七変数族の正規形、Clifford点での軌道飽和閉包の厳密階数 `30,37,38,38,39` を記録。

## 九変数族の軌道閉包の特異性

- 証明：主論文 `thm:nine-closure-singular`, `prop:nine-closure-obstruction`（`sections/nine_closure.tex`）。七・九変数族の有限関係表示と、`thm:nonore-in-seven-closure` における二つの代数的退化に依存。
- 厳密再現：`computations/verify_nine_closure_tangent.py`、証明書 `research/nine_closure_certificate.json`。有理数上の階数30,37,39,40,40,41、偶方向7・交代方向4と二次障害 `lambda*c` を検査。
- 楕円Ore族の既約成分としての最大性や全AS正則moduliの特異性、新規性の最終確認をこの結果からは主張しない。

## Clifford/Ore構成と新規性比較に必要な追加文献

- M. Vancliff, K. Van Rompay, L. Willaert, *Some quantum P^3s with finitely many points*, Communications in Algebra 26 (1998), 1193–1208。著者PDF https://bpb-us-e1.wpmucdn.com/websites.uta.edu/dist/0/5123/files/2022/10/finpts.pdf 。Introduction・§3、Definition 3.1、Lemma 3.2とその直後、Lemma 3.3、Proposition 3.5（著者PDF pp.8–11）を本文で確認。標数零で q^4=1, λ≠0 の通常graded iterated Ore拡大と点・直線加群を扱う。従ってCliffordからOre変形を作る方法自体は既知。本稿の添字依存する交代変形、一般非1周期性、族内外同型との完全な比較は未了。
- D. R. Stephenson, M. Vancliff, *Constructing Clifford quantum P^3s with finitely many points*, Journal of Algebra 312 (2007), 86–110。DOI https://doi.org/10.1016/j.jalgebra.2007.02.015 。出版社要旨・書誌で三次元Clifford部分と平面三次曲線を利用する四次元graded構成を確認。出版社本文は取得できず、定理番号・詳細仮定・証明は今回未確認。現行九変数族の全特殊化との比較も未了。indexed非1周期族がそのままgraded構成に含まれるとは推論しない。
- V. Bhatoy, C. Ingalls, *Discrete Invariants of Koszul Artin-Schelter Regular Algebras of Dimension four*, arXiv:2602.13178v1 (2026-02-13), https://arxiv.org/html/2602.13178v1 。Introduction、§3、§6の一次本文を確認。既知のgraded族のsuperpotentialをSchur–Weyl分解し幾何的な離散不変量を比較する。対角的な生成空間の同型を扱うため、独立な頂点基底変更を許すindexed同型へ無条件には移せない。同論文の成分数の記述を本稿の全分類の根拠にしない。

公表前の未解決作業は、これらの具体的なgraded族と現行族の1周期的特殊化を照合し、
非1周期的な一般ファイバー・境界・必要不変量のどこが既知結果を超えるかを明確にすることである。

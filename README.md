# nc-p3 — 非可換射影三空間の完全分類

**共同研究：** 本リポジトリは大川新之介・植田一石の共同研究プロジェクトです。この表示は個々の論文の著者を確定するものではありません。

**2026-11-28 / v555。標準14641型の完全分類は未完成。正本は GitHub main。**

## 固定する最上位目標

**数学的正確性を保ち、標準14641型 quadratic Koszul AS正則 indexed Z-algebra の添字保持同型による完全分類までの最短コースを進む。この最上位目標は変更しない。**
完全分類前の途中公表は考えず、部分結果の投稿・論文分割・公表用の仕上げを目標や工程に入れない。
Codexは公表しやすさ、分類の巨大さ、局所問題の難しさを理由に、目標を部分分類へ変更したり対象範囲を狭めたりしない。
手法・中間課題・優先順位は、完全分類への寄与に基づいて変更する。

完成には、全対象の被覆、全パラメータのAS正則性、特殊値を含む添字保持同型の必要十分条件、周期層、既約成分・次元・特殊化の決定が必要である。
分類の枝数・族数・原稿量に上限を設けず、有限検査を通った候補の列挙や、未決定のAS条件を残した一般方程式への還元を完成扱いしない。
既存の有効な定理・証明・族・反例・計算資産は、完全分類の道具と検証入力として保持する。

四辺の点対応が有限で、ある点の接空間が三次元なら、その点は一意であり、この枝の添字保持同型類は次の四種類で尽きる。各行に一般最小周期2の型と4の型がある。[全被覆の証明](sections/maximal_tangent_classification.tex)・[三生成基底の復元](sections/maximal_tangent_ore.tex)

| 三生成基底の点曲線 | 各型の同型類の次元 | 正規表示・特殊値での同型と周期 |
| --- | ---: | --- |
| 滑らかな平面三次曲線 | 9 | [位数2・4の平行移動と次数12の因子](sections/elliptic_tangent_windows.tex) |
| 節点三次曲線 | 8 | [節点評価と平滑化](sections/singular_tangent_smoothing.tex) |
| 横断する直線と二次曲線 | 7 | [二成分の評価と倍率・反転](sections/reduced_ternary_point_graphs.tex) |
| 非共点三直線 | 6 | [九係数の量子基底表示](sections/quantum_tangent_windows.tex) |

後の三種類は対応する楕円族の特殊化である。三次元接点を持たない全ての有限四辺でも、局所点環と三つの四次元空間から添字保持同型と周期を判定できる。埋込み次元2の局所環では、最初の二元関係から高々六つの直線を選び、局所環の定義イデアルと三空間を同時に正規化する。[同型の完全不変量](sections/planar_point_flags.tex)

Walsh基底を持つ隣接階数6の四窓は、符号付き平方行列 $N$ の階数が2のとき、かつそのときに限り標準Koszul左右AS正則である。四つの正規二次式の二乗上で各始点が自由階数256になる。指定なし商の正規化は六平面に分岐する二重被覆の三文字置換商であり、特殊値の同型と最小周期も三切断の特異値から有限候補を復元して比較できる。[全点の正則性](sections/walsh_normal_quadrics.tex)・[正規化と特殊点の判定](sections/walsh_spectral_fibers.tex)

この枝の階数12は一つの量子窓、階数14は二点除去直線、階数15は境界四点を含む有理曲面の開集合である。階数14の周期1は一点、階数15の周期1は二点である。[テンソルからの認識](sections/walsh_symmetric_windows.tex)・[階数14の同型曲線](sections/walsh_rank_fourteen_orbits.tex)・[階数15の全点分類](sections/rank_fifteen_all_ordinary.tex)

Walsh三次元族を含む十五次元成分と、そこから出る周期2の成分では、次のことが分かっている。

- この閉包の有限四辺では全正則点が中心有限な合成整域になり、正則性を例外値も含めて16次までで判定できる。実現した窓の添字保持同型は偏極中心とKummer印で決まる。[全正則点の証明](sections/main_regular_specialization.tex)
- 十五次元成分では、指定の窓に法53で還元される持上げは全て正則で、最小周期4の正則な代数的閉点が稠密にある。[良い還元の判定](sections/unrestricted_good_reduction.tex)
- 一般Walsh点から出る周期2の局所枝の次元は9、5、5、5である。後の三枝は大域では一つの成分に属し、二倍率と二つの二次式を満たす平面から一般窓を復元できる。[四つの周期層](sections/walsh_full_component.tex)・[成分の正規化](sections/half_period_components.tex)・[五次元成分の復元](sections/half_period_two_quadrics.tex)
- 九次元成分の一般窓は滑らかな二次曲面、四直線枠と三自由重みで決まる。三重みを変えても中心商は一定で、枠のGram行列と逆行列の対角が全て非零なら全ての許容重みで正則である。[枠の復元](sections/half_period_frame.tex)・[共通中心商](sections/half_period_weight_fibers.tex)・[正則性の証明](sections/half_period_semisimple.tex)
- 十五次元成分が量子多項式族と交わるのは反対階数12の一同型類だけである。この点へ退化する全窓は正則で同成分に属し、十八座標の台とトーラス・偶置換で分類できる。対角量子Ore枝の所属も量子背景で決まる。[共通部分の証明](sections/quantum_component_intersection.tex)・[Ore枝の判定](sections/quantum_ore_component_incidence.tex)・[全台の同型・周期・退化](sections/quantum_attraction_fiber.tex)

量子点への退化は、四つの向き付き三角積でも判定できる。四値が相異なればその量子点だけが残る。一対だけが非単位の重複値なら、三係数の台と一つの比で全同型類が決まる。三面が一致する場合は九係数表示のうち平方係数の積が消える三平面で尽くされ、周期は全て1である。[相異なる場合と一対重複](sections/quantum_face_resonance.tex)・[三面一致の場合と九係数表示の正則性](sections/quantum_triple_face.tex)

二つの非単位値が二面ずつに現れる場合は、全窓が周期1のOre表示を持つ。二値の積が負なら四係数が自由、正なら三次重なりが混合係数と反対側平方係数の同時非零を禁じる。[二対重複の特殊値・同型・退化](sections/quantum_double_pair.tex)

可換点を含む任意の量子点への退化は、文字の置換を除いて三種類の余指標で尽くされる。各表で重み19までの線型復元を行うと、全次数で正則な窓を生成し、障害のある入力を棄却できる。特殊値の同型・周期・同じファイバー内の退化も有限判定できる。[三表と復元手順](sections/quantum_universal_chambers.tex)・[厳密検算](computations/verify_quantum_weight_recursion.py)

正則正規な一次直線列を持つ窓の最小周期は1である。滑らかな楕円三次商ならその直線列は窓から一意に戻り、三生成商・線型補正・定数補正から全ての同型を判定できる。固定した滑らかなSklyanin商では、分裂点と重み付き射影空間 $\mathbb P(1,1,1,2,2,2)$ の有限群商が全特殊値を与える。節点・尖点の整三次商でも直線列は一意である。一般の四窓からは四直線の階数検査で全ての印を取り出し、非被約商を含めて正則性と印なし同型を判定できる。[中心斉次化の分類](sections/regular_normal_linear.tex)・[四窓からの有限認識](sections/normal_line_recognition.tex)

重み付き反復Ore環の四生成元が窓に入り、窓がその上で左有限なら、四次数の和までの経路で左右AS正則性が決まる。パラメータの正規性や中心性は不要で、無限位数の量子例では非正規なパラメータを使って次数5で認証できる。[片側有限パラメータの判定と証明](sections/onesided_parameter_cutoff.tex)

一辺が滑らかな二次曲面のグラフで全辺の射影ファイバーが零次元なら、全標準窓の同型類はアフィン直線から一点を除いたものとなる。全て周期1で、量子極限も正規一次直線列も持たない。これらの各窓に退化する全同型類は一般に五つで、そのうち二つの周期が2である。例外値では三つに減り、全て周期1となる。[四辺の復元と二次曲面枝](sections/smooth_quadric_graphs.tex)・[各閉軌道への全退化](sections/quadric_attraction_fibers.tex)

余分な直線や非被約構造を許しても、滑らかな二次曲面商を持つ全窓は二行二列の行列対で分類できる。各辺の射影ファイバーが零次元なら、真の曲面因子の総次数は高々2である。[全二次曲面商と次数上界](sections/smooth_quadric_quotients.tex)

整二次錐商を持つ窓も、双対数上の二行二列行列で全特殊値まで分類できる。平面成分があれば、射影ファイバーが一点という条件の下で正規一次直線列が復元される。これにより、この条件下の二次元成分は滑らかな二次曲面・二次錐・平面商の三通りで覆われる。[二次錐の正規形と曲面被覆の証明](sections/quadric_cone_quotients.tex)

平面商については、三生成商の移動の全Jordan型で特殊値を含む正規形が決まった。可換多項式窓以外で印の候補は高々二平面で、三つの二次式の共通因子から第二の印も回収できる。[平面を忘れた全同型判定](sections/plane_quotient_normal_forms.tex)

一次元の点グラフでは、埋込み多重度を含む曲線次数は高々7で、各成分の次数は全位置で一定になる。楕円四次成分の偏極の二階差は12で消える。[曲線の次数と偏極](sections/point_graph_curves.tex)

滑らかな楕円四次成分を持つ点グラフでは、切断環への二本の二次核の積を四次で検査すれば、全次数のAS正則性が決まる。二窓と楕円曲線の有限群から同型と周期も判定できる。[楕円四次点成分の分類](sections/elliptic_point_quartic.tex)

次数7に達する点グラフは純粋な算術種数5の曲線となる。滑らかな場合の候補は位数4の差を持つ二点からの標準射影で決まり、三つの可換文字による自由性で全次数の正則性を判定できる。[次数七の点グラフ](sections/septic_point_curves.tex)

滑らかな七次点曲線の枝は八係数の三つの標準二次式で尽くされ、許容条件は完全交叉の滑らかさだけである。同型類は5次元、一般最小周期は2で、一周期層は混合係数の符号変更と三次元トーラス・交換で判定できる。[七次枝の全正規形](sections/septic_bielliptic_forms.tex)

四次完全交叉の印が有限軌道を持てば、特異・可約・非被約な場合も移動する二次核の四次積で全次数AS性を判定できる。被約な印はすべて列挙できる。[種数一の四次部分曲線](sections/genus_one_point_marks.tex)

正次元射影ファイバーを許しても、四移動の決定式対応の反復可能な支配的枝は双有理で一般重複度1に限られる。一つの有向道が通る相異なる巡回成分は高々四つである。[反復する点対応の成長制限](sections/recurrent_point_correspondences.tex)

有限点グラフでは六関係に階数1の元が存在せず、反対階数6・7の有限Frobenius候補は全て標準AS正則となる。より一般に三支持が全て三次元の縮約があれば、三生成部分代数への有限判定が使える。[三次縮約による全次数判定](sections/ternary_point_subalgebras.tex)

有限な四辺を持つ全Frobenius候補で、標準四移動が長さ20の点スキーム上の恒等射になることを、反対階数を問わず証明した。三行三列の行列式は反対側の対応点の局所socleを分離する。[有限枝の標準四移動](sections/finite_cubic_separators.tex)

有限点枝の被約な印では、全四テンソルが四組の跡零三行行列と一つの交換子方程式で決まる。点評価の核は三文字・八本の三次交換子関係だけで表示できる。[被約な印の行列表現](sections/marked_point_commutators.tex)

曲線上の反対二点縮約の階数低下は、曲線の平方上で両側次数2の因子になる。楕円四次の偏極二階差は6で消え、滑らかな七次枝の標準四移動は二重被覆の対合そのものである。[曲線上の階数低下因子](sections/opposite_curve_divisor.tex)

非超楕円・種数3の滑らかな六次点グラフは存在しない。三割線とJacobian作用からFermat平面四次曲線に限定し、残る30本の偏極列を三次の厳密階数で全て排除した。[非超楕円六次の排除](sections/nonhyperelliptic_sextic_exclusion.tex)

滑らかな整六次曲線と孤立残余から成る点グラフは全て排除される。超楕円の場合は二つの鉛筆を用いて三次商の次元を18以下と評価する。[六次曲線の排除](sections/smooth_sextic_exclusion.tex)

滑らかな整五次点グラフと孤立残余の場合、種数2は存在せず、残るのは楕円五次曲線そのものとなる。[五次の偏極と種数2の排除](sections/genus_two_quintic_ruling.tex)

楕円五次候補では楕円曲線の j 不変量は1728、標準四移動は平行移動となる。各偏極で射影中心は高々十点、四位置で同じ位数5の群を保ち、窓は有限個の一パラメータ族に限定される。[三割線による限定と乗法核](sections/elliptic_quintic_trisecants.tex)

楕円五次候補の三次重なりを全パラメータで検算し、二つの投影それぞれで二捩れ点を交互に移す関係列だけが残った。[全幾何点の検算と標数零への持上げ](sections/elliptic_quintic_cubic_overlaps.tex)

滑らかな整五次曲線と孤立残余を持つ点ファイバーの枝は、添字保持同型までちょうど二類である。どちらも孤立部のない楕円五次で、裸周期1、全次数でKoszul左右AS正則となる。[二類の実現と全次数の証明](sections/elliptic_quintic_regularity.tex)

滑らかな整四次曲線と孤立残余の場合、有理四次は二次syzygy束とHankel行列式により排除され、楕円四次と長さ4の孤立部に限られる。[有理四次の排除](sections/rational_quartic_exclusion.tex)

滑らかな平面三次と孤立残余の枝は、三捩れ以外の平行移動を持つ三生成Sklyanin商の全中心PBW持上げに一致する。裸周期は1で、孤立部は明示三次式の長さ8の臨界スキームである。[平面三次枝の被覆](sections/plane_cubic_point_graphs.tex)


滑らかな点曲線と孤立残余の枝は、平面三次・楕円四次・楕円五次・種数5の七次の四行で尽きる。四移動の固定イデアルを余接層から計算すると、滑らかな有理成分はどの次数でも生じない。[固定イデアルと全四行の証明](sections/point_graph_fixed_differentials.tex)

次数3以下の被約な点曲線と孤立残余は、滑らかな三次のほか三角形だけである。三角形の全正則窓は六係数の中心持上げで裸周期1、全て滑らかな三次窓から平坦特殊化する。[係数・同型と境界](sections/reduced_cubic_point_curves.tex)

被約四次点曲線と孤立残余の特異な場合は、二つの二次曲線の二角形か空間四角形に限られる。各模型で全次数のAS性、同型、周期を係数から有限判定できる。[被覆と多角形係数の判定](sections/reduced_quartic_point_curves.tex)

整な被約五次点曲線と孤立残余では特異曲線も全て除かれ、正則窓は滑らかな楕円五次の二類に限られる。[三次鉛筆・固定長と錐頂点からの排除](sections/singular_integral_quintic_exclusion.tex)

整な七次点グラフの全窓は、特異点を許しても八係数の三二次式が整な完全交叉を定める場合に一致する。全て滑らかな五次元層の平坦特殊化で、周期と同型は同じ倍率・交換判定となる。[被覆と尖点の実例](sections/integral_septic_boundary.tex)

整六次点曲線では三次元接点と正規化種数2も排除され、孤立部のない特異平面四次の正規化種数0,1だけが残る。非恒等な四移動の位数は2または4である。[三割線と固定点による排除](sections/singular_sextic_normalization.tex)

孤立残余と交わらない整な被約点曲線の全窓は次数3・4・5・7の四行で尽きる。六次の特異曲線も全て排除した。特異な整曲線は七次の許容係数だけに現れる。[全被覆と六次排除の証明](sections/integral_sextic_exclusion.tex)

埋込み点を許しても六次点グラフは純な種数3の曲線となる。整な純部分を持つ五次以上では有限構造が消える。[純性と固定長の証明](sections/pure_curve_finite_defect.tex)

被約な種数一の五次点曲線は滑らかな楕円五次の二類で尽きる。五角形は各節点での反対縮約の階数で排除される。[種数一の成分と五角形の証明](sections/genus_one_polygon_exclusion.tex)

直線・二次曲線成分があれば全次数AS性と添字保持同型を有限に判定できる。これで埋込み点・非被約構造を含む五次点グラフは全窓が被覆された。[直線成分](sections/line_component_as_test.tex)・[二次曲線成分と五次の分類表](sections/conic_component_as_test.tex)

点グラフ内に平面三次があれば、埋込み点や重複成分を許しても中心一次PBW表示と最小周期1が強制される。[平面三次の印からの復元](sections/plane_cubic_mark_normality.tex)

有理正規三次の二重曲線を全て排除し、七次点グラフの全窓を四行の有限アトラスで被覆した。六次で既出の印を持たない場合に残るのは二本の有理正規三次だけである。[二重曲線](sections/double_twisted_cubic_exclusion.tex)・[七次の被覆](sections/reducible_septic_coverage.tex)

六次の最後の候補であった二本の有理正規三次の和を、非被約な交点も含めて排除した。点ファイバーで一次元の点グラフは次数4以上の全窓が被覆された。[六次の証明と全次数の残余](sections/binary_sextic_exclusion.tex)

有理正規三次上に基点を持つ階数二関係では、点加群による短完全列から直線加群の全次数の次元と関係を決めた。この窓の正則性は全次数のHilbert関数に帰着する。[証明と割線の行列式](sections/twisted_cubic_line_modules.tex)

残余三次の二項行列式が整なら、節点・尖点を含めて種数一曲線の加法から高々四つの鉛筆を結び、全次数正則性を証明した。射影の縦横の成分でも三文字部分代数から正則性が従う。[証明と残る二種類の因子](sections/elliptic_binary_cubic_pencils.tex)

二重の一次対一次グラフは全て正則である。三次残余の全域で割線の対応が平面三次どうしの同型となり、正則性を証明する両方向の割線列の存在を有限回で判定できる。[二重化・二面体作用・閉集合の判定](sections/doubled_binary_pencils.tex)

残っていた二本のメビウスグラフを、二節点・接点の両方で解決した。これで点ファイバーを持つ一次元点グラフの全枝に、全次数正則性と添字保持同型の有限判定が得られた。[混合型・全被覆・三次残余の同型と周期](sections/split_binary_pencil_charts.tex)

有限点グラフの反対階数8では、長さ1の局所成分が一つあれば、三支持による既知の場合以外は必ず両対称表示になる。従って既知の有限判定で全次数の正則性が決まる。[交換子の退化を含む証明](sections/rank_eight_finite_symmetry.tex)

長さ1の局所成分を持つ反対階数9の有限窓は、少なくとも一方の反対位置が対称になる。その対称核の二次形式が階数2以下なら全次数正則である。[四次元行列空間と階数9の制限](sections/rank_nine_finite_symmetry.tex)

片側対称な有限窓で対称核が滑らかな二次形式なら、二つの二項方向の階数から全次数正則性を判定できる。[二項切断と非全射の場合の増大](sections/smooth_kernel_binary_transfer.tex)

反対階数9で長さ1の局所点を含む有限窓は、中支持3の点があるとき、かつそのときに限り正則である。[二次錐と必要十分性](sections/cone_kernel_binary_transfer.tex)

交代・対称混合型の印付き有限窓は長さ4の可換Frobenius環と三変数二次形式で分類され、全て階数10・最小周期2の正則窓である。非有限の場合も局所socleの二次式で判定できる。[混合型の全パラメータ](sections/mixed_point_frobenius_algebras.tex)

混合型の全正則窓は最小周期2で、長さ4の環の正規形と残余群、点グラフの有限性が決まる。八次元の閉包は九次元半周期成分の因子となる。[正規形と特殊化](sections/mixed_point_orbits.tex)

全標準正則窓で射影平面以上の点ファイバーを排除した。残る直線ファイバーは支持平面ごとに高々一点で起こる。[端支持と経路成長](sections/endpoint_plane_growth.tex)

直線ファイバーの中心が曲面をなす場合を一般に排除した。曲線をなす中心の直線は、次の位置で既約曲面を掃く。[分岐中心の次元](sections/branching_curve_bound.tex)

全標準正則窓で直線ファイバーの中心の台は有限本の直線と有限個の点に限られ、各直線が掃く曲面は平面か二次曲面である。[二項切断と曲線の排除](sections/branching_lines_only.tex)

動く直線ファイバーは平面の直線束に限られ、反対辺に純関係を強制する。純関係のない窓では分岐中心は有限個である。[平面束への帰着](sections/plane_branch_pencils.tex)

純関係のない窓では、孤立した直線ファイバーを許しても全ての曲面射影を既知の分類で覆う。[二次曲面間の双有理対応](sections/birational_quadric_chains.tex)

純関係を許しても、直線ファイバーの中心が有限なら曲面射影を持つ全窓を既存の分類に帰着できる。[孤立ファイバーを許す曲面枝](sections/isolated_surface_fibers.tex)

曲面射影を持つ全ての標準正則窓を被覆した。可動直線束の存在は二重対角を含めて指数増大を強制する。[可動束の排除](sections/exclude_moving_plane_pencils.tex)

孤立分岐の直後の動く曲線は次数三以下の有理正規曲線である。両側分岐には四文字の自己経路が一本あり、前後の例外ファイバーを二項二次式で検出できる。[孤立分岐の二側切断](sections/isolated_branch_blocks.tex)

両側に直線ファイバーを持つ点も排除した。存在すれば二つの中間境界の非自明な一次変形を固定端点の間で反復できる。[両側分岐の排除](sections/exclude_two_sided_centers.tex)

片側分岐中心は有限な有理曲線列で逆側の分岐中心と一対一に対応する。次数の一次モーメントから、分岐を通る周期点列が必ず存在する。[曲線列の次数公式](sections/finite_branch_corridors.tex)

曲面のない片側分岐があれば、長さ四十以下の周期点加群が存在し、右・左とも無限生成の表示イデアルが必要である。[周期と表示イデアル](sections/periodic_branch_obstruction.tex)

反対階数七以下は曲面・曲線・有限射影を含めて全正則窓を被覆した。純関係のない残余分岐には反対階数九以上が必要である。[低階数の分岐排除](sections/low_rank_branch_exclusion.tex)

純関係のない二直線分岐には反対階数十二以上が必要で、両側の例外二点は反対三次曲線との割線・接線交叉である。[二直線の切断](sections/two_line_branch_blocks.tex)

単純な有限点で一組の交換子三平面が一致すれば、反対階数の上限を仮定せず片側対称となる。[共通三平面の対称性](sections/marked_common_triple_symmetry.tex)

交換子が二つの共通三平面では可逆な束の元が存在し、Sylvester方程式から三平面を復元できる。[二交換子束の復元](sections/marked_two_commutators.tex)

両側の印付き線形包が四次元なら交換子階数は4か5で、両側で二つの三列の差を同時に階数1にできる。[四平面対](sections/paired_marked_four_planes.tex)

被約点を仮定せず、反対階数8の印付き行列の核と接方向から反対位置の対称性を判定する十分条件を得た。[非被約の印からの対称性](sections/finite_marked_projection_symmetry.tex)

階数1切断の長さが2以下なら、曲線状の印では例外割線があっても反対対称性が従う。[割線と支持長](sections/finite_marked_secant_pencil.tex)

階数1切断の長さ2以下では接方向の例外も解決した。階数1点が三つの被約点の場合には曲線状の印で対称となる。[接方向の割線](sections/finite_marked_tangent_secants.tex)

印付き切断の階数1軌跡が正次元なら被約な一本のSegre直線で、非被約の印でも反対対称性が従う。[正次元の階数1切断](sections/finite_marked_ruling_pencil.tex)

長さ3のGorenstein階数1切断では三重方向を含め曲線状の印から対称性が従い、長さ4の切断は有限環の積跡で正規化できる。[Artin割線](sections/finite_marked_artin_secants.tex)

曲線状の印では非Gorensteinの三重点も二平面の長さ6の完全交叉によって対称になり、長さ3の全型を覆う。[三重点と余因子](sections/finite_marked_planar_triple.tex)

全ての有限窓で反対座標積は10次元か16次元であり、双一次式が一つでも消えればその側が対称になる。[反対座標積の次元](sections/finite_opposite_product_gap.tex)

被約四頂点の切断では反対座標の六つの非対角積が一致し、対称でなければ六つの差が互いに交わらない長さ1、2の部分に残る。[四頂点の積](sections/finite_marked_tetrahedral_products.tex)

純関係も二直線列もない分岐では順逆中心が相対する二位置に一つずつあり、列長対は三通りだけである。[分岐位置と列長](sections/branch_centers_opposite_phases.tex)

純関係のない曲線分岐には反対階数11以上が必要で、階数11の交点の束の分裂型は (2,1,1) である。[分岐の階数下界](sections/branch_rank_eleven_floor.tex)

純関係を許しても反対階数8の曲線分岐には中央支持2の中心が必ずある。階数9にそれがなければ全中心の中央支持は3となる。[中央支持の次数障害](sections/branch_middle_defect_bound.tex)

反対階数8では二直線列と中央支持3の中心を排除した。各位置の順中心・逆中心は各々高々一つである。[第二曲線への制限の証明](sections/branch_rank_one_residual.tex)

純関係を許す反対階数8でも、曲面射影のない正則窓の全ての辺のファイバーは空か一点である。[分岐の排除と二次ジェットの証明](sections/rank_eight_branch_exclusion.tex)

四つの被約頂点を持つ共通切断では、印の接次元が1以下なら行・列に当たる反対二位置が対称になる。他の局所成分は非被約でもよい。[外積の支持と三次分離](sections/finite_tetrahedral_kernel_symmetry.tex)

長さ4の非被約切断にも行・列側の対称性を拡張した。反対階数8の有限窓に接次元1以下の印があれば三文字認証か、少なくとも片側の対称表示に入り、全次数の正則性が有限に決まる。[Artin環上の証明](sections/finite_artin_kernel_symmetry.tex)

反対階数8の有限枝では、極大イデアルの二乗が単項の平面局所点も全次数で判定できる。未判定候補の閉点は高々三つで、局所長は全て6以上である。[平面局所点の証明](sections/finite_stretched_planar_marks.tex)

反対階数8の有限枝で長さ6の局所因子も全次数判定に入った。未判定窓の台は一点または二点である。全有限枝の直線切断長は3以下である。[節点・尖点と直線切断](sections/finite_nodal_cuspidal_marks.tex)

分離三次式の重平面因子から切断の正規形を決定した。台が一点なら長さ4の全切断で接次元によらず行・列側の対称性が従う。[分離式と切断の証明](sections/finite_split_separator_slices.tex)

純関係があっても二直線列は反対階数11以上でしか起きない。階数11では両端中央支持3・端点縮約階数1に限る。[直線への全制限](sections/two_line_full_middle_restriction.tex)

全有限窓で平面切断長は7以下である。反対階数8の二点間の接平面入射は反対位置で一致し、同じ位置で双方向にはなれない。[平面切断と二点の入射](sections/finite_plane_sections_and_tangent_incidence.tex)

局所定義式が四次以上の未判定有限階数8窓には、四座標が線型同定される長さ8の三二次完全交叉が必ずある。そのイデアルは接平面の式の二乗の零化イデアルである。[共通八点切断](sections/finite_high_order_common_octad.tex)

共通の三二次完全交叉の網から四窓の両対称部分と四形式への分解が従う。有限な反対階数8の未判定窓では、各平面局所点の定義式に二次または三次で始まるものがある。[共通網と対称性](sections/finite_common_net_symmetry.tex)

**四辺が有限で反対階数8の全候補について、全次数の正則性を有限判定できる。** 一点または二点に集中した非被約な候補も含む。平面局所点では重複度2・3の曲線を一回ブローアップし、共通切断の誤差を完全跡で消す。[有限階数8の全窓](sections/finite_rank_eight_all_marks.tex)

**反対階数8以下では全窓の正則性と添字保持同型を有限判定できる。** 四辺が有限で反対階数9なら正則性は中央支持3の点の存在と同値で、階数9の被覆の残余は純関係を持つ曲線分岐である。二点支持では双方の長さ4切断が強いる接平面入射を使い、一点の高重複度では共通八点切断から四形式の階数下界を得る。[有限階数9の全窓](sections/finite_rank_nine_all_marks.tex)

中央支持2の順分岐中心は全て、直前位相の中央支持2の逆中心と一意に対をなす。分岐列で一般中央支持3の成分は対の端の直線だけである。階数9では連続する非直線成分のうち非平面的なものは高々一本となる。[分岐中心の対応](sections/branch_binary_center_pairing.tex)

**反対階数9以下の全窓で、全次数の正則性と添字保持同型が有限判定になった。** 階数9の曲線分岐では三次元中心が平面二次曲線と純関係を強制し、逆向きとの両立が不可能になる。[階数9の全分岐の排除](sections/rank_nine_branch_exclusion.tex)

四辺が有限で全中央支持4なら、局所二式がともに四次以上となる印一つから、反対階数によらず窓全体の両対称性が従う。正則なら階数は8以下である。[高重複度局所点からの全域への延長](sections/finite_high_order_half_turn.tex)

被約二十点で二つの印付き四空間が異なり五空間を張るなら、階数一切断には必ず曲線がある。六点の割線のうち十三本だけでも、反対座標の平方の関係から二超平面に矛盾が生じる。[五空間の割線判定](sections/finite_marked_five_space_secants.tex)

被約二十点・全中央支持4では、五空間の母線を除外でき、平面二次曲線からも反対対称性が従う。このため反対階数10は少なくとも片側対称、階数11は生じない。[五空間の階数一曲線](sections/finite_five_space_rank_one_curves.tex)

階数10で中央支持4の順中心の第二曲線がねじれ三次なら、その上に別の順中心があり、第三次数は4以下である。第三曲線が非平面なら次は元の中心を通る平面の曲線で、直線には戻れない。[三次曲線に続く分岐](sections/rank_ten_cubic_after_center.tex)

階数10の四次元分岐中心の第二曲線が二次なら、第一直線上の次の順中心は特異三項鉛筆に限られる。第三曲線とその次は平面に入る。特異三項中心の射影同定から三次曲線の後の次数にも上界がある。[二次曲線と特異三項鉛筆](sections/rank_ten_conic_after_center.tex)

二直線列を持たない曲線分岐では、中央支持2の順・逆中心は存在しない。全ての動く成分の一般中央支持は4となる。二つの位相で全可動成分が二次錐に入ることを使う。[二項中心の二次錐](sections/binary_centers_quadric_cones.tex)

階数10の四次元中心からは高々三辺で三次元中心に達する。第二曲線がねじれ三次なら第三曲線は直線または二次曲線で、二次なら第四曲線は直線である。[極限の核直線](sections/rank_ten_flat_kernel_lines.tex)

階数10で中央支持4の分岐中心の第二曲線は必ず平面二次である。二つの核直線族で辺の非被約因子の次数が次辺の小行列式の次数を超えるため、ねじれ三次は生じない。[核直線族の因子による排除](sections/rank_ten_cubic_kernel_divisor.tex)

階数10の分岐では全列が第二二次の場合と、全列が (1,3,1) で正則三項鉛筆を両端に持つ場合に分かれ、両者は混在しない。[三次スクロールによる分岐列の二分](sections/ternary_center_cubic_scroll.tex)

階数10の (1,3,1) 分岐は五次欠損の因子から排除できる。二次曲線列でも特異三項中心の反対消去点が二位置先の順中心に一致し、ほかの中心の支持が2以下に落ちて矛盾する。従って同階数の曲線分岐はなく、全次数判定の未解決候補は有限辺だけである。[五次欠損と二中心周期](sections/rank_ten_quintic_centers.tex)

有限辺・全中央支持4で単純点十八個と二重点一個なら、反対階数10で片側対称となり、階数11は生じない。二重点を印に取り直すと単純印での割線の一致を対称性に変換できる。[一つの二重点の判定](sections/finite_one_double_five_space.tex)

十七単純点と一つの曲線状三重点でも、全中央支持4なら階数10で片側対称、階数11は不可能となる。曲線状有限窓の単純印の五空間には一般に母線がない。[付値と三重点の証明](sections/finite_curvilinear_pencils.tex)

支持が十八点以上の有限窓全体で、全中央支持4なら階数10の片側対称性と階数11の排除が成立する。[合流割線の平方関係と二重点での判定](sections/finite_two_double_five_space.tex)

片側対称な有限窓の対称三平面の和が四次元なら、指定交換子の判別三次式は直線と二次曲線の積になる。三交換子の復元を一つの巡回行列と線型系まで縮約した。[判別式と復元の証明](sections/marked_three_commutators.tex)

十七支持の曲線状有限窓でも階数10で片側対称となる。階数11も排除できる。[埋込み接点の長さと三つの一次方向](sections/finite_seventeen_site_pencils.tex)

有限窓の支持が十七点以上、全中央支持4なら、局所環の型を問わず階数10で片側対称となり、階数11は生じない。[平面四重点の証明と母線の判定](sections/finite_planar_four_seventeen.tex)

既約な三行列の三独立交換子が固定二平面を固定直線に送る場合を排除した。片側対称有限窓の三交換子枝には常に可逆な交換子があり、四次元和の共通二平面は可換にならない。[入射・非入射の排除と判別式](sections/marked_commuting_intersection.tex)

可逆交換子を含む鉛筆では復元の線型階数は3か4となる。可逆な例外方向の全解を二つの平方根と一つの自由係数で表示した。[鉛筆の階数と例外の復元](sections/marked_bracket_pencil_rank.tex)

一般方向にも復元階数3が残る三平面を二行二列の非対称行列で分類した。この例外は四次元和・五次元和のどちらでも有限辺を持てず、一般方向では線型解が一意になる。[一般復元と例外模型の証明](sections/marked_persistent_bracket_defect.tex)

任意に指定した判別直線と跡零二次曲線に対し、接線・重複因子・非半単純型も含め、幾何学的一般点が全次数で正則な階数10の有限窓を構成した。両因子を動かす族の同型類の像は十二次元である。[正則変形と例外三型の検算](sections/marked_prescribed_discriminant_as.tex)

半単純な判別法線と任意の直線を指定した族では、印を忘れても一般最小周期は4で、各位置の射影自己同型は自明である。[二十点配置の変形による証明](sections/marked_semisimple_period.tex)

有限混合窓と指定判別式の全構成を既知の十五次元成分に位置付けた。混合点から階数16へ出る明示曲線では、標準正則性が壊れる四値を決定した。[局所環・成分・特殊値の証明](sections/marked_mixed_component.tex)

混合窓の横断曲線の添字保持同型パラメータは $x=t^2$ で、除外値は $-3,-1/3$、最小周期1は $x=1,-1$ に限る。階数14の点では中心文字を持つ一周期表示を得た。[全同型と周期](sections/marked_mixed_curve_isomorphisms.tex)

混合印の二次鉛筆の四基点を通る全三次式を判別式とする正則変形を構成した。分裂四基点と滑らかな三次曲線を指定すると、相手三平面の全解はその曲線の次数二Picard多様体である。[判別三次式と楕円五次曲線](sections/marked_cubic_linear_system.tex)

指定三次式からの相手三平面を、特異・可約曲線と非被約四基点を含め二本のスカラー方程式で全復元した。[核の二列と全特殊値](sections/marked_companion_plane_reconstruction.tex)

計量も固定した滑らかな三次曲線では、独立交換子を持つ有限混合印が存在するための必要十分条件を三つのVeronese曲面で与えた。例外外の各曲線で正則変形の同型類の像は八次元である。[四割線判定と実在する例外](sections/marked_fixed_metric_elliptic.tex)

滑らかな判別曲線の可換二平面側でも、相手三平面をすべて一本の二次曲線に復元した。剰余点が三次因子に重なるときと離れるときの境界も決まる。[可換側の全ファイバー](sections/marked_commuting_plane_fibers.tex)

滑らかな判別式を持つ三交換子の全組を、楕円曲線の負二倍写像の四逆像と各六次元線型空間で復元した。固定計量例外でも対称三平面の基底を選べば持上げが存在する。[三交換子の全復元と四つの選択](sections/marked_elliptic_commutator_doubling.tex)

基底まで固定した三交換子を滑らかな判別式の全標的で解いた。一般に八解、例外では六解または六孤立解と一つの滑らかなアフィン二次曲面である。[四つの余因子式と全ファイバー](sections/marked_fixed_commutator_adjugate.tex)

滑らかな判別式を持つ片側対称印付き表示の全体は既約で、すべての交換子持上げが既知の十五次元成分に入る。計量例外の各固定曲線では位数4の四点で表示成分が分かれる。[全持上げの成分と特殊化](sections/marked_smooth_discriminant_components.tex)

節点・尖点を持つ整三次曲線について、基底を固定した三交換子の全逆像が特異点の階数と余因子条件で決まる。[特異整三次曲線の全逆像](sections/marked_singular_integral_commutators.tex)

直線と横断二次曲線の全交換子持上げを、余因子と高々二つの階数一補正から復元した。一般ファイバーは四孤立点と二つのアフィン二次曲線である。[横断する二成分の全逆像](sections/marked_transverse_conic_commutators.tex)

接する二成分・三角形・共点三直線を含む被約判別式の全交換子持上げを復元した。接する一般ファイバーの二つの一次元部分はアフィン直線となる。[全被約型の復元](sections/marked_reduced_cubic_commutators.tex)

非被約型も含め、非零判別式を持つ順序付き三交換子を全て復元した。有限窓の片側共通三平面では指定対称側に対する全交換子方程式を覆う。[非被約判別式の復元](sections/marked_nonreduced_cubic_commutators.tex)

両側一般・反対階数16の十五次元成分では、一つの交叉交換子テンソルから窓が復元される。八固有直線の順序付き積表で、開部分全体の添字保持同型を有限判定できる。[高階数窓の復元と全同型](sections/marked_mixed_bracket_reconstruction.tex)

十五次元成分で全二十印と四開始位置を分離した。九語の有理積表で無印同型を判定でき、合同近傍の全持上げが全次数正則である。[二十印の有理不変量と全次数の実現](sections/marked_twenty_rational_invariants.tex)

階数五の交叉交換子の軸を復元し、全持上げの正則性・最小周期・無印同型を決定した。同じ交換子のもとで非正則な二候補と全次数正則な二候補が共存する特殊値もある。[五次元像の全ファイバー](sections/marked_five_dimensional_bracket_fiber.tex)

四元群次数を持つ全係数非零の印付き窓で、主枝の全次数正則性が単一条件 $u^6\ne1$ に一致する。対角正規化の欠損・付加中心枝・無印同型も扱った。[四元群混合トーラスの全閉点](sections/marked_klein_mixed_torus.tex)

四元群の対角交換子非零枝は零係数も含めて全て判定できた。疎な十八の場合のうち、全点正則の三つ、空の二つ、有限次数で排除される十三に分かれる。[全台の表と証明](sections/marked_klein_mixed_supports.tex)

対角交換子が零の疎な四元群窓では、指数増大商で全候補を一つの台に絞り、残る二変数表示の全点で左右AS正則性を証明した。[正規平方と全次数基底](sections/marked_klein_zero_bracket_boundary.tex)

四元群次数を持つ印付き交換子窓では台・符号・全特殊値のAS判定が揃い、正則点の周期は全て2を割る。[対角差零の全係数非零窓](sections/marked_klein_zero_bracket_torus.tex)

四文字の置換に支えられた全窓は、標準正則なら全特殊値を含め量子窓に一致する。単純印を仮定せず、二つの欠色を持つ窓でも未処理の対角台をグラフ三種に絞った。[置換台と欠色](sections/squarefree_window_coverage.tex)

無印の三色階数2枝では、座標の異なる二点を結ぶ辺がなければ、反復対の一成分でも消える全窓が三角PBW極限を持ち正則である。相異なる四文字の成分の消失にも仮定を置かない。[二直線配置とPBW極限](sections/unmarked_klein_pbw.tex)

両切断に欠色がある四元群窓を、印・台を仮定せず全て判定した。文字別Euler漸化式で階数を強制し、零行・零列は全て指数増大となる。[欠色二切断の全被覆](sections/klein_missing_color_coverage.tex)

片切断だけが欠色を持つ場合、反対の $(3,1,1,1)$ は全て指数増大し、$(2,2,1,1)$ は十個の明示台に絞られる。[片側欠色の候補台](sections/klein_single_missing_supports.tex)

片側欠色の十台の方程式を全て解き、特殊値も含め九台を四次で排除した。最後の一変数曲線の全点はPBW基底を持ち正則である。[唯一の正則曲線](sections/klein_single_missing_windows.tex)

片側欠色曲線の無印添字保持同型類は $X∈k^*$ の各点に対応し、一移動は $16/X$、最小周期1は二点に限る。[同型・周期と六直線](sections/single_missing_curve_isomorphism.tex)

両切断で三色が階数1なら、全ての正則窓が対称完全交叉かその二コサイクル捻りになる。[全台の排除と反対階数の特殊化](sections/klein_three_lines.tex)

二切断の色階数が $(3,1,1,1)$ と $(2,2,1,1)$ の組となる正則窓は存在しない。[例外値と疎な台を含む証明](sections/klein_mixed_positive.tex)

二重色の方向が異なる稠密窓では、48成分から正則な二変数族が残る。全点の正則条件は $u\ne v$、最小周期は1である。[二次完全交叉への変換と例外値](sections/klein_transverse_dense.tex)

二つずつの二重色を持つ疎な窓の全台から、平行方向148代表、異方向51代表を抽出した。[軌道表と全数検査](sections/klein_two_doublet_supports.tex)

異方向の二重色の窓は全て、射影直線上の異なる二点の族で被覆される。疎な二表示は無限遠の境界を尽くす。[全63成分と特殊化](sections/klein_sparse_transverse.tex)

異方向二重色族では、反対階数16での無印同型は二点の交換と同時符号だけである。階数12の二曲線も例外値を含めて分離した。[全同型判定](sections/transverse_doublet_isomorphism.tex)

平行な二重色で薄い関係が片方の文字対に限られる全67台の正則性が決まった。実現する14台の全係数はPBWである。[文字順と例外値の全検査](sections/klein_parallel_flags.tex)

三つの薄い色を持つ対称族の無印同型は、四係数の順序交換と、二対反復直線の一つの分数変換で尽くされる。[六節点と例外の五平面](sections/three_lines_symmetric_isomorphism.tex)

三薄色の捻りの無印同型は、全例外を含めて四係数の多重集合で決まる。三薄色の二族・異方向二重色族・一欠色曲線の間に無印の重なりはない。[楕円曲線による復元と枝間の分離](sections/three_lines_twisted_isomorphism.tex)

平行二重色の三台が与える対称・捻りの二曲面で全点の無印同型を決定し、異方向族の無限遠境界を捻り曲面の対角線に貼り合わせた。[例外の二次式と貼り合わせ](sections/parallel_symmetric_surface.tex)

平行二重色の台では、53台について全係数を排除し、未決定部分を持つ台を16台まで絞った。新しい対称曲面と三薄色族も例外を含めて分離した。[文字別の経路次元と六節点の網](sections/parallel_excluded_shapes.tex)

可逆交差ブロックを持つ疎長方形二形の全係数を解いた。主成分は次元3と4で一般に周期4、四つの特別成分では不正則集合を一次条件で決定した。[補完の式と全特殊値](sections/parallel_invertible_rectangles.tex)

第三の疎長方形の全係数を三変数で表示した。四次中心商と長さ14の自由性証明により、三主成分のそれぞれに全閉点で正則な開集合を得た。[標数零への持ち上げと中心](sections/parallel_third_rectangle.tex)

第三長方形の二零条件の曲線で、全点の正則性・無印同型・最小周期を決定した。計64個の残余根も中心商の検算に含めた。[全点の証明](sections/parallel_double_zero_curve.tex)

第三長方形の三変数族全体で、文字座標を指定せずに同型類と最小周期を決定した。一零条件の面は二零条件の曲線に沿って局所的に二枚が交差する。[全点の不変量と特殊化](sections/third_rectangle_orbits.tex)

可逆交差長方形の七つの台で、正則になる全係数を特定した。主成分の二切片では七次の次元が少なくとも134となる。[全係数の判定表](sections/parallel_support_fibers.tex)

台143の二曲線で、正則性・指定始点の同型・周期を全点で決定した。不正則な二根では四次元数が36となり、別の二根は周期1の一類に合流する。[二曲線の完全分類](sections/parallel_septic_curve.tex)

台145の二曲面の無印商は二次錐の開部分であり、二位置移動は二座標を交換する。不正則な因子を特定し、その外では四次中心の基底が全点で延びる。[曲面の同型・周期・中心](sections/parallel_rectangle_surface.tex)

台145の両曲面では積が負の単位となる因子だけが不正則である。残る全点を中心商で検証し、二次錐の商上で正則な同型類を決定した。[曲面の全点の証明](sections/parallel_rectangle_surface_regularity.tex)

第三長方形の全点が正則である。奇数始点で必要な二つの非零条件の和は単元になる。三変数族では四元群の印も全点で一意である。[全点の証明と安定群](sections/third_rectangle_regularity.tex)

台144の一般の添字保持同型類はPicard数20・判別式8の楕円K3曲面で記述される。台146,147の無印商は有理的である。[商と特殊化](sections/parallel_horizontal_moduli.tex)

台82,85の特別成分の無印同型を全点で決定し、二形に共通する不正則因子を求めた。[商と四次障害](sections/first_rectangle_orbits.tex)

第二長方形の特別成分の八代表台で全同型と周期を決定した。対角では符号により有限群と一次元群に安定群が分かれる。[全点の商](sections/parallel_vertical_orbits.tex)

水平台85では四次が標準となる全点について、トーラス三変数が同型を一意に分類する。[全点の枠、周期、低次中心](sections/first_rectangle_horizontal_frames.tex)

水平台85の五つの反対階数とその特殊化を決めた。階数10の二面の全点が四次中心上有限で、正則性を長さ13までで必要十分に判定できる。[階数と中心有限性](sections/first_rectangle_rank_strata.tex)

水平台143から147の全点で四次・五次の判定を得た。第一長方形の階数10面の曲線では26個の数体も含めて全正則性を確認した。[水平台の証明](sections/horizontal_rectangle_low_fibers.tex)

台85の階数10の二面とその符号捻りは全点で正則である。五つの中心元上の表示を使い、開集合と残る四数体で全始点を覆った。[証明と計算](sections/first_rectangle_edge_regularity.tex)

第二長方形で反対階数14を与える交面はPicard数20、格子判別式16の楕円K3である。一般無印商は有理曲面で、周期2の分岐曲線も有理である。[層・曲面・正則開集合](sections/second_rectangle_opposite_ranks.tex)

第一長方形の反対階数14曲線を八数体も含めて全点で正則と確定した。第二長方形の五層全域では五つの四次中心を同じ文字移動で求められる。[証明](sections/rectangle_quartic_bundles.tex)

対称性を外した変形を計算し、三長方形の水平族と正・負の垂直族が、それぞれ同型を除いて次元15、14、4の相異なる既約成分に入ることを決定した。[成分と特殊化の証明](sections/rectangle_unrestricted_components.tex)

第二長方形の階数14・周期2曲線は例外数体を含む全点で正則であり、指定始点の無印商は十二点を除く射影直線である。[全点の証明](sections/rank_fourteen_conic_complete.tex)

負の四次元成分では各辺が一般に楕円四次と四孤立点のグラフで、最小周期は2である。有限点対応を持つ窓はこの成分に属さない。四元群を固定した両対称表示は全ての疎な台を含め既約である。[点配置と特殊化の証明](sections/rectangle_component_point_geometry.tex)

負交換の $(a-d)(zb-c)≠0$ の全正則点で、二円錐曲線や四角形への退化も含めて無印同型と安定群を決定した。商は有理四次元で、周期1の二平面が商上の一因子になる。[商と全周期の証明](sections/negative_pairing_unmarked_quotient.tex)

負交換の全係数表示で例外を含む全正則点の無印商、周期、反対階数を決定した。二重二次曲線では埋込み点から枠を復元し、二次曲面境界では既知の一変数枝への写像を決めた。[境界の証明](sections/negative_pairing_boundary.tex)

疎な台を含む全四元群両対称型で、反対階数と二つの射影入射だけによる全次数判定を得た。負交換の $z=-1$ の商は重み付き射影三次元の明示開集合、$z=0$ は第二長方形の負垂直表示となる。[端点の証明](sections/negative_exchange_end_values.tex)

第二長方形の負垂直族の無限遠を含む全射影商を求めた。錐頂点を保つ加法を除いて全同型を比較し、滑らかな二次曲面枝が無限遠で既知の交代二次錐に特殊化することを示した。[射影境界の証明](sections/zero_exchange_projective_quotient.tex)

薄い全係数が非零で交差ブロックの階数が一の両対称型では、全正則窓と捻りは周期1で、負成分に別の三次元因子がある。[特異角の分類](sections/klein_singular_exchange_corners.tex)

交換対称性を仮定しない $z=1,\lambda^2\ne1$ では、全正則候補を一つの有理超曲面とその移動に絞った。四次経路の障害、五次元族への平滑化、一般周期4を証明した。[混合補完](sections/klein_mixed_exchange_resolution.tex)

薄い二切断がともに退化する繊維で、対称性を仮定せず全正則窓を二つの型とその捻りに分類した。第二型の負表示の二十点は例外値でも全て相異なる。[同時退化の分類](sections/klein_double_exchange_fourier.tex)

平行二重色で四つの薄い長方形に欠損がない全正則窓を、五次元楕円族と正負の交換族の閉包で覆った。最後の可逆繊維には一般周期4と2の階数10の層がある。[全被覆の証明](sections/full_thin_parallel_cover.tex)

等行列式の第二補完では全評価障害が八平面で尽き、残りの二十点と長さ三の衝突を明示できる。四元群の印を復元して一般の無印商を求めた。[第二補完の八平面](sections/equal_determinant_second_sheet.tex)

行列式比が異なる円錐曲線族でも評価で排除される全点を決定し、残る二十点の配置、三重衝突、三次元の有理無印商を求めた。[行列式比の円錐曲線族](sections/unequal_determinant_conic.tex)

点と二平面の巡回から全次数の正則性を証明した。二つの階数10表示で巡回を関数体上で検算し、最後の可逆な薄い繊維の全閉点で正則性が決まった。[点と二平面の巡回](sections/point_line_flag_cycles.tex)

反対階数16を通る七次元の対合成分を構成し、五次元主枝の全閉包で純関係がなく四辺の評価が適合する全点の全次数正則性を証明した。[対合の四点表示](sections/involution_quartet_component.tex)

水平台85と143–147の全点で正則性の必要十分条件を決め、第三長方形の台88–90は全点が正則と証明した。[三つの水平長方形の全点](sections/horizontal_rectangle_flags.tex)

七次元対合成分では四点列を順序付けた商が有理で、無印の一般安定群は位数2、最小周期は4である。[有理座標と安定群](sections/rational_ordered_quartet.tex)

十五次元成分の一般窓は、どの印でも二平面の積内の二枚の三次曲面から復元できる。交点は長さ27で、19個の点列のほかに長さ8が残る。[三次曲面対と復元](sections/marked_cubic_surface_pairs.tex)

有限な一辺に単純点があれば、残る長さ19の配置だけで一辺を復元できる。全補完は48未知数の線型核で表され、十五次元成分の一般点では一意である。[十九点からの一辺と補完](sections/nineteen_point_edge.tex)

十五次元成分の印付き窓は、双三次超曲面の71節点を保つ変形の一成分を支配し、一般ファイバーは有限である。[71節点と局所分類](sections/marked_bicubic_nodes.tex)

五次元主族で両交差ブロックが可逆なら、全点で左右評価の下界3が正則性に必要十分である。除外される二型には指数増大商がある。[判定と例外表示](sections/elliptic_evaluation_obstructions.tex)

混合交換境界の正則性は一式で決まる。これで四つの薄い長方形に欠損のない二重色窓の全てに全次数判定が揃った。[混合境界と被覆](sections/mixed_exchange_regularity.tex)

四元群が各位置で正則表現として作用する全窓で、正則性・無印同型・周期の有限判定が揃った。安定群が正次元でも印の共役類は有限である。[全枝と無印不変量](sections/klein_unmarked_atlas.tex)

全正則表現型四元群窓は半安定であり、唯一の閉軌道も正則である。疎な平行17台の閉軌道は二つの台に集約される。[閉軌道と商の正規化](sections/klein_polystable_limits.tex)

四元群を仮定しない全窓で、トーラス付きEuler行列式は実際の四文字に分解する。各位置が二重みなら七型だけが残る。[指標公式と全候補](sections/torus_euler_characters.tex)

各矢が二平面に分かれる窓の可逆部は全点が二次PBWである。全表示と同型は二つの逆行列の移送と安定群の共役に帰着する。[全正則性と分類](sections/split_inverse_windows.tex)

二平面分解で対角化可能な全スペクトルを曲面・直線・六次式で分類した。六次式の二軌道は印を忘れると一つになり、周期1・安定群が乗法群だけの正則窓を与える。[スペクトルと孤立型](sections/split_bell_locus.tex)

二平面可逆部では四固有値の条件が全Jordan型にも必要十分である。全補完が同じスペクトルの半単純補完に退化し、閉軌道も正則である。[重根と退化](sections/split_jordan_transport.tex)

どの標準窓でもSL二群の矢表現は四位置で一致する。非自明四型のうち三型の全正則窓は反対階数6・8・12の三つの周期1類に尽きる。[分類と残る型](sections/sl2_window_profiles.tex)

非自明な連結半単純対称性を持つ全正則窓を印なしで分類した。全て周期1で、固有の多重安定極限も正則である。[正規形と特異境界の排除](sections/one_standard_sl2.tex)

各矢に高々二重みの連結トーラスを持つ全窓を分類した。四重み差と重複度は一致し、三対一では中心多項式拡大、二対二では八純関係全てが可逆である。[分類定理](sections/binary_torus_windows.tex)

有限点グラフの各辺のトーラス階数は高々一で、非自明な窓作用を持つ点配置は長さ十の二点、または両端が長さ九で中間が長さ二の二型に限る。[重みと局所長](sections/finite_torus_profiles.tex)

長さ十の二点を持つ四単重みのトーラス窓は有理曲面と有理曲線で尽き、ただ一点で交わる。周期と有限グラフの二つの判別式を決定した。[式と正規形](sections/simple_weight_finite_windows.tex)

有限点グラフと非自明なトーラス作用を併せ持つ全窓の分類が決まった。中央二重型は反転を除く一曲線で、全点が周期一である。[鉛筆・正規形と例外値](sections/double_middle_finite_windows.tex)

有限トーラス窓を閉極限に持つ全窓を、対称性を仮定せず重み付き射影空間で分類した。二つの異なる次元の全体成分が周期一の窓で交わる。[局所成分と全ファイバー](sections/finite_torus_slices.tex)

周期一の共通窓で、全軌道横断局所環の方程式は二本の積に尽きる。76次元と75次元の成分および共通部分は全てここで滑らかである。[交叉の局所環](sections/torus_component_crossing.tex)

対称性を仮定しない吸引ファイバーで、重み付き座標の全ての消失場合を含めて最小周期を決定した。[周期一・二・四の判定](sections/torus_fiber_periods.tex)

有限点対応を保つ加法群のJordan型を制限し、四位置が推移変換なら全次数で正則性を証明した。[加法群と閉極限](sections/finite_additive_windows.tex)

正次元の安定群を持つ有限Frobenius窓の全てについて、全次数の正則性と半安定性を証明した。[Jordan型と各還元先](sections/finite_connected_regularity.tex)

有限点グラフの加法群対称性を全て排除した。正次元の安定群がある有限候補はトーラスの一曲面と二曲線に尽きる。[固定点の接空間による排除](sections/finite_unipotent_exclusion.tex)

量子点のうち、有限点グラフの軌道閉包に入るのは正確に二類である。通常の多項式環は二次障害で除かれる。[量子閉極限の必要十分条件](sections/finite_quantum_closed_limits.tex)

量子点へ退化する有限窓の残る全てを、六次元アフィン空間一枚と直線を除いた二枚で分類した。[六座標・例外直線・全周期](sections/finite_quantum_planar_fibers.tex)

任意の有限点グラフで四文字評価は階数20、両対称の場合だけ19になる。通過成分の次元下界も得た。[部分スキームへの制限と余法線](sections/finite_four_slot_evaluations.tex)

有限Frobenius窓はすべて半安定であり、厳密半安定な候補は例外なく全次数で正則となる。[重み不等式・等号条件・PBW還元](sections/finite_window_semistability.tex)

有限窓の安定性を、長さ9以上の局所点から復元する高々二点の旗で判定した。[局所環による必要十分条件](sections/finite_intrinsic_stability_test.tex)

非正則な有限両対称窓の全体を二元十次形式と八次形式から復元した。一般の例は安定で相異なる20点を持つ。[割線被覆と判別式](sections/binary_secant_obstructions.tex)

二元商を持つ窓では、二行三列の階数だけでFrobenius性が全例外を含め決まる。商の向きを固定した同型も求めた。[差分式と商の復元](sections/binary_apolar_rank_test.tex)

対称性のない四窓でも、二辺に二元商があれば中央の三次曲線は自動的に揃い、非正則な四窓全体が一意に決まる。[非被約交叉を含む復元証明](sections/binary_two_edge_rigidity.tex)

二元商の二方向が重なる条件は、二次式と四次式が一次因子を共有することに等しい。逆向きの表示と最小周期も計算できる。[十四次因子と座標変換](sections/binary_reverse_orientation.tex)

安定な非正則窓の一般点を通る全成分の局所環を決定した。14次元と15次元の成分が13次元で交わり、15次元成分の一般点は最小周期4の正則窓である。[形式局所環と全次数の証明](sections/binary_component_crossing.tex)

15次元成分で反対行列式は片側対称な二因子の各六乗と滑らかな残余因子の積になる。一般階数は16、10、残余では15である。[局所因子と階数の証明](sections/binary_opposite_determinant.tex)

二成分の13次元の交叉を、二つの束の差が非自明な二捩れとなる楕円鉛筆の閉包と同定した。[四対角の切断と二角形からの平滑化](sections/elliptic_two_torsion_crossing.tex)

有限Frobenius窓では一辺の二元商だけで指数成長と反対位置の片側対称性が従う。長さ29で標準Hilbert関数を超える。[核束の分裂と三文字ブロックの証明](sections/single_binary_edge_obstruction.tex)

一辺だけの二元障害を持つ有限窓を、一つの直線複体の族にまとめた。反対階数10、9、8を全例外で分け、非両対称窓の最小周期4を証明した。[直線複体の階数層](sections/binary_linear_complex.tex)

新しく局所的に得た15次元成分を既知の15次元成分と同定した。単一二元辺の非正則族と階数15残余因子の区別も確定した。[滑らかな印付き三次による大域的同定](sections/binary_global_component.tex)

階数15残余因子の一般点では左右の核がともに非退化で、添字保持同型は二つの内部共役に帰着する。[核の四次計算と共役判定](sections/rank_fifteen_opposite_kernels.tex)

二方向の二元点で残余行列式の四次接錐を求め、階数15残余因子の一般点と指定の片側対称切片の一般点が全次数で正則であることを証明した。[四次接錐からの全次数判定](sections/double_binary_residual_regularity.tex)

点グラフが正次元でも、一辺が二元のHankel六空間に同型に射影する窓を全て排除した。[基点を持つ評価層の全場合](sections/binary_edge_basepoints.tex)

二元Hankel空間への射影が真部分空間になる場合も、正則候補の射影階数を4・5、隣接評価の捩れ長を2・3に限定した。[非全射な二元射影の必要条件](sections/deficient_binary_support.tex)

一辺の射影が二元Hankel六空間に含まれる全窓を、射影階数・点グラフの次元・重根の有無を問わず排除した。[両側の点ファイバーによる排除](sections/deficient_binary_quintic_pencil.tex)

任意の正則辺に含まれる純関係の正次元部分は、各向き高々一本の母線、または唯一の滑らかな二次曲線に限る。後者に加わる孤立純関係は高々三点である。[純関係の曲線の支持](sections/pure_relation_curve_supports.tex)

一点ファイバーの曲線グラフでは一辺の純関係は高々四点である。四点のとき曲線部は五つの座標直線配置に限り、各配置の有限長部分も決まる。[長方形切断と四純点の座標枠](sections/curve_graph_pure_frames.tex)

曲線グラフのAS正則窓では純関係は各辺で被約な一点以下である。一点あれば全曲線成分は全位置で直線となり、反対辺の直線グラフの縮約にこの純関係が現れる。[四次鉛筆と次数零の束による排除](sections/curve_graph_single_pure.tex)

一次元点グラフの正則窓には純関係がなく、すべて半安定である。直線成分があれば純関係の不在で正則性が決まる。四移動で保たれる直線上で両評価が階数3なら、外側に分岐があっても四本の二文字関係で全次数のPBW性と正則性を判定できる。[曲線グラフの排除と直線基底](sections/curve_graph_pure_exclusion.tex)

一次元点グラフでは純関係のない厳密半安定窓はすべて正則となる。四移動不変直線上で評価階数3の正則窓の閉軌道極限は、二組の相補直線を持つ正則な二次PBW窓である。[等号方向と閉軌道極限](sections/curve_graded_limits.tex)

相補直線を持つ任意の正則窓への指定重みの全持上げは一次条件と二次障害で決まる。四本の法線関係の合成が恒等なら障害空間は五次元、最後の係数に三自由度がある。それ以外では障害空間は二次元で最後の係数は一意となる。[全持上げの式と指定付き同型](sections/invariant_line_filtered_lifts.tex)

六次方程式の例外窓は、四テンソルの全256係数を動かしても剛直である。射影軌道の閉包は独立の59次元成分をなし、半安定商は被約な一点である。この窓に特殊化する別の標準同型類はない。[全接空間と標数零の証明](sections/isolated_split_component.tex)

一般ファイバーで同型な二つの正則模型では、四文字の行列式付値の和が最上テンソルの付値の四倍になる。多重安定な正則窓から別の正則窓への特殊化は起こらない。六次式の59次元成分の正則類は、不安定境界も含めてこの一類で尽きる。[全次数格子の比較](sections/regular_orbit_valuations.tex)

可逆二平面分解を持つ窓とその正則な全持上げは全点で分岐を持たず、二つの四行四列鉛筆で点スキームが決まる。六次式の孤立正則窓では、二本のねじれ直線と三本の結線が六点で交わり、さらに長さ一の埋込み部分が二点にある。[点グラフの式と代数体上の検算](sections/split_point_pencils.tex)

Bell曲面上の四補完の一般点に特殊化する全窓を、対称性や直線を指定せずに分類した。閉軌道以外は重み1が四個、重み2が三個の重み付き射影空間の四元群商で、曲面を動かすと同型類は八次元である。[全二次障害の消滅と無印の比較](sections/split_surface_full_fibers.tex)

Bell曲面の一般の各閉軌道に至る全窓の最小周期を決めた。非閉窓で周期2の部分は次元4,2,2,2の四既約成分で、それ以外は全て周期4となる。[全座標での固有空間判定](sections/split_surface_periods.tex)

Bell曲面の四補完はすべて既知の十五次元成分に属する。一般の閉軌道とその全流入窓で成分は滑らかで、全点が正則で中心上有限となる開近傍がある。[五中心・量子点への特殊化・全接方向](sections/split_surface_component.tex)

十五次元成分の正則窓では分岐を通る周期は長さ12以下で、二周期間に二文字以上の切替経路があれば周期の組は(4,4)か(4,8)である。一般の正則窓でも分岐周期では正次数の中心元が全て消える。[行列恒等式と周期点加群](sections/pi_periodic_point_paths.tex)

Bell直線の半単純部は全量子類と四つの非量子行で尽き、例外値を含め全点が最小一周期である。[全補完と添字保持同型](sections/split_spectral_line_classification.tex)

Bell直線から十五次元成分に入るのは位数四の量子窓と反転なしの一点の二類である。片側反転族は別の射影59次元既約成分をなし、その商は有理曲線である。[成分交差と全中心跳躍](sections/split_line_component_incidence.tex)

片側反転成分では不安定境界を含む全正則点を分類した。正則類はアフィン直線をなし、無限遠は全て非正則である。全量子族にも三角積が極を持つ正則特殊化はない。[経路格子と動く境界](sections/pbw_valuative_boundaries.tex)

量子点を閉軌道に持つ三表の和は、背景と変形係数を同時に動かしても正則特殊化で閉じる。特殊窓の量子三角積は一般背景の剰余で決まる。[全量子流入軌跡の閉性](sections/quantum_regular_closure.tex)

二固有値のBell直線で非半単純な全窓を、量子三角積と一辺または二辺、および片側反転に流れる例外一類で分類した。全て最小一周期である。[Jordan窓の無印同型](sections/split_two_eigenvalue_jordan.tex)

Bell直線の単一固有値の全Jordan分割について、無印同型・周期・閉軌道を決定した。これで直線上は非半単純部を含めて尽くされる。[単一固有値の全分類](sections/split_single_eigenvalue_jordan.tex)

Bell曲面の固有値衝突部は全Jordan型と例外値を含む三行で尽きる。最小周期・無印同型と十五次元成分との全交差を決定した。[重複曲面の分類](sections/split_repeated_surface.tex)

可逆二平面分解の全窓で無印同型と最小周期を決定した。単純曲面の例外曲線は独立な射影60次元成分を与える。[単純曲面と可逆分解の全被覆](sections/split_simple_surface.tex)

片側反転の一般の閉軌道には他の窓は流入しない。最初の例外では閉軌道とJordan一類だけで、二成分の横断完備環は二つの二次単項式で決まる。[全流入と局所環](sections/mixed_flip_slice.tex)

片側反転の最後の例外への流入は五つの座標空間で尽きる。全点で正則性、同型、最小周期と軌道境界を決定した。[正規形と同型](sections/mixed_flip_fiber.tex)

可逆二平面分解のある閉軌道へ特殊化する窓は、正則な特殊化の後も同じ範囲に留まる。追加単純成分の正則点は境界の二類を含めて決定した。[体積公式・全正則特殊化](sections/split_regular_limits.tex)

通常Bell曲面の単純点と反転なし固有値衝突の全点で、流入同型類の六次元商と全最小周期を決定した。閉軌道が周期1となる衝突の上に五つの周期1成分がある。[全パラメータの証明と式](sections/split_surface_all_fibers.tex)

反転を伴う衝突に至る同型類は一般に五類または三類である。例外値の六次元族まで含め、全窓、軌道境界と最小周期を決定した。[全流入の証明](sections/reflected_collision_fibers.tex)

二平面分解を持つ閉軌道へ特殊化する全窓の被覆・正則性・同型・最小周期判定を得た。反転なし直線の周期1軌跡は通常八成分、特殊値五成分である。[全アトラスの適用範囲](sections/noflip_line_fibers.tex)

二直線列・曲面射影・純関係のない分岐では、三次曲線の交叉長から長い二配置を排除した。残る二列はともに次数1,3,1で、同時に38係数の線型表示に入る。[交叉の証明と同時座標](sections/branch_cubic_contact.tex)

最後の1,3,1型配置も、二切断の補完式から生じる共通四次根により排除した。純関係自由・曲面射影なしの分岐があれば、二直線列が存在し反対階数は12以上である。[全ての四次式の退化を含む証明](sections/branch_pencil_resultant.tex)

分岐を許す全曲線射影でもスキーム次数は7以下、分岐を通る周期点列は長さ28以下に取れる。純関係自由な窓では全位置で各向きの中心は一つ以下で、二直線列の両端の割線はともに接線に限る。[次数評価と接線の証明](sections/branch_degree_seven.tex)

純関係自由な全分岐の次数配置は、連続する二つの二直線列と一つの1,3,1列、または四つの二直線列に限る。始点を保つ有限列挙と三次曲線の射影で証明した。[全列の配置](sections/pure_branch_degree_patterns.tex)

曲線・点射影で純関係のない標準正則窓には、どの反対階数でも分岐がない。最後の二配置を四窓の同時座標で排除した。[両配置の矛盾](sections/final_pure_branch_frames.tex)

分岐には実際の中心の三次元縮約像に純関係が必要である。純関係が有限なときの中心鉛筆を標準形に分け、二直線列の始点への過剰接触が長さ2の純関係を強制することを示した。[中心の局所判定](sections/branch_local_pure_pencils.tex)

全四辺にただ一つの被約な純関係を持つ標準正則窓は存在しない。純関係の個数を問わず、各位相一中心で四つの二直線列が巡回する場合も排除した。[全射影次元での排除](sections/single_pure_branch_exclusion.tex)

全正則窓で純関係スキームは有限であり、分岐中心の中央支持2も存在しない。固定因子を持つ純直線は平面四次因子の全重複度から、純二次曲線は二項中心の曲線列から排除した。[純関係の有限性](sections/finite_pure_relation_locus.tex)

純関係があれば二辺以上に現れ、そのうち一辺は長さ2以上となる。一辺だけの任意長の純関係、および全辺各一点以下の場合を排除した。[辺分布の証明](sections/sparse_pure_edge_exclusion.tex)

二次曲線の平面方程式が純関係の因子になる条件を示し、次数列121と1221の端での過剰接触を制限した。[証明](sections/conic_corridor_excess.tex)

正次元射影を持つ全窓の正則性と添字保持同型を有限判定に帰着した。全ての分岐を排除し、全正則窓に純二次関係がないことも従う。[被覆の証明](sections/unrestricted_corridor_exclusion.tex)

有限点グラフの三文字評価の階数は16,19,20に限られる。対称性のない両側で点列イデアルは四次までの生成で、十五次元混合成分の一般点では四次生成・五次関係である。[証明](sections/finite_three_slot_evaluations.tex)

四頂点切断の点加群の相互Extを反対評価で決めた。二十点配置の分離三次曲面の特異点は四位置に共通の無向グラフを作り、各点の次数は4以下である。[証明](sections/finite_cut_point_pairings.tex)

有限片側対称窓の可換三段加群と極小自由分解を決めた。線型シジジーが部分束を与える領域は電荷4のインスタントン束から復元できる。[証明](sections/finite_symmetric_instanton_monad.tex)

有限片側対称窓のモナドは常に捩れなし層を与える。反対階数10で局所自由でない全境界は二元商を持ち非正則である。[証明](sections/finite_rank_ten_instanton_boundary.tex)

階数10の有限片側対称な全正則窓は、位相ごとに一つの14次元既約閉包に入る。空間楕円八次曲線と四次Grassmann写像から窓を復元し、添字保持同型と二移動を決めた。[証明](sections/finite_quartic_grassmann_windows.tex)

両局所自由な候補では四次K3曲面上の四つの線束から有限Frobenius代数を復元し、四つの球面twistの積を残余な二つの楕円鉛筆で分解した。[証明](sections/finite_quartic_k3_cut.tex)

四つのtwistの積が各滑らかな楕円繊維上の一般の次数マイナス四線束を二移動を除いて固定することを、線型行列式で示した。[証明](sections/finite_elliptic_fiber_transfer.tex)

曲面のNéron–Severi群が超平面と楕円繊維で生成される場合、二つの相対Picard曲面を同定し、四移動の逆を自己交叉4の豊富な線束による捩りにした。[証明](sections/finite_residual_picard_surfaces.tex)

超平面と繊維でPicard格子が生成される四次切断を持つ反対階数10の有限片側対称窓では、二鉛筆への制限から全次数の左右AS正則性を証明した。[証明](sections/finite_minimal_picard_regularity.tex)

全辺有限・片側対称・反対階数10では、二つのモナドの局所自由性が全次数のKoszul左右AS正則性の必要十分条件である。楕円八次曲線表示の全許容値に成立し、曲面のPicard格子に制限はない。[証明](sections/finite_rank_ten_full_regularity.tex)

十五次元の非対称成分の一般窓に一文字を印すと一般の四行四列四層テンソルが得られ、拡張は生成的有限である。行列式K3上の同じ二十点を切る二鉛筆から窓を復元した。[証明](sections/finite_determinantal_marked_slice.tex)

非対称十五次元成分の既知の中心開集合では滑らかな平面四次曲線と四つの半安定階数四束から窓を復元した。逆向きには元の環で中心に持ち上がる切断が二本でも有限判定できる。[証明](sections/finite_genus_three_curve_cut.tex)

十五次元成分の一般の射影中心に滑らかな四次K3超平面切断が存在し、特異点は有限である。既知の二十節点だけでは充満例外列のEuler数と両立せず、頂点次元二の単純が実在する。[証明](sections/finite_central_k3_general_section.tex)

全辺有限の窓で、各頂点二次元の単純表現を四つの階数三・四の二次曲面と罫線指定で分類した。後続は一意または空間三次曲線を含む二次曲面の網となる。[証明](sections/finite_two_dimensional_quadrics.tex)

十五次元成分の閉包で二十点が被約かつ互いにExt直交なら、正則性を仮定しなくても全窓に頂点次元二の単純が存在する。存在しないときは非被約点または分離三次の特異辺が必要である。[証明](sections/finite_two_dimensional_specialization.tex)

十五次元成分の一般中心は通常二重点をちょうど三十個持つ。追加の十点に二種類ずつ頂点次元二の単純があり、一般中心には指数四の不分岐Brauer類がある。[証明](sections/thirty_central_nodes.tex)

中心写像は三十節点四次式の一成分に支配的で生成的有限である。三十点の配置から中心の式が一意に決まり、中心を固定した窓に一次変形はない。[証明](sections/thirty_node_center_map.tex)

十五次元枝の一般中心の不分岐類は周期と指数がともに四である。中心の同定とその類を固定すると、窓の違いは添字の四通りの移動に尽きる。[証明](sections/thirty_node_brauer_windows.tex)

十五次元枝の一般窓は中心の四次式から添字移動・逆転を除いて一意であり、添字保持では八枚を区別する。三十点の二種類は中心の巡回四次Brauer群から読める。[証明](sections/thirty_node_dihedral_cover.tex)

頂点次元二の単純を印すと、矢の全十六階数型を線型空間と残余基底群で記述できる。一般次元は81で、全辺階数三で一周の同定が恒等のときだけ88になる。[証明](sections/two_dimensional_marked_tensor_forms.tex)

中心切断が認証された有限窓で、中心が節点四次式または節点を持つ四次分岐二重二次曲面なら、特殊値も含めて偏極中心と八つのKummer類で同型と周期を判定できる。[証明](sections/nodal_center_special_values.tex)

四つの四次中心パラメータを持つ素な窓では、Hilbert重複度の整性だけで各頂点の一般隅が可除と分かる。正規な種数三中心の特異曲線の一般横断型は二重被覆型を含めてA1に限られる。[証明](sections/central_cut_index_and_transverse_type.tex)

四つの四次中心パラメータに通る有限窓はすべて素で、中心は正規Gorensteinである。中心の次数と頂点指数は(4,4)または(16,2)に限られ、後者の一般超平面は滑らかな種数九K3である。[証明](sections/central_cut_center_dichotomy.tex)

四次中心列を持つ全有限窓の添字保持同型は、偏極中心と零番Azumaya隅の特殊線型持ち上げ障害で判定できる。非孤立特異点も含む。[証明](sections/central_cut_kummer_classification.tex)

頂点指数二の有限窓はすべて通常の周期一Clifford代数である。二次式の四空間の基点無し条件と低階数交叉の有限性が全例外値を含む判定を与え、同型は四空間の射影同値に一致する。[証明](sections/index_two_central_clifford.tex)

十五次元閉包の全有限候補で、正則かつ一隅が整域であることは経路長48の四中心パラメータを使う有限検査と同値である。四次中心列を持たない場合も含む。[証明](sections/finite_domain_specialization_parameters.tex)

十五次元閉包で有限辺を持つ正則点は、一隅が半素なら既に整域枝の有限判定に入る。残余の正則点があれば全隅に非零冪零両側イデアルが必要である。[証明](sections/finite_semiprime_specialization.tex)

成分を指定せず、行列恒等式のサイズを固定した左右Noether正則窓を有限な中心切断で判定した。中心の正規性と合成整域性も従う。[証明](sections/noetherian_pi_window_cutoff.tex)

PI正則窓が中心上有限になる必要十分条件を、非頂点単純と頂点単純との両向きの第一拡張で与えた。半素性や成分の指定を要しない。[証明](sections/pi_vertex_extension_criterion.tex)

単純表現の反対二辺または一辺を補う線型空間から、全成分で使える指数成長の障害を得た。頂点次元二で戻り辺が二次元動く場合は二元三次式の偏極に一致し、全て排除される。[証明](sections/simple_linear_completion_growth.tex)

四辺が有限で一つの頂点隅が二行二列行列の恒等式を満たす正則窓を、例外も含めて四二次式系で完全に分類した。中心有限な場合には頂点指数3も存在しない。[証明](sections/one_corner_clifford_classification.tex)

指数4で四次中心が一斉に消える繊維を三型に絞り、非Azumayaの二型では八次中心が消えないことを示した。[証明](sections/index_four_quartic_basepoints.tex)

多項式成長の全窓で、安定な文字作用の周期列は同じ位相の一辺を共有すれば一致する。向き付き二次曲面の原始周期から頂点次元と倍率安定化群を復元した。[証明](sections/periodic_stable_arrow_uniqueness.tex)

有限中心の次数安定化群を頂点指数で抑え、指数4の非自明安定化群を二重周期単純と逆符号対で決定した。十五次元整域閉包で二倍偏極はCartierで、因子上の安定化群は自明である。[証明](sections/central_character_periods.tex)

全成分の指数4以下の有限中心正則窓で八次の中心四元が取れることを示した。十五次元閉包の半素正則枝の検査上限は経路長564から84に下がる。[証明](sections/index_four_octic_normalization.tex)

一隅がPIの正則窓では、左coherent性が右片側イデアルの四次元成長を強制し、左右coherentなら全成分で合成整域となる。[証明](sections/coherent_pi_regular_domains.tex)

十五次元閉包に非整域の正則窓が残るなら、経路長4か8に中心零因子があり、どの頂点隅にも一様有界次数の消去元があることを示した。[証明](sections/finite_short_trace_pole.tex)

特殊経路環から相対trace環への核が冪零であることを示した。残余の正則点は四次または八次の中心冪零元を必要とし、一つの有限次元表現多様体でこの障害を検査できる。[証明](sections/relative_trace_nilpotence.tex)

二次元印の四辺全てで同方向に行解が増える全作用窓は、四つの二次写像の長さ17の固定スキーム上で表示でき、汎関数の最小支持は長さ15以上となる。[証明](sections/directed_binary_defect_cycle.tex)

十五次元閉包の正則整域窓では次数安定化点が孤立する。一般八次曲面の四つの剛な束から窓が戻り、四次系に基点があっても同じ中心上の添字保持同型類は有限である。[証明](sections/octic_surface_fixed_center.tex)

五つの四次中心切断の像を四型に絞り、うち三型で重み付き完全交叉表示を得た。残る三次スクロール型の標準基点は一般標準曲線上の異なる二つのWeierstrass点である。[証明](sections/quartic_center_image_atlas.tex)

三次スクロール型の中心の最小生成・関係次数を半標準曲線から決定した。十五次元閉包の正則整域判定は経路22次・剰余21次・複体24次までに短縮できる。[証明](sections/scroll_spin_center_generators.tex)

任意の成分の中心有限正則窓でアフィン中心の因子類群は有限で、その指数は頂点指数の16倍を割る。射影中心の有理因子類群は一次元である。[証明](sections/finite_center_divisor_torsion.tex)

指数4の可除類の分岐を調べ、四次デルペッツォ曲面像を排除した。滑らかな三次元二次曲面上では非正規・可約な極因子を含め基点は存在できない。[証明](sections/exclude_del_pezzo_center_images.tex)

基点を持つ三次元の四次中心像は階数3の二次曲面に限られる。特異二次式二本は共通頂点を持ち、第2 Hirzebruch曲面上の反標準楕円曲線または節点巡回列で四次分岐する。[証明](sections/singular_quadric_common_vertex.tex)

十五次元閉包の全ての正則整域窓で四次中心パラメータが存在する。各例外を含め中心とKummer印が同型を決め、整域正則性の検査上限は経路14次・剰余13次・複体16次である。[証明](sections/main_domain_quartic_parameters.tex)

四行列恒等式の有限正則窓では二・三次元の均衡単純が頂点と第一拡張を持つ向きは高々一方である。三行列恒等式の場合、一隅半素なら任意の成分で通常Clifford窓になる。[証明](sections/pi_yoneda_bridges.tex)

全四隅が四行列恒等式を持つ有限正則窓では、全非頂点単純が均衡し、根基を割った環は左右Noetherである。[次元型の排除と半素商](sections/four_matrix_balanced_simples.tex)

十五次元閉包の有限四辺では、全正則点が合成整域になる。標準候補の正則性は経路14次・中心剰余13次・Koszul複体16次までの四次中心検査と同値で、正則軌跡は開である。[整域仮定を除く証明](sections/main_regular_specialization.tex)

どの成分でも、有限四辺の全隅が四行列恒等式を満たし一隅が三行列恒等式を満たせば、全正則窓は通常Cliffordである。半素性を仮定しない。[周期尾部を使った証明](sections/periodic_tail_vertex_bridges.tex)

長さ二十の有限辺では、点イデアル層の第m冪に次数2mの非零曲面はない。非被約な辺でも成立する。四行列恒等式の非整域候補では、二次元単純による六次行列式が全位置で恒等的に消える。[多重点不等式と行列判定](sections/finite_edge_symbolic_powers.tex)

全成分の有限四辺で、全隅が四行列恒等式を満たす正則窓は必ず中心有限な合成整域になる。この枝の正則性は経路30次・中心剰余29次・複体32次と一般隅指数で決まる。[trace指標の無限遠と全正則窓](sections/four_matrix_all_regular.tex)

一隅の三行列恒等式だけで有限四辺の正則窓は通常Clifford窓に限られる。他の三隅への行列次数の仮定も不要である。[一隅からの全分類](sections/single_three_matrix_windows.tex)

反対二隅の四行列恒等式だけで有限正則窓は中心有限な合成整域となり、32次までで判定できる。任意の二隅なら全非頂点単純の次元が均衡する。[普遍拡張と二段積の有効次元](sections/opposite_four_matrix_corners.tex)

多項式成長から、任意の非頂点単純の一般次数倍率間で第一Extが消え、頂点単純経由の二段積が単射となる。有限辺・行列恒等式を要しない。[経路代数による証明](sections/polynomial_growth_extension_bridges.tex)

非恒等の平面グラフを基底に持つOre窓では、指定点を忘れた同型類も決定した。[正規形と平坦特殊化](sections/plane_graph_ore_windows.tex)・[検算](computations/verify_plane_graph_ore.py)

正次元の辺でも、三生成部分を四点の階数条件から検出できる。滑らかな楕円基底を持つこの枝では、有限四辺の部分以外はすべて最小周期1である。[全補正と特殊値](sections/elliptic_ore_recognition.tex)・[検算](computations/verify_elliptic_fixed_ore.py)

## 最初に読む

[研究方針](research/plan_ja.tex) → [主原稿](papers/nc-p3.tex)の対象節 → [必要な研究ノート](research/notes_ja.tex)。
Cloud Codexの環境準備・権限・止まりにくい研究ループ・再開の手順は [Cloud Codexセットアップ](docs/cloud_codex_ja.md) と
`bash scripts/cloud_preflight.sh` を参照する。
実行規則は [AGENTS.md](AGENTS.md)、終盤を見据えた分類の方法・nc-p3自体の巨大化に対する三層アクセスと差分検証は [分類アトラス設計](research/classification_architecture.md)、
一次文献と外部依存は [BASELINE.md](research/BASELINE.md)（as4は必要箇所だけ取得）、
主張と課題の索引は [STATE.json](research/STATE.json)、今回の検証範囲は [validation.json](research/validation.json) を参照する。

対象は、標数零の代数閉体上の標準14641型 quadratic Koszul AS-regular indexed Z-algebra。
同値関係は**添字保持同型**であり、graded同型、周期構造を指定した同型、qgr同値とは区別する。

## 成果と限界

| 主題 | 原稿に証明がある結果 | まだ結論していないこと |
| --- | --- | --- |
| 有限四窓 | Frobenius標準4周期、幾何学的同型類の四窓への還元、認証済みAS対象の同型・周期の有限方程式 | 四次の必要階数条件と五次までの正確性の二階数判定は証明済み。一般四窓のKoszul/AS条件の有限な必要十分判定、明示的な正則族による全被覆は未解決 |
| 量子多項式 | 一般の局所成分でgraded同型類6次元、indexed同型類3次元 | 非1周期的対象が増える現象ではない。三つのtwist方向を忘れる比較 |
| Clifford点 | 選んだgraded表示の接空間9、indexed変形の接空間15。全15変数の二次障害は、六つの交代変数上で非退化な二次式 Q 一つ | 成分次元14、高次障害の消滅、全1周期構造の分類 |
| 楕円Clifford–Ore族 | 九変数の実際の代数的族。全ファイバーAS正則、同型類の像9次元、一般の最小周期2 | 九次元像が既約成分であること、全特殊値の周期、完全な族内同型判定 |
| 同型不変量 | 内在的な偏極付き楕円曲線と非零2等分点。次数4の必要不変量 U(th) | Uの一致による十分な同型判定 |
| 境界と特異性 | 五変数の別のOreチャートは同型類像2次元、一般最小周期2。七変数族の軌道閉包に含まれる。閉包 Z7, Z9 はClifford点で特異。指定11方向の滑らかな同時積分は不可能 | 独立した新成分、全モジュライの特異性 |

「原稿に証明があること」「今回の論理監査」「厳密計算の再実行」「形式化・査読・新規性確認」は別の状態である。
新規性と第三者査読は未確定。具体的な検証範囲は validation を参照する。

### 分類に用いる既知の代数的族

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

最上次数テンソルによる [Frobenius候補の復元](sections/frobenius_tensor.tex) と、
第三の行列階数による全対象の二分岐を証明した。
行列束の全場合分けにより階数2を除外し、有限Frobenius候補の階数範囲を3以上16以下へ絞った。AS対象の下界は6であり、可換多項式環が達成する。
[枝台帳](research/classification_branches.json) の両枝の一般正規形・正則性判定は未完成。
量子部分では階数6・8・10・12の全条件を求めた。さらに [量子部分の全indexed分類](sections/quantum_indexed.tex) として、
全ての非零パラメータを三次元トーラスの明示的なS4作用へ還元した。特殊値も含み、階数6は二種類となる。
対称テンソル部分は四二次式の完全交叉によるチャートで全ファイバー正則性を認証し、
その全特殊値を含むindexed同型を33未知数41方程式へ還元した。
このチャートへの所属は三つの線型解空間内の可逆行列の存在で必要十分に判定できる。
符号を変えた線型判定は可換多項式のindexed類を認識する。有限Frobenius候補では一つの可逆な負符号証明書で十分である。残る低階数領域の非空性・排除は未解決。

可換三生成部分上の [18係数Oreチャート](sections/polynomial_ore.tex) では、
全ファイバーの正則性と `rho = 6 + rank Q` を証明した。`Q` は6×6交代行列で、
低階数部分は可換多項式のindexed類だけとなる。18は係数空間の次元である。

位置依存の可換・反可換Ore四窓でも、三次の線型条件が全AS/Koszul正則性に必要十分となる
（[可換側](sections/polynomial_ore.tex)、[反可換側](sections/anticommutative_ore.tex)）。
反可換側の三つの消せない交代方向は階数8となり、両チャートの低階数部分は既知の二枝へ入る。

さらに [対角量子Oreチャート全体](sections/quantum_ore.tex) でも30本の線型導分条件が正則性に必要十分。
量子背景への特殊化と上の正規化から、低階数なら可換型または対称型に限られ、
未分類低階数領域はこの表示を持たないことを証明した。

[Clifford点の低階数近傍](sections/clifford_rank_slice.tex) も対称チャートで被覆した。
階数6以下の有限Frobenius候補表示空間はこの点で滑らかな69次元となり、
既知の二次障害Qは交代一次変形の4×4行列のPfaffianとして現れる。
69は表示空間の次元で、全Clifford成分や粗モジュライの次元ではない。

全反可換点でも対称チャートの局所被覆を証明し、可換多項式の軌道は低階数表示空間で開となる。
三つの行列階数で切った**基底付きテンソル表示空間**には、69次元と45次元の二つの既約成分を同定した。
他の成分や境界は未分類であり、全ASモジュライの成分とは主張しない。

### 次の研究

**最優先は、二窓の全次数AS条件を実効化し、対称閉包外を含む網羅的な分類木を構築すること。**
対称閉包内のAS対象は対称完全交叉チャートに尽きると証明した（`thm:symmetric-closure-as-exact`）。
階数5を仮定すると、各係数空間の三次小行列式に共通因子が生じる。動く平面の条件で五次元の行列空間を一つの正規形に絞ると、反対二位置で同じ評価平面が現れ、四辺の関係が消える指数成長商が得られる。従って低階数で未分類なのは階数6の残余である。階数6で四辺の点対応が有限なら点列の四移動は恒等になる。有限枝では長さ20の座標環と三つの座標四空間から二窓と添字保持同型を復元できる（[証明](sections/finite_point_chart.tex)）。二十点例の第三空間は符号倍で尽き、有限枝でのAS延長は一意である。一般の第一窓は四辺有限なAS対象に延びない（[復元と三次障害](sections/finite_window_rigidity.tex)）。この二十点を通る反転窓の成分は44次元で、基底変換による有理商は14次元となる。対称部分に横断する特異性は六変数二次錐で、幾何学的一般点の最小周期は2、第三階数はちょうど8となる。印付きの二つの三行対称行列系から一般窓を構成でき、14次元の有理商は単有理的である（[構成と階数](sections/marked_opposite_rank.tex)）。一般点の全次数のAS正則性も示した（[成分の局所方程式](sections/reversible_component.tex)）。二窓のExt反復は四項の微分と経路の乗法から計算でき、指定基底で各窓の先頭語の重なりが四つなら、有限Frobenius候補を全次数で認証できる。

滑らかな楕円曲線の二つの次数4の直線束と、その積の基点自由鉛筆から、14次元の商の稠密開部分を構成した。このチャートでは全ファイバーがAS正則で、異なる二束の場合は曲線・順序付き二束・鉛筆を復元できる。節点・尖点および同じ正の多重次数の環状境界でも全次数認証と同型判定が成立する。環状境界では異なる正の多重次数は指数成長商で排除できる。正次元の点対応を持つ一周期例の被約支持は楕円曲線を正規化とする曲線と四孤立点である。[証明と境界条件](sections/elliptic_reversible_pencil.tex)

正の多重次数を持つ被約二二次式完全交叉では、二束の多重次数の一致が鉛筆構成の標準AS正則性の必要十分条件になる。接点・三重点・四枝の楕円多重点でも、異なる二束の同型は曲線データから判定できる。[多重点の証明](sections/additive_pencil_boundary.tex)

二重二次曲線の鉛筆も全ファイバーがAS正則で、異なる二束の9次元の層では冪零方向まで同型を復元できる。反転テンソルの隣接階数が6、反対階数が8で二つの反対核が滑らかな楕円四次曲線なら、第一窓の一つの単純零点から楕円表示とAS正則性を認識できる。[二重曲線と逆判定](sections/ribbon_and_graph_recognition.tex)
非被約な二二次式完全交叉の全ての正次数鉛筆についてAS正則性を証明した。異なる二束では冪零方向を含む曲線の同型が添字保持同型の必要十分条件になる。滑らかな反対核を持つ有限窓では、第一窓の零点の台が11点以上なら楕円表示に入る。[非被約曲線の証明](sections/nonreduced_quartic_pencils.tex)・[小行列式による逆判定](sections/ribbon_and_graph_recognition.tex)

隣接階数6・反対階数8で両反対核が滑らかな四次曲線を切る反転窓では、純テンソル関係の不在が標準AS正則性の必要十分条件となり、全AS点が楕円鉛筆の表示を持つ。印付き零点の接空間次元0・1・2を全て扱った。さらに両対称テンソルから得る標準AS反転窓は全て楕円鉛筆の閉包に属し、反対階数は8以下である。階数9・10は反対核の因子型によらず排除された。階数8のAS点では片側の反対核が滑らかなら他側も同じ滑らかな曲線を切り出す。[変形と判定式](sections/marked_gram_smoothing.tex)

片側の極限モデルが整な節点・尖点四次曲線のときにも階数7のAS窓を排除した。非可逆な束の極限では他側はねじれ三次曲線と割線または接線になる。この候補の基点自由な全鉛筆で二次双対の有限Frobenius性と同型判定を得た。節点・尖点の各一例は六次まで標準次元で、七次で121となり標準値120を超える。節点型の一パラメータ族全体でも七次が121以上となる。さらに全ての非可逆鉛筆で、二変数二次式・三次式からの共通商が十七次で1152次元となり標準値1140を超える。整な極限曲線を保つAS点はその上の可逆な二束の表示に尽きる。環状・楕円多重点模型では零次数を許しても二束の多重次数一致が必要十分である。[全鉛筆の排除](sections/binary_product_obstruction.tex)・[正規化模型と三次計算](sections/integral_apolar_limits.tex)

両側が被約で一方が節点多角形となるAS極限では、他側も同じ多角形となり多重次数が一致する。異なる次数配分174通りは、退化で乗法の階数が落ちる場合も含めて排除した。[二側の特殊化](sections/polygon_apolar_limits.tex)

二側がともに被約四次完全交叉となる全てのAS特殊化で、二曲線は同じ曲線の等しい正の多重次数の鉛筆から得られる。従って階数8の被約反対核ではグラフと鉛筆で認識でき、階数7のAS極限では少なくとも片側の極限鉛筆が被約完全交叉でなくなる。[多重点グラフと統一定理](sections/reduced_additive_limits.tex)

反対核に重平面がある両対称窓では、軸の評価階数と三変数の三二次式の共通零点で全次数のAS正則性を判定できる。正則な反対階数は6または8で、8では両核が同じ重平面付き鉛筆となり、順序付き同型を係数で判定できる。[重平面のOre表示と因子対応](sections/square_apolar_ore.tex)

両対称AS反転窓の反対階数は **6または8** に限られる。6では完全対称な四二次式完全交叉の既存分類に帰着し、8では両核が四次完全交叉となる。異なる平面因子を指定し、横断する二軸では二つの三次元可換積を、共通根基が直線のときは二つの三次形式を比較して階数7を排除した。[因子の一意性と階数二分法](sections/distinct_factor_rigidity.tex)

**両対称反転窓の枝ではAS正則性と添字保持同型の判定が完結した。** 相異なる二平面の評価像が横断すれば三次元可換単位積の非結合性が必要十分で、共通根基が直線ならそこから得る二行二列関係の可逆性が必要十分である。二つの三生成Clifford代数の自由積を共通部分で貼り合わせ、全次数のTorを計算した。平方・被約反対核の結果と合わせて全入力を覆う。[因子対の正則性と分類手順](sections/transverse_factor_amalgam.tex)

一組だけの反対対称性を仮定しても、AS窓では四つの評価が全て階数3以上で、二つの窓は同じ点スキーム上の恒等四移動を与える。第三窓は第一窓の24元の線型系から復元できる。有限な場合は共通の最小解消から両奇位置の二乗像が決まり、それを含むVeronese三次元多様体の一意性で両対称性が従う。[片側対称の成長障害と相手窓](sections/one_opposite_symmetric_evaluations.tex)

相異なる二つのVeronese三次元多様体に二乗像が入る被約二十点を構成した。対応する片側対称窓は反対階数10、最小周期2であり、全次数で標準AS正則となる。四つの正規二次式と多項式環上の自由性を使い、この例と同じ仮定の窓を次数6までの検査で認証できる。[二十点の相手窓と全次数証明](sections/twenty_point_companion.tex)

三軸族では全係数でAS正則性と添字保持同型が決まった。正則条件は三つの平方の二次基本対称式が非零であること、同型条件は三つの平方の順序を忘れた射影組の一致である。反対階数10の同型類の像は二次元で、点スキームに直線が現れる特殊値を含む。[証明と判定](sections/cross_product_parameters.tex)

四元群で斉次な片側対称窓のうち、一つの文字次数が欠ける全入力で正則性を判定できた。反対階数10で点スキームにねじれ二直線を持つ例は三軸族の外にある。次数次元が3,1,1,1の場合は全対称性が必要となる。[証明と判定](sections/klein_companion_windows.tex)

残る次数次元2,2,1,1では、二行列の非対角成分が一致する。第二の行列式解から、同型類の像が四次元の新しい族が得られた。幾何学的一般点は反対階数10の標準AS正則窓で、被約二十点のどの三点も共線でない。[証明と判定](sections/klein_two_doublets.tex)

この四次元族は一般に最小周期4で、三つの二十点配置は互いに射影同値でない。各配置の射影自己同型群は四元群に一致し、明示した開集合では四つの有理不変量と十六通りの置換で印なし同型が決まる。[証明と判定](sections/klein_generic_isomorphisms.tex)

次数次元2,2,1,1の全正則窓を、両対称枝と四次元族の閉包で被覆した。交差行列が退化すると正則性には両対称性が必要となる。薄い次数に零がある場合と、補完の解が無限にある例外値も含む。さらに具体的な有理係数の最小周期4の例で四つの中心四次式を構成し、全次数の正則性を証明した。[証明と判定](sections/klein_doublet_boundary.tex)

中心四次式を独立な二変数の有理関数体に持ち上げた。この断面の一つの開集合では、全ての閉点が四つの中心四次式により標準AS正則となる。第二補完の全閉包でこの自由表示を認識する有限判定も得た。[証明と判定](sections/klein_quartic_center_lift.tex)

四変数全体で五つの中心四次式が持ち上がった。一つの稠密開集合の全閉点が正則であり、閉包内の全ての標準ASファイバーには共通の多項式恒等式が成立する。[証明と判定](sections/klein_generic_central_parameters.tex)

偶位置の対称性を外して、反対階数16・最小周期4の五次元族を構成した。最後の変数は楕円曲線上を動く。有理閉点で中心四次式による全次数の正則性を証明し、四つの二十点配置が全て互いに射影非同値であることも検証した。[証明と判定](sections/klein_elliptic_completion.tex)

楕円曲線の根と別の組合せは四つの倍率でのみ生じる。それらは両対称窓またはそのコサイクル捻りで判定できる。明示的楕円曲線上でも五つの中心四次式が延び、一般の印なし同型は32通りの置換で判定できる。[証明と判定](sections/klein_second_pairings.tex)

四元群を仮定しない片側対称窓に、一般の射影自己同型が自明な十四次元族を得た。階数10と16の両有理点で対称性も外した接空間は十五次元となり、全二次障害が消える。[証明と判定](sections/unrestricted_doublet_deformations.tex)

二十点から階数方程式の間の独立な二十の恒等式を構成し、十五の接方向を全て実現する成分を証明した。一般点は反対階数16・最小周期4で自己同型は自明である。片側対称点の近くでは階数10の滑らかな因子の外は全て階数16となる。[証明と判定](sections/twenty_conormal_component.tex)

行列式二次式が両方とも消える八種類の基底を決定した。各ファイバーは一つの双一次式で定まり、全て楕円曲線族の閉包に入る。これで薄い成分と主ブロックが非退化な範囲の階数条件を覆った。[証明と判定](sections/klein_zero_quadratics.tex)

十五次元成分では一位置の二十点配置だけで一般の窓が有限候補に絞られ、この配置自身の一般射影安定群も自明である。逆向きの窓も同じ配置を持つため復元次数は偶数である。階数10の点で偶位置の二つの平方根が重なる局所長さは2となる。[証明と判定](sections/single_configuration_projection.tex)

一つの二十点配置からの一般復元次数はちょうど2で、二つの向きが全候補である。二十個の三重点を持つ四次鉛筆から、退化解を含めて全平方根を決定した。[証明と判定](sections/configuration_quartic_pencil.tex)

偶奇それぞれの二つのVeronese像の交叉が二十点そのものであることを示した。鉛筆を一つ残すと四次の四元群拡大、Gale双対な二鉛筆の順序を忘れると八次の二面体拡大となる。奇対称因子で二つの向きは合流する。[証明と判定](sections/veronese_pencil_quotient.tex)

四次の導分を計算し、十五次元成分の全方向に五つの中心四次式が延びることを証明した。一つの開集合の全幾何点が正則となり、閉包の全標準AS点は共通の多項式恒等式を持つ。[証明と判定](sections/quartic_euler_lifting.tex)

全中心は五つの四次中心元で生成され、一つの四次関係を持つ正規な整超曲面である。十五次元成分の開集合でPI次数は16、各頂点の商は次数4の可除代数となる。[証明と判定](sections/quartic_center_order.tex)

片側対称因子では全中心は二次超曲面の二重被覆で、四次分岐式を持つ。四次中心元の積写像の行列式がこの因子で一重に消える。[証明と判定](sections/symmetric_center_cover.tex)

十五次元成分の一般窓で、射影中心には最小単純加群に対応する相異なる二十の特異点がある。各指標上の単純加群は二種類であり、滑らかな軌跡とAzumaya軌跡は一致する。[証明と判定](sections/minimum_central_characters.tex)

二十個の特異点は全て通常二重点であり、各点の単純加群の次元は4と12である。完備化した折り畳み環を二行二列行列の余核から復元した。[証明と判定](sections/twenty_central_nodes.tex)

四つの中心四次式で有限次元商が得られる全ての有限Frobenius候補について、経路長16までの必要十分なAS正則性判定を与えた。[証明と判定](sections/central_window_cutoff.tex)

二つの行列式二次式が恒等的に消える境界で、八表示の全点をAS分類から除外した。長さ6の一般経路次元92は、二直線に分かれる特殊ファイバーでも下界になる。[証明と判定](sections/zero_pairing_defect.tex)

行列式二次式の片方だけが消える点は四つの例外倍率に限られ、全32表示が経路長7の次元超過で非正則となる。[証明と判定](sections/single_zero_pairings.tex)

真の重根の判別面は十六の通常二重点を持つKummer四次曲面である。この判別族上に反対階数16の正則な有理点と、全幾何点が正則な開集合を得た。[証明と判定](sections/kummer_ramification.tex)

重根での二行列の交換はKummer曲面の極写像である。十六の二重接触平面を明示し、その上の許容窓を全て非正則と判定した。[証明と判定](sections/kummer_polar_completion.tex)

符号付き判別被覆の可逆部分で、射影化した対根のグラフが八節点を小解消する。各例外直線の第二行列像は対応する接触二次曲線である。[証明と判定](sections/paired_root_resolution.tex)

二行列交換は全ての通常倍率で解消Kummer曲面上に固定点を持たず、印を固定した二位置交換の商はEnriques曲面となる。[証明と判定](sections/kummer_exchange_quotient.tex)

四つの例外倍率では判別面は滑らかな二次曲面の二重化となり、可逆窓は二つの式がともに真の二次式であるとき、かつそのときに限り標準AS正則となる。[証明と判定](sections/exceptional_ramification.tex)

例外倍率の四つの交換表示全体で非正則集合は八直線に一致し、単根窓は全て正則となる。[証明と判定](sections/exceptional_exchange_sheets.tex)

全ての両対称反転窓について、二始点の長さ3と21の経路次元だけで標準Koszul左右AS正則性を判定できる。[証明と判定](sections/bisymmetric_finite_path_test.tex)

純関係のない六平面で、二変数形式の商を階数二の有理正規五次曲線から復元できる。四辺の曲線が頂点ごとに合えば、一周の座標変換が恒等でなくても指数成長を検出する。[証明と判定](sections/binary_window_transfer.tex)

一辺の二変数商を持つ六平面は余次元13の既約軌跡をなし、一般の核は一意である。階数二曲線の残りは次数15、正規化の種数10となる。[証明と判定](sections/binary_incidence_sextic.tex)

片側対称な第二補完の階数二曲線を三つの超平面で切り、良い還元で成分次数を4と16に決定した。この辺は標数零の開近傍で両方向とも二変数商を持たない。[証明と判定](sections/companion_rank_two_sections.tex)

階数16の単根・重根の両近傍で四辺全ての階数二曲線は次数4と16に分かれる。四次成分から二文字ブロックを、第一辺から印なし不変量を復元できる。四辺に整合する順序なし文字直線の候補は高々三通りである。[証明と判定](sections/rank_two_block_recovery.tex)

薄い二色と四次成分の分岐因子で文字直線の分解が一意に決まる。二十点配置を計算せず、32置換と係数の積で印なし添字保持同型を判定できる。[証明と判定](sections/unique_character_frames.tex)

四位置の同時文字交換が自己同型となるのは二直線で、平方の計算により別の文字分解が生じ得るのは正の直線に限られる。その上で二つの正則表示を持つ有理例を得た。[証明と判定](sections/flip_character_locus.tex)

二つの固定直線上の楕円曲線の分岐八点を同時に正則と認証し、単根でも長さ4で非正則な十二表示を得た。[証明と判定](sections/fixed_line_elliptic_curves.tex)

両楕円曲線上の二位置移動は位数2の平行移動で、明示した商曲線に長さ4の障害が降りる。[証明と判定](sections/fixed_exchange_isogeny.tex)

倍率(3,2)の両固定曲線で長さ4の障害は全幾何点を通じて十六点で尽きる。新しい四点は全通常倍率で反対階数16の非正則表示に延びる。[証明と判定](sections/fixed_geometric_quartic_locus.tex)

全通常倍率で二つの固定曲線上の長さ4の障害を決定した。一般の十六点は特別な倍率で十二点に合流する。[証明と判定](sections/uniform_fixed_quartic_support.tex)

倍率(3,2)で長さ4を通過した全幾何点は長さ5も通過する。二つの八次数体に残る候補を別の小行列式で除いた。[証明と判定](sections/fixed_quintic_locus.tex)

負の曲線では全通常倍率、正の曲線では(3,2)を含む開集合で、長さ5に初めて現れる障害がないことを証明した。[証明と判定](sections/varying_fixed_quintics.tex)

両固定曲線では全通常倍率で、長さ5の障害点が長さ4の障害点と一致することを証明した。[証明と判定](sections/uniform_positive_quintics.tex)

倍率(3,2)の正の曲線で、端点二点、相方四点、八次閉点二組の全共役十六点について全次数のAS正則性を証明した。[証明と判定](sections/fixed_residual_regularity.tex)

反対階数14が現れる倍率と四点を全通常倍率で決定した。二つの成分の完備化はともに j=1728 の楕円曲線で、具体的な倍率で階数14・15の全点の全次数正則性を証明した。[証明と判定](sections/fixed_opposite_rank_drop.tex)

階数14の二位置移動を二重被覆の共役と同定し、全域の表示から中心の計算を一変数に降ろした。101で有理な全192剰余点に還元される標数零の窓はAS正則である。[証明と判定](sections/rank_fourteen_descent.tex)

階数14の全通常点で四次中心の文字次元を (3,2,0,0) と決定し、中心小行列式の例外に属する全280幾何点の全次数正則性を証明した。[証明と判定](sections/rank_fourteen_center.tex)

階数14の全通常点で四中心商の長さ6までの次元と、元の道空間の長さ6の次元84を確定した。例外22閉点の上にある4768幾何点は全次数でAS正則となる。[証明と判定](sections/rank_fourteen_central_cut.tex)

標数101の全通常幾何剰余点への持上げと、倍率座標が1の冪根である全通常窓の全次数正則性を証明した。[証明と判定](sections/rank_fourteen_geometric_reduction.tex)

斉次表示の全48境界パラメータを判定し、32点の全次数正則性と16点の非標準性を決めた。通常域で未判定の閉点は101で全ての還元が六つの値に入る場合に限られる。[証明と判定](sections/rank_fourteen_boundary.tex)

四つの分岐点で生成元を重み付き変換し、そこに101進還元される通常窓の全次数正則性を証明した。階数14の未判定閉点は、全ての還元が19または80となる場合だけである。[証明と判定](sections/rank_fourteen_ramified_disks.tex)

反対階数14の二つの楕円成分の全通常点で標準Koszul左右AS正則性を証明した。最後の二剰余類では局所変数で関係式を割り、長さ13までの四中心商を全16組で認証した。[証明と判定](sections/rank_fourteen_all_ordinary.tex)

正の固定曲線の反対階数15を二つの四重被覆で表示し、全通常倍率での根数と分岐を決定した。分岐曲線は j=35152/9 の楕円曲線で、その一般窓と具体的な16点が全次数で正則となる。[証明と判定](sections/rank_fifteen_covers.tex)

階数15の二分岐曲線の全通常点で四次中心と長さ6までの中心商を決定した。消去に残る5336点は全次数で正則。[証明と判定](sections/rank_fifteen_branch_center.tex)

階数15の二分岐曲線で、101で境界または有理点に還元される全通常点の正則性を証明した。第一座標が有理数の全点を含む。[証明と判定](sections/rank_fifteen_boundary_disks.tex)

階数15の二分岐曲線の全通常点で正則性を証明した。残る階数15面を種数5の曲線族に分け、二つの種数3の特殊値を八次超楕円式で表示した。[証明と判定](sections/rank_fifteen_all_ramification.tex)

階数15面の一般水準と種数3の二特殊水準で、無限遠八点の正則な重み付き極限を求めた。種数3で109におけるKの極を持つ全通常窓が正則となる。[証明と判定](sections/rank_fifteen_pencil_infinity.tex)

種数3の四曲線の全通常点で四次中心を決定した。消去に残る7104点と、109で通常な有理剰余または二つの追加境界に還元する全通常窓が正則である。[証明と判定](sections/rank_fifteen_octic_center.tex)

種数3の四曲線で109の有限境界 K=0,±2 に還元する全通常点を認証した。残る K=±1 では二次関係の極限の経路次元を数体上で決定した。[証明と判定](sections/rank_fifteen_octic_boundary.tex)

種数5の各水準の二位置移動商を平面四次式で与えた。一般水準の有限非分岐境界24点では、その近傍の通常窓が全次数で正則である。[証明と判定](sections/rank_fifteen_pencil_unramified.tex)

種数5の一般水準で K=±s±1 上の16境界点も認証した。通常近傍を認証した境界は56点中48点で、残る八点は K=±1 上にある。[証明と判定](sections/rank_fifteen_pencil_ramified_ends.tex)

階数15の二つの種数3水準の両表示について全通常点の標準Koszul左右AS正則性を証明した。各始点で四次中心の多項式環上の自由階数は256である。[証明と判定](sections/rank_fifteen_octic_all.tex)

種数5の水準 ξ=9/16 で、101における56境界点のどれに還元する通常窓も、両表示で正則となる。[証明と判定](sections/rank_fifteen_pencil_trace_disks.tex)

階数15で水準を支配する非正則曲線の特殊化先を、各種数3水準の K=±1 上の八点に限定した。K=0 の節点にも曲面の二変数で正則な関係格子を構成した。[証明と判定](sections/rank_fifteen_horizontal_specialization.tex)

階数15の残る八点で水準を変動させたまま、長さ5までの経路、四次中心と四式による商を認証した。[証明と判定](sections/rank_fifteen_relative_trace_centers.tex)

階数15の両表示で、有限個の代数的水準を除き、種数5曲線の全通常点の標準Koszul左右AS正則性を示した。[証明と判定](sections/rank_fifteen_generic_levels.tex)

階数15の両曲面で、標準Koszul左右AS正則性が失敗し得る通常窓を有限個の代数点に限定した。全ての種数5水準の一般点を無限遠で認証した。[証明と判定](sections/rank_fifteen_vertical_curves.tex)

階数15の全通常窓で二次双対の完全対合と次元列14641を示した。孤立候補でもこの条件は成立する。[証明と判定](sections/rank_fifteen_quadratic_dual.tex)

階数15の両表示で添字保持同型による像が各々二次元であることを証明した。二つの独立な不変量と、行列積を用いる同型の必要十分条件を得た。[証明と判定](sections/rank_fifteen_adjoint_invariants.tex)

同じパラメータで階数15の二表示は反対代数となり、六次の階数判定は四行列に減る。一般窓の最小周期は2で、添字を固定すれば二表示を三次形式で区別できる。[証明と判定](sections/rank_fifteen_word_symmetries.tex)

階数15の全座標端で、両表示・全始点の六次経路次元は84、四次中心の文字次元は (3,2,0,0) である。二つの曲線の残余18因子上では全次数の正則性も得た。[証明と判定](sections/rank_fifteen_coordinate_ends.tex)

倍率を反転すると階数15の二表示は添字保持同型類の同じ集合を与える。各座標端成分内の一般同型類は相互反転する二点に限られ、楕円曲線と有理曲線の像は異なる。[証明と判定](sections/rank_fifteen_endpoint_orbits.tex)

楕円座標端 n=±1 の全通常窓で標準Koszul左右AS正則性を証明した。座標端の三つの添字保持像は相異なり、正規化は j=35152/9 の楕円曲線一本と有理曲線二本である。[証明と判定](sections/rank_fifteen_coordinate_regularity.tex)

座標端の三曲線上で、同時反転以外の同型が残り得る水準を明示した。楕円・有理間の重複には両側でそれぞれ特定の三次式の消滅が必要である。[証明と判定](sections/rank_fifteen_coordinate_collisions.tex)

階数15の両表示で、全座標端の全通常点が標準Koszul左右AS正則である。有理成分の残余12軌道も解決し、正則性の失敗候補は四座標がすべて非零の内部に限られる。[証明と判定](sections/rank_fifteen_all_coordinate_ends.tex)

階数15の一般窓の添字保持同型を48置換の有限判定に還元し、射影安定化群の位数16を決定した。倍率が101で (4,9) に還元される全窓は正則である。[証明と判定](sections/rank_fifteen_pauli_comparison.tex)

階数15の各通常倍率の添字保持同型類は非分岐なら二つ、分岐上なら一つである。非分岐全点の最小周期は2で、座標端の二水準の同型も確定した。[証明と判定](sections/rank_fifteen_root_exchange.tex)

分岐曲線上でも周期1を排除し、正の固定曲線の通常な階数15窓の最小周期が全点で2であることを示した。[証明と判定](sections/rank_fifteen_branch_period.tex)

剰余標数7,11,13,17,19の全480有理通常倍率に還元される階数15窓の全次数正則性を証明した。[証明と判定](sections/rank_fifteen_small_residues.tex)

階数15の一般添字保持同型類を与える第一表示は12組の倍率上の24点と確定した。二積の生成次元16の全点で同型な全パラメータを192変換から復元できる。[証明と判定](sections/rank_fifteen_parameter_orbits.tex)

倍率平面の三本の新しい内部曲線の全通常点で全次数正則性を証明した。変換先が非通常となる有限例外も直接のテンソル検査で扱った。[証明と判定](sections/rank_fifteen_regular_orbits.tex)

階数15の一般添字保持同型商は有理曲面で、二つの明示量が一般類を決定する。全通常点の左右正則性は二変数 (x,y) だけに依存する。[証明と判定](sections/rank_fifteen_rational_quotient.tex)

正の固定曲線の通常な階数15窓すべてで、三量 (x,p,q) の一致が添字保持同型の必要十分条件となる。特殊倍率・座標端も含む。[証明と判定](sections/rank_fifteen_global_invariants.tex)

階数15の正則性平面の剰余標数3,5,7の全68有理点（yが単元）で、全ての標数零の持上げを全次数正則と認証した。[証明と判定](sections/rank_fifteen_quotient_residues.tex)

通常倍率の像が商曲面から欠く四点を決定し、全四点で正則窓を認証した。二点は最小周期1、他の二点は周期2である。[証明と判定](sections/rank_fifteen_formal_boundary.tex)

剰余標数3で x が零に落ち、y が単元となる直線の全幾何点について四中心切断を検証し、標数零の持上げの正則性を証明した。[証明と判定](sections/rank_fifteen_residue_line.tex)

Pauli射影境界の32点が一つの周期1の量子多項式環を与えることを証明し、その極限で階数15正則性平面の付値領域全体を全奇標数で認証した。[証明と判定](sections/rank_twelve_quantum_limit.tex)

Pauli射影閉包の全境界で正則性を判定した。階数14・12は正則、階数9・4は有限障害で排除され、付値の等号でも剰余標数零・101で正則となる。[証明と判定](sections/pauli_projective_boundary.tex)

三次曲線 p=0 の全点で正則性を証明した。二点は周期1、残りは周期2である。境界で増える経路を3進の飽和計算で除き、曲線近傍にも判定を拡張した。[証明と判定](sections/rank_fifteen_cubic_zero_disk.tex)

x=0 の全窓を正則と証明し、添字保持同型類を穴あき二次曲線で表示した。全点の最小周期は2である。[証明と判定](sections/rank_fifteen_x_zero.tex)

3で跡と局所変数がともに消える領域で、全始点の第六経路次元84と四次中心の文字別次元を証明した。[証明と判定](sections/rank_fifteen_variable_zero_disk.tex)

3進の可変跡円板で選んだ四中心式に一次の関係がなく、全始点の四中心商の長さ5の次元が40であると証明した。[証明と判定](sections/rank_fifteen_first_center_products.tex)

3進の可変跡円板で四中心商の長さ6は文字0,1で各12次元、残る二文字で各10から12次元となる。[証明と判定](sections/rank_fifteen_saturated_center_products.tex)

傾き59319/10000の全点の正則性を証明し、周期2の添字保持類を種数2曲線の四点補集合で分類した。[証明と判定](sections/rank_fifteen_gaussian_pencil.tex)

可変跡円板で残る両文字の小行列式を確定し、全始点の六次中心商の次元44を証明した。[証明と判定](sections/rank_fifteen_sixth_center_minors.tex)

整な傾きで六次経路84、六次中心商44となる範囲を確定した。剰余標数零と全ての奇素数で剰余傾きに例外はない。[証明と判定](sections/rank_fifteen_integral_slope.tex)

3進の可変跡環帯で左右の積を照合し、全次数のKoszul・左右AS正則性と中心上の自由階数256を認証した。[証明と判定](sections/rank_fifteen_overlap_annulus.tex)

水平切片の添字保持曲線、移動商、二つの節点退化と周期層を決定し、切片y=3の明示的な正則領域を得た。[証明と判定](sections/rank_fifteen_horizontal_fibers.tex)

可変跡円板全体の左右AS正則性と自由階数256を証明し、3の上で正付値の傾きを持つ全切片を分類した。[証明と判定](sections/rank_fifteen_three_residue_classes.tex)

整な傾きの3進環帯を全剰余類で認証し、両座標が素体の単元へ還元する四領域も全次数で決定した。[証明と判定](sections/rank_fifteen_integral_annulus.tex)

整な傾きで正付値の全領域を深さ制限なく認証し、有理座標の整傾き全点と傾きが極になる別の環帯も決定した。[証明と判定](sections/rank_fifteen_integral_depth.tex)

極をもつ傾きの八重被覆で、単元の剰余座標 $x$ が1以外の環帯全体を認証した。[証明と判定](sections/rank_fifteen_pole_fibers.tex)

剰余座標が1の分岐方向でも、正規化した接線が衝突する値を除き、六桁で正則性を証明した。[証明と判定](sections/rank_fifteen_ramified_pole.tex)

単元座標の非分岐極と非衝突の接方向で、正付値の深さの制限を外した。[証明と判定](sections/rank_fifteen_unit_pole_depth.tex)

外側の極を二変数で表示し、4ν(x−1)=ν(y)>0 の全点を標数零で決定した。[証明と判定](sections/rank_fifteen_outer_pole.tex)

二変数のSmith指数を両辺と対角線で比較し、外側の極の全ての傾きを決定した。[証明と判定](sections/rank_fifteen_outer_pole_fan.tex)

内側衝突の深い範囲と中間の等号を認証。3進単元の有理数 x が1 mod 9以外なら全ての有理水準が決まった。[証明と判定](sections/rank_fifteen_inner_collision.tex)

内側衝突の残る三つの錐を認証し、x が3進単元なら全ての非零有理水準で正則性を確定した。さらに x が1に還元される全縦切片を、非零の全幾何水準で分類した。[証明と判定](sections/rank_fifteen_inner_collision_fan.tex)

非単元座標の極のうち、0 < 5ν(x) ≤ ν(y) の全域で正則性を確定した。[証明と判定](sections/rank_fifteen_nonunit_pole.tex)

x が正付値の全縦切片を全幾何水準で分類し、x が3整な全有理通常点で正則性を決定した。[証明と判定](sections/rank_fifteen_nonunit_pole_strip.tex)

3整な有理数 x で、3上に剰余次数4以下の素点を持つ任意の数体の非零水準 y を認証した。[証明と判定](sections/rank_fifteen_negative_unit_levels.tex)

x が負付値の等号境界を全幾何傾きで認証し、対角線近傍の正則性と残る錐の細分条件を得た。[証明と判定](sections/rank_fifteen_negative_x_boundary.tex)

負付値の二曲線の近傍を認証し、二変数の指数比較に残る四次中心の差を特定した。[証明と判定](sections/rank_fifteen_negative_x_arcs.tex)

円錐曲線の極で平方座標を取り、10ν(y) ≤ 19ν(x) < 0 の閉領域を全て認証した。[証明と判定](sections/rank_fifteen_negative_x_fold.tex)

負付値 x の全縦切片を全水準で分類し、全有理座標の通常な階数15窓を認証した。[証明と判定](sections/rank_fifteen_negative_x_fan.tex)

3上の剰余次数が2以下の素点を持つ数体で全座標を認証した。一般単元の四窓は種数13の関数体に降下する。[証明と判定](sections/rank_fifteen_unit_quartic_towers.tex)

全ての有理数 x の通常な階数15窓を任意の非零の幾何水準で認証した。剰余が −1 の縦切片は全て添字保持同型まで分類した。[証明と判定](sections/rank_fifteen_genus_two_unit_line.tex)

横座標が3上に剰余次数2以下の素点を持つ数体に属せば、任意の非零幾何水準で通常窓を認証した。[証明と判定](sections/rank_fifteen_quadratic_verticals.tex)

全ての非零有理倍率 λ で、放物線 y=λ(1−x)^2 の全通常窓を任意の非零幾何水準で認証した。[証明と判定](sections/rank_fifteen_prime_parabolic_pencils.tex)

二次剰余の全倍率で放物線上の全通常窓を認証した。語の偶奇式により、一般の経路と文字0の中心商行列を二変数有理関数体へ降下した。[証明と判定](sections/rank_fifteen_quadratic_parabolic_pencils.tex)

標数3の非分岐単元曲面全体で、五次までの全始点の経路と四次中心の同一の五本の基底を認証した。[証明と判定](sections/rank_fifteen_bivariate_word_paths.tex)

非分岐単元域全体で任意の独立な四中心の五次商を認証し、全幾何単元対に六次44次元の四中心商を与えた。[証明と判定](sections/rank_fifteen_bivariate_sixth_cut.tex)

六次の固定四列の例外点対を解決し、九次で23次元になる表示と、一列の交換で十三次まで回復することを証明した。[証明と判定](sections/rank_fifteen_fixed_center_obstruction.tex)

全幾何単元対で七次四中心商が四始点とも40次元となること、および十三の剰余軌道で十三次までの証明書を得た。[証明と判定](sections/rank_fifteen_bivariate_seventh_cut.tex)

標数零の全通常階数15窓で七次までの中心単項式基底を証明した。八次から十一次の商次元には下界31,20,10,4があり、各次数には上界の証明だけで足りる。[証明と判定](sections/rank_fifteen_truncated_central_freeness.tex)

符号付き六表示で始点を保ち、全通常階数15窓の四始点で八次の未切断165次元・四中心商31次元を証明した。[証明と判定](sections/rank_fifteen_oriented_eighth_cut.tex)

九次では固定四列の商が両符号の全十二表示・四始点で23次元となる。中心四平面の交換が必要な点を特定した。[証明と判定](sections/rank_fifteen_ninth_orbit_obstruction.tex)
同じ剰余表示では、三中心を共有する鉛筆の例外値はちょうど三つである。五つの中心四次式の上で有限な標準AS窓では、悪い四平面を一つの隅の次数64の行列式で記述できる。正の固定曲線の通常階数15窓は全て四次中心上で有限であり、四中心商の長さ13までの階数条件がAS正則性の必要十分条件となる。[証明と適用範囲](sections/three_center_pencil.tex)

縦切片 $x=1$ の正付値領域を全て認証し、全有理水準で添字保持同型類と最小周期を決定した。[証明と判定](sections/rank_fifteen_vertical_one.tex)

座標1の単元剰余を種数1の48重被覆で表し、縦切片の三次以下の全代数的水準を決定した。[証明と判定](sections/rank_fifteen_vertical_unit_levels.tex)

縦切片 $x=1$ の全水準を分類した。種数1の被覆から残る四つの剰余軌道を全て検証した。[証明と判定](sections/rank_fifteen_vertical_all_levels.tex)

二窓から作る56次元代数で、低次Extの消滅を逆Serre反復が保つ構造条件を探す。
「全次数を証明してAS認証」「有限障害で排除」「未決定」を分け、有限検査の成功を全AS認証としない。
正規形の枝、特殊値、枝間同型の条件を証明し、未分類の残余を明示する。
既存の九変数族・量子成分・Clifford点・境界族は、全分類を目指す方法の回帰試験として使う。
九変数族の最大性・完全族内同型やClifford点の高次障害は有用な補助課題だが、
全分類の前提条件にはしない。詳細は研究方針と分類アトラス設計を参照する。
機械可読な巨大分類台帳と自動生成器は、新しい被覆分岐と実際の規模に合わせ最小構成で導入する。

## 再現

Python 3 と SymPy 1.14.0。LuaLaTeX・latexmk・luatexja・Harano Ajiフォントを用いる。
楕円曲線の関数体と二十点配置の行列計算には python-flint 0.9.0、中心の法二乗計算には NumPy 2.3.5 も使う。すべてリポジトリのルートから実行する。

~~~sh
python3 -m pip install -r requirements-symbolic.txt
python3 scripts/check_project.py
python3 computations/verify_exterior_trace.py
python3 computations/verify_period_witness.py
python3 computations/verify_quartic_obstruction.py
python3 computations/verify_finite_square_chart.py
python3 computations/verify_tangent_pencil_window.py
python3 computations/verify_finite_window_rigidity.py
python3 computations/verify_reversible_component.py
python3 computations/verify_marked_opposite_rank.py
python3 computations/verify_elliptic_reversible_pencil.py
python3 computations/verify_additive_pencil.py
python3 computations/verify_ribbon_trace.py
python3 computations/verify_nonreduced_polarizations.py
python3 computations/verify_smooth_apolar_defect.py
python3 computations/verify_conductor_pencil.py
python3 computations/verify_binary_product_obstruction.py
python3 computations/verify_polygon_multiplication.py
python3 computations/verify_additive_graph_limits.py
python3 computations/verify_square_apolar_ore.py
python3 computations/verify_transverse_apolar_factors.py
python3 computations/verify_factor_regular_windows.py
python3 computations/verify_single_opposite_block.py
python3 computations/verify_twenty_point_companion.py
python3 computations/verify_cross_product_parameters.py
python3 computations/verify_klein_companion_windows.py
python3 computations/verify_klein_two_doublets.py
python3 computations/verify_klein_generic_isomorphisms.py
python3 computations/verify_klein_doublet_boundary.py
python3 computations/verify_klein_period_four_center.py
python3 computations/verify_klein_center_lift.py
python3 computations/verify_klein_generic_quartics.py
python3 computations/verify_klein_elliptic_completion.py
python3 computations/verify_klein_exceptional_pairings.py
python3 computations/verify_unrestricted_doublet_deformations.py
python3 computations/verify_twenty_conormal_identities.py
python3 computations/verify_klein_zero_quadratics.py
python3 computations/verify_single_configuration_projection.py
python3 computations/verify_configuration_quartic_pencil.py
python3 computations/verify_configuration_intersections.py
python3 computations/verify_quartic_euler_derivations.py
python3 computations/verify_quartic_center_order.py
python3 computations/verify_symmetric_center.py
python3 computations/verify_minimum_central_characters.py
python3 computations/verify_central_nodes.py
python3 computations/verify_central_window_cutoff.py
python3 computations/verify_zero_pairing_defect.py
python3 computations/verify_all_zero_pairing_defects.py
python3 computations/verify_single_zero_pairings.py
python3 computations/verify_kummer_discriminant.py
python3 computations/verify_kummer_polar.py
python3 computations/verify_paired_root_resolution.py
python3 computations/verify_kummer_exchange_quotient.py
python3 computations/verify_exceptional_ramification.py
python3 computations/verify_exceptional_exchange_sheets.py
python3 computations/verify_bisymmetric_finite_test.py
python3 computations/verify_binary_window_transfer.py
python3 computations/verify_binary_incidence_sextic.py
python3 computations/verify_companion_rank_two_sections.py
python3 computations/verify_rank_two_block_recovery.py
python3 computations/verify_unique_character_frames.py
python3 computations/verify_flip_character_locus.py
python3 computations/verify_fixed_line_elliptic_curves.py
python3 computations/verify_fixed_exchange_isogeny.py
python3 computations/verify_fixed_geometric_quartic_locus.py
python3 computations/verify_uniform_fixed_quartic_support.py
python3 computations/verify_fixed_quintic_locus.py
python3 computations/verify_varying_fixed_quintics.py
python3 computations/verify_uniform_positive_quintics.py
python3 computations/verify_fixed_residual_regularity.py
python3 computations/verify_fixed_opposite_rank_drop.py
python3 computations/verify_rank_fourteen_descent.py
python3 computations/verify_rank_fourteen_center.py
python3 computations/verify_rank_fourteen_central_cut.py
python3 computations/verify_rank_fourteen_geometric_reduction.py  # Singular 4.4 (PATH または SINGULAR_BIN)
python3 computations/verify_rank_fourteen_cyclotomic.py  # Singular 4.4
python3 computations/verify_rank_fourteen_boundary.py
python3 computations/verify_rank_fourteen_ramified_disks.py
python3 computations/verify_rank_fourteen_all_ordinary.py
python3 computations/verify_rank_fifteen_covers.py
python3 computations/verify_rank_fifteen_branch_center.py
python3 computations/verify_rank_fifteen_boundary_disks.py
python3 computations/verify_rank_fifteen_all_ramification.py
python3 computations/verify_rank_fifteen_pencil.py
python3 computations/verify_rank_fifteen_pencil_infinity.py
python3 computations/verify_rank_fifteen_octic_center.py
python3 computations/verify_rank_fifteen_octic_residues.py
python3 computations/verify_rank_fifteen_octic_boundary.py
python3 computations/verify_rank_fifteen_octic_path_jump.py
python3 computations/verify_rank_fifteen_pencil_unramified.py
python3 computations/verify_rank_fifteen_pencil_ramified_ends.py
python3 computations/verify_rank_fifteen_octic_ordinary.py
python3 computations/verify_rank_fifteen_octic_trace_disks.py
python3 computations/verify_rank_fifteen_pencil_trace_disks.py
python3 computations/verify_rank_fifteen_horizontal_specialization.py
python3 computations/verify_rank_fifteen_relative_trace_centers.py
python3 computations/verify_rank_fifteen_generic_levels.py
python3 computations/verify_rank_fifteen_vertical_curves.py
python3 computations/verify_rank_fifteen_quadratic_dual.py
python3 computations/verify_rank_fifteen_adjoint_invariants.py
python3 computations/verify_rank_fifteen_word_symmetries.py
python3 computations/verify_rank_fifteen_coordinate_ends.py
python3 computations/verify_rank_fifteen_endpoint_orbits.py
python3 computations/verify_rank_fifteen_coordinate_regularity.py
python3 computations/verify_rank_fifteen_pauli_comparison.py
python3 computations/verify_rank_fifteen_root_exchange.py
python3 computations/verify_rank_fifteen_parameter_orbits.py
python3 computations/verify_rank_fifteen_rational_quotient.py
python3 computations/verify_rank_fifteen_global_invariants.py
python3 computations/verify_rank_fifteen_formal_boundary.py
python3 computations/verify_rank_fifteen_residue_line.py
python3 computations/verify_rank_twelve_quantum_limit.py
python3 computations/verify_pauli_projective_boundary.py
python3 computations/verify_rank_fifteen_cubic_zero_disk.py
python3 computations/verify_rank_fifteen_x_zero.py
python3 computations/verify_rank_fifteen_variable_zero_disk.py
python3 computations/verify_rank_fifteen_first_center_products.py
python3 computations/verify_rank_fifteen_saturated_center_products.py
python3 computations/verify_rank_fifteen_sixth_center_minors.py
python3 computations/verify_rank_fifteen_integral_slope.py
python3 computations/verify_rank_fifteen_overlap_annulus.py
python3 computations/verify_rank_fifteen_horizontal_fibers.py
python3 computations/verify_rank_fifteen_overlap_bound.py
python3 computations/verify_rank_fifteen_integral_annulus.py
python3 computations/verify_rank_fifteen_integral_depth.py
python3 computations/verify_rank_fifteen_smith_slopes.py
python3 computations/verify_rank_fifteen_pole_slope_annulus.py
python3 computations/verify_rank_fifteen_pole_fibers.py
python3 computations/verify_rank_fifteen_ramified_pole.py
python3 computations/verify_rank_fifteen_unit_pole_depth.py
python3 computations/verify_rank_fifteen_outer_pole.py
python3 computations/verify_rank_fifteen_outer_pole_fan.py
python3 computations/verify_rank_fifteen_inner_collision.py
python3 computations/verify_rank_fifteen_inner_collision_fan.py
python3 computations/verify_rank_fifteen_nonunit_pole.py
python3 computations/verify_rank_fifteen_nonunit_pole_strip.py
python3 computations/verify_rank_fifteen_negative_unit_levels.py
python3 computations/verify_rank_fifteen_negative_x_boundary.py
python3 computations/verify_rank_fifteen_negative_x_arcs.py
python3 computations/verify_rank_fifteen_negative_x_fold.py
python3 computations/verify_rank_fifteen_negative_x_fan.py
python3 computations/verify_rank_fifteen_unit_quartic_towers.py
python3 computations/verify_rank_fifteen_genus_two_unit_line.py
python3 computations/verify_rank_fifteen_quadratic_verticals.py
python3 computations/verify_rank_fifteen_prime_parabolic_pencils.py
python3 computations/verify_rank_fifteen_quadratic_parabolic_pencils.py
python3 computations/verify_rank_fifteen_bivariate_word_paths.py
python3 computations/verify_rank_fifteen_bivariate_sixth_cut.py
python3 computations/verify_rank_fifteen_fixed_center_obstruction.py
python3 computations/verify_rank_fifteen_bivariate_seventh_cut.py
python3 computations/verify_rank_fifteen_truncated_central_freeness.py
python3 computations/verify_rank_fifteen_oriented_eighth_cut.py
python3 computations/verify_rank_fifteen_ninth_orbit_obstruction.py
python3 computations/verify_rank_fifteen_center_pencil.py
python3 computations/verify_rank_fifteen_full_center_cut.py
python3 computations/verify_rank_fifteen_vertical_one.py
python3 computations/verify_rank_fifteen_unit_residues.py
python3 computations/verify_rank_fifteen_vertical_unit_levels.py
python3 computations/verify_rank_fifteen_vertical_all_levels.py
python3 computations/verify_rank_fifteen_gaussian_pencil.py
python3 computations/verify_rank_fourteen_wall_centers.py
python3 computations/verify_rank_fifteen_quotient_residues.py
python3 computations/verify_rank_fifteen_regular_orbits.py
python3 computations/verify_rank_fifteen_small_residues.py
python3 computations/verify_rank_fifteen_branch_period.py
python3 computations/verify_rank_fifteen_coordinate_collisions.py
python3 computations/verify_rank_fifteen_pencil_linear_residues.py
python3 computations/verify_rank_fifteen_rational_residues.py
python3 computations/verify_varying_opposite_rank.py
python3 computations/verify_fixed_rank_fourteen.py
python3 computations/verify_uniform_quintic_companions.py
python3 computations/verify_uniform_fixed_quartic_residuals.py
python3 computations/verify_opposite_fixed_determinants.py
python3 computations/verify_fixed_quartic_ranks.py
python3 computations/verify_fixed_branch_regular.py
python3 computations/verify_klein_elliptic_center.py
python3 computations/verify_klein_full_rank_configurations.py
python3 computations/verify_rank_five_pencils.py
python3 computations/verify_matched_pencils.py
python3 computations/verify_incidence_cubics.py
python3 computations/verify_two_window_ext.py --output research/two_window_ext_certificate.json
python3 computations/verify_frobenius_tensor.py --output research/frobenius_tensor_certificate.json
python3 computations/verify_quantum_tensor_rank.py
python3 computations/verify_symmetric_congruence.py
python3 computations/verify_rank_two_exclusion.py
python3 computations/verify_quantum_indexed.py
python3 computations/verify_polynomial_ore.py --output research/polynomial_ore_certificate.json
python3 computations/verify_indexed_polynomial_ore.py --output research/indexed_polynomial_ore_certificate.json
python3 computations/verify_indexed_anticommutative_ore.py --output research/indexed_anticommutative_ore_certificate.json
python3 computations/verify_quantum_ore.py --output research/quantum_ore_certificate.json
python3 computations/verify_clifford_rank_slice.py --output research/clifford_rank_slice_certificate.json
python3 computations/verify_pure_anticommutative_slice.py --output research/pure_anticommutative_slice_certificate.json
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
latexmk papers/nc-p3.tex
latexmk research/plan_ja.tex
latexmk research/notes_ja.tex
~~~

check_project.py は参照・入力・JSONの構造検査であり数学的証明ではない。
最後の三次サンプル検算は三方向の観察に限る。PDFは build/ に生成しGit追跡しない。
文書・コード・入力・証明書をGitHub mainから再生成する。通常のfast-forwardで直接保存し、PR・force push・GitHub Actionsは使用しない。
GitHub保存成功時はLibraryへ重複保存しない。

対称チャート全体の表示像は `PGL4^3 × U` と同型な滑らかな69次元局所閉部分である。
四二次式完全交叉のテンソルcentroidが常にスカラーであることから、一意な三つの相対フレームを回収できる。
個々のindexed同型類は正則領域内でも非閉になり得る明示例を得た。
対称表示像の和集合は有限Frobenius領域内でも相対的に閉ではない。
[明示的な非Koszul境界](sections/nonkoszul_boundary.tex)は全ての四次検査を通るが、
ある始点の五次元数が56ではなく63となる。五次の7次元Koszul障害は、平方関係による6次元と明示的な7項の混合サイクルによる1次元で完全に生成され、3座標の線型汎関数で後者の独立性を検出できる（`prop:boundary-complete-quintic-cycles`）。一般の基点三因子対称境界も、左右反転対合によるTor下界5から非Koszulと証明した。対称閉包の非AS境界の全列挙は未解決だが、AS対象は対称チャートに尽きる。

三因子対称な有限Frobenius候補を、四層の切断代数のtrivial extensionとして記述した。
また一般点の有限次数障害を族全体の排除へ拡張する階数補題を得た。
二変数係数族では一般点の五次階数が77・36であることを有理関数体上で厳密計算した。
半連続性により、非零係数の特殊値も含めて標準型Koszul AS点を一様に排除した。
これは全ての一基点四二次式系の排除ではない。

二変数境界の厳密計算には追加の多項式演算バックエンドを使う：

```bash
python3 -m pip install -r requirements-symbolic.txt
python3 computations/verify_generic_one_point_boundary.py --output research/generic_one_point_boundary_certificate.json
```

係数は有理関数体 `QQ(u,v)` の元として扱い、数値標本や補間で一般点を認定しない。

交差項を含む、固定Boolean完全交叉上の全一基点二次式には、五係数域の全体で正則な評価フレームを構成した。
可換商のHilbert列は全域で `(1,4,6,4,1,1,...)`。Frobenius条件を通る点は対称チャート外にある。
この五係数域の一般五次$d_2$階数も77と厳密計算し、全域で標準型Koszul AS点を排除した。

```bash
python3 computations/verify_boolean_quadric_frame.py --output research/boolean_quadric_frame_certificate.json
python3 computations/verify_boolean_quadric_homology.py --output research/boolean_quadric_homology_certificate.json
```

後者は五変数の厳密消去で、上の二変数検算より重い。

相異なる八点の任意の三二次式完全交叉にも、一基点境界の評価フレームを拡張した。
点配置は列スカラーと射影変換の後に `[I | M]`, `M M^T = I` と表示できる（古典的な自己付随性に対応）。
一般位置の有理数標本でも、有限Frobenius・四次正確・五次障害7を検算した。
一般の八点配置全体も、左右反転対合による普遍的Tor下界によりAS対象から排除された。

```bash
python3 computations/verify_reduced_octad_frame.py --output research/reduced_octad_frame_certificate.json
```

基点を固定した四二次式系 `Gr(4,9)` の一般点を、積型の三二次式と第四式からなる8係数族と基底変換で被覆した。
接写像の20次小行列式 `-696/78125` により支配性を証明した。これはAS正則性や同型類次元の判定ではない。
八変数の五次階数が84未満なら、フレーム束上の半連続性で基点域全体を排除できる移送命題も得た。
八変数の正確な階数は未確定だが、一般の基点域全体のAS排除は独立の構造的証明で完成した。
一方、長さ2の基点部分スキームは構造的に排除した。相異なる二点も二重点も四次元数42以上を強制する。
四次検査を通る三因子対称な基点候補は一縮約点だけを持ち、縮約八点netによる表示と可換最小分解が適用できる。

左右反転対合により、戻り辺一回の商複体の負固有部分は次元 `60,200,135` を持つ。
Euler標数 `-5` が非線型 `Tor_(4,5)` を強制するため、基点を持つ三因子対称候補はすべて非Koszulである。
したがって、認証済みAS対象の対称枝所属は二つの可逆な正符号証明書で判定できる。

```bash
python3 computations/verify_product_net_dominance.py --output research/product_net_dominance_certificate.json
python3 computations/verify_product_net_homology.py --output research/product_net_homology_certificate.json
# One-return bimodule Tor: checked positive control and exact two-variable case
python3 computations/verify_trivial_extension_bitor.py --field ci --output research/trivial_extension_bitor_ci_certificate.json
python3 computations/verify_trivial_extension_bitor.py --field two --output research/trivial_extension_bitor_two_certificate.json
# Point/jet quartic obstruction and exact point-quotient regression
python3 computations/verify_basepoint_length_two.py --output research/basepoint_length_two_certificate.json
python3 computations/verify_trivial_extension_bitor.py --field two --frame ternary --point-quotient --last-only --output research/trivial_extension_bitor_quotient_certificate.json
# Signed-reversal theorem regression, also at a non-Boolean product-net point
python3 computations/verify_one_return_involution.py --output research/one_return_involution_certificate.json
python3 computations/verify_one_return_involution.py --fixture generic --output research/one_return_involution_generic_certificate.json
# All cuts through degree seven, exact differentials and two controls
python3 computations/verify_cut_tor_dimensions.py
python3 computations/verify_general_cut_tor.py --output research/general_cut_tor_certificate.json
# Tangent spaces at a square-free pencil (no rank-one quadric), and the four-square recognition fixture
python3 computations/verify_ci_middle_rigidity.py --output research/ci_middle_rigidity_certificate.json
python3 computations/verify_square_kernel_recognition.py --output research/square_kernel_recognition_certificate.json
# Boundary divisor and exponential growth from a single marker
python3 computations/verify_boolean_boundary_divisor.py --output research/boolean_boundary_divisor_certificate.json
python3 computations/verify_cubic_slice_growth.py --output research/cubic_slice_growth_certificate.json
# Finite Frobenius two-pair degeneration and both cyclic orders
python3 computations/verify_two_pair_boundary.py --output research/two_pair_boundary_certificate.json
python3 computations/verify_transfer_plane_patterns.py --output research/transfer_plane_patterns_certificate.json
# Control for the residual single-adjacent-negative linear matrix slice
python3 computations/verify_single_negative_slice.py --output research/single_negative_slice_certificate.json
# Exhaustive partial-alternation models and their Frobenius exclusions
python3 computations/verify_negative_witness_rigidity.py --output research/negative_witness_rigidity_certificate.json
# Optional exact eight-variable rank; exclusion theorem no longer depends on it
python3 computations/verify_trivial_extension_bitor.py --field eight --frame ternary --point-quotient --last-only --composition-check structural --output research/trivial_extension_bitor_eight_certificate.json
```

四周期左辞書順の二次PBW表示は、位置ごとに番号を付け直すと正規な二文字が $i\le j$ の一種類に限られる。四辺が有限なら、局所長3以上の点からなる四周期点列が必要である。対称チャートの共通順序基底での三角形表示は、二段の平方条件で判定できる。[証明と被約二十点の例](sections/periodic_pbw_graph.tex)

任意の有限Frobenius標準四窓も、一辺切断した32次元四層代数のtrivial extensionとして表せる。
四切断の戻り辺一回Torは内部次数7以下を覆い、最大の次元列は120,1040,2400,1600となる。
符号付きOre対照点とBoolean境界点で、七次までの全微分の階数とホモロジーを有理数上で検算した。
連続二窓の三次交叉が四次元のとき、そこから作る56次元代数が3-representation infiniteであることは標準AS延長の必要十分条件である。延長の同型類はこの二窓で決まる。この全次数条件を有限個の反復で判定する方法と、二窓の明示的な軌道分類は未解決である。

対称完全交叉チャート全体は三つの二対二階数を6以下にした射影表示空間の滑らかな69次元開部分である。対称類と可換多項式類を除いた低階数部分は有限Frobenius候補内で相対閉となる。この閉部分に標準AS対象が存在するかは未解決である。

対称成分の境界には、四位置の三対一の分割から四つ、二対二の分割から三つの、計七つの異なる68次元因子を構成した。それぞれの一般点の近傍で、低階数のAS点は対称チャートに尽きる。さらに任意の標準AS四窓では、連続二つの二次関係が非零の純テンソルとその共通の中央文字を持つことはない。Boolean境界では長さ $4m$ の次元が少なくとも $10^m$ となる。対称閉包内の標準AS点が対称チャートに属するための必要十分条件は、六つの位置対のどれか一つに可逆な正の対称化写像があることである。向かい合う一対だけが可逆だと仮定すると、二つの相補平面を使った経路商が指数成長を強制する。全ての正の対称化写像が特異な場合、同時極限から隣接階数1のpencilと二変数表示が強制される。その係数超平面を消す汎関数と両端支持による商では二辺の全関係が零となり、三文字ブロックの次元が4以上残る。長さ24で4096以上となって標準値2925に反するため、この残余も排除できる。したがって **対称閉包内のAS対象は対称完全交叉チャートに尽きる**（[完全な証明](sections/singular_transfer_planes.tex)、`thm:symmetric-closure-as-exact`）。ほかの非AS境界因子の有無と、閉包外の全分類は未解決である。

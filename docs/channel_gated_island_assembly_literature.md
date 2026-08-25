# Channel-gated island assembly hypothesis: primary-literature audit

## Bottom line

場所によって island floral syndrome が出たり出なかったりする理由を、文献に最も忠実に表すと次のようになる。

> 島嶼隔離が花形質を一方向へ変えるのではない。source region で実際に利用されていた送粉機能が、島への移住・存続過程で不足し、その機能が他の送粉者によって十分に代替されない場合にだけ、そのチャネルへ依存する植物が島フローラから選択的に欠落する、または繁殖保証型の植物が相対的に増える。

概念的には、島 *i* の floral filtering は次の積として考えるのがよい。

```text
filtering_i
  ~ source-channel exposure_i
  x isolation-associated functional deficit_i
  x (1 - effective replacement_i)
  x plant dependence / available trait variation_i
```

これは推定済みの数式ではなく、因果仮説の分解である。重要なのは、`Northern`、`Tropical`、`Southern` という地域名自体を原因にしないことである。地域差は、source pool、送粉者の進化史・分散能力、島での絶滅、代替guild、植物系統と繁殖系の違いをまとめた effect modifier である。

## 三つの channel gate

### Gate 1 — source region に機能チャネルが成立していたか

欠損を論じるには、まず比較対象となるsource regionで、その送粉者が存在するだけでなく、対象植物群へ有効な送粉機能を提供していた必要がある。Bombus は寒冷適応した主にHolarctic・山岳性の進化史をもち、地域的source opportunityが強く異なる（[Hines 2008](https://doi.org/10.1080/10635150801898912)）。一方、熱帯低地を含む例外的なNeotropical系統も存在するため、緯度だけからBombusチャネルを二値分類することはできない。

**支持される主張:** Bombusを候補チャネルとする生物地理学的根拠は地域依存である。

**未支持の飛躍:** 「北温帯ではBombusがすべての植物の主要送粉者」「熱帯・南半球ではBombusが構造的にゼロ」。

### Gate 2 — isolation に伴う実際の機能欠損があるか

海はBombusの移動を制約し得る。California沿岸の *Bombus vosnesenskii* では海洋障壁が遺伝子流動を強く制限した（[McFrederick & LeBuhn 2015](https://doi.org/10.1111/mec.13090)）。しかしBombusはCanary IslandsやMadeiraへ歴史的に海を越えて定着しており、島嶼隔離は自動的な不在を意味しない（[Widmer et al. 1998](https://www.nature.com/articles/6884070)）。したがって必要なのは、source poolと調査努力を考慮した `expected functional availability - observed functional availability` である。

Izu–Honshu系ではこの矢印に近い直接証拠がある。島で長舌pollinatorの機能多様性が低下すると、残存送粉者の花利用が変わった（[Hiraiwa & Ushimaru 2017](https://doi.org/10.1098/rspb.2016.2218)）。さらに40ネットワークの解析では、本島からの距離とともに送粉者の機能多様性とproboscis–corolla matchingが低下し、trait matchingは柱頭への花粉到達と対応した（[Hiraiwa & Ushimaru 2024](https://doi.org/10.1111/1365-2435.14527)）。*Campanula punctata* では、BombusがいないIzu諸島の多くでhalictid beesが優占し、Honshu・Oshimaの主な自家不和合から外側の島の自家和合へ変化した（[Inoue 1988](https://doi.org/10.1111/j.1442-1984.1988.tb00178.x)）。ただし後者は一系統の自然比較であり、Bombus欠損以外の歴史差を完全には除いていない。

**支持される主張:** 海洋隔離が特定の長舌・大型送粉機能を低下させる系は実在する。

**未支持の飛躍:** continent distanceだけから、各島で歴史的Bombus lossが起きたと推定すること。

### Gate 3 — 残存・代替送粉者が機能を埋めたか

送粉者の喪失後にネットワークは静止しない。Bombus一種を実験的に除去すると、残存pollinatorのresource overlapとconnectanceは増えたが、訪問される植物種は減少した（[Brosi et al. 2017](https://doi.org/10.1098/rsbl.2017.0243)）。すなわちrewiringは起きても、機能が完全に置換されるとは限らない。

他方、Galápagosでは昆虫食・種子食の陸鳥が広範な花資源を利用し、高度にgeneralizedなbird–flower networkを形成した（[Traveset et al. 2015](https://doi.org/10.1038/ncomms7376)）。これは島で新しい送粉チャネルが成立し得る直接証拠だが、同研究自身が鳥を「effectivenessにかかわらずpollinator」と扱っているため、完全な機能代替の証明ではない。Canary Islandsの鳥・トカゲでは訪花頻度と一訪花あたりの質が異なり、pollinator間の機能等価性は低かった（[Rodríguez-Rodríguez et al. 2013](https://doi.org/10.1007/s00442-013-2606-y)）。一般にも訪花者の約40%が有効な送粉者でなかった実測例があり、visitationだけでは代替を判定できない（[King et al. 2013](https://doi.org/10.1111/2041-210X.12074)）。

熱帯・南半球は固定的な「逆syndrome」ではなく、代替チャネルが残り得る反証領域である。熱帯海洋島ではhawkmothが長花筒・自家不和合の植物を送粉する直接例がある（[Xu et al. 2018](https://doi.org/10.1038/s41598-018-32143-5)）。Réunionでは鳥への送粉者交代と特殊な花形態が実測されている（[Micheneau et al. 2006](https://doi.org/10.1093/aob/mcl056)）。ただし、これらは「熱帯では蝶」「南半球では鳥」という地域一律の割当を支持しない。

**支持される主張:** 代替・rewiringの有無と質が、同じ隔離度でも植物側のfilteringを変え得る。

**未支持の飛躍:** 鮮明色・深い花・specialized formが残れば、それだけで鳥・蝶・蛾による機能代替が証明されたとすること。

## なぜ植物trait compositionが場所で異なるのか

### 1. source pool が場所ごとに異なる

島フローラは世界共通の母集団からランダム抽出されない。178 oceanic islandsと735 mainland regionsを用いた解析では、島ごとに推定したsource poolに対して植物科の過不足が大きく異なり、biotic-pollinated familiesは平均的に過少だった（[König et al. 2020](https://doi.org/10.1111/ecog.05174)）。Galápagosでも約39,000種の地域系統樹と気候適合性を用いた解析は、単純な長距離分散能力よりhabitat filteringがフローラ組成を強く説明した（[Carvajal-Endara et al. 2017](https://doi.org/10.1111/ele.12753)）。したがって同じ距離でも、sourceに小花・SC種が多い地域と、鳥・蛾送粉型が多い地域では、得られる島フローラが違う。

### 2. 送粉者と植物で分散・定着filterが異なる

全球の島フローラでは、animal-pollinated plantsなど相利共生者を要する植物が島で過少で、そのspecies deficitは小面積・遠距離と組み合わさって増えた（[Delavaux et al. 2024](https://doi.org/10.1038/s41586-024-07110-y)）。ただしこの研究はmutualist statusを文献・上位分類群から割り当てたassemblage研究であり、欠けた具体的pollinatorや花形質進化を同定していない。

### 3. reproductive assurance の必要性と実現可能性が植物系統で異なる

三つの大きな植物科・1,500種超の比較では、island speciesの66%がself-compatible、mainland-only speciesの41%がself-compatibleだった（[Grossenbacher et al. 2017](https://doi.org/10.1111/nph.14534)）。一方、1,752種・161科でself-compatibilityとautofertilityを連続量として扱った解析では、island colonistへの利点は「moderate」で、island endemicityは増やさなかった（[Razanajatovo et al. 2019](https://doi.org/10.1111/geb.12854)）。これはreproductive assuranceが一つのassembly filterであることを支持するが、SCをautonomous selfingや個々の島での進化と同一視できないことも示す。

### 4. 花の単純化は普遍的な進化則ではない

Pacific islandsの556種・136 sister contrastsでは、島の花が平均して小さくなる傾向はなく、縮小・拡大・無変化が島群と科によって異なった（[Hetherington-Rauth & Johnson 2020](https://doi.org/10.1086/709018)）。これは「隔離そのものが小花化を起こす」という普遍仮説に反する。同時に、島フローラで小花種が多く見える場合には、post-colonization evolutionよりecological assembly filteringが有力になり得る、というChapter 1の組成的解釈と整合する。

## Bombus environmental suitability だけで deficit を推定できるか

**できない。** 環境ニッチ適合度はsource-channel成立より手前の
「recipient environmentで成立可能か」という補助的opportunityの近似であり、
Gate 2のobserved deficitではない。

1. 本リポジトリの実装はMaxentではなく、出現地点の9環境軸から作るwinsorized・standardized・ridge-regularized Mahalanobis ellipsoidである。95%経験境界で適合度0.5になる連続スコアは再現可能だが、出現データ由来の環境類似度であって、調査済み不在、局所占有、個体数を直接推定しない。一般にpresence-onlyデータではsampling biasが環境変数と共変すると推論が歪む（[Yackulic et al. 2013](https://doi.org/10.1111/2041-210X.12004)）。
2. Habitat suitabilityと局所abundanceの関係は一般に弱い。246 mammal speciesと158 tree speciesの比較では、niche-model suitabilityはabundanceと概して無関係だった（[Dallas & Hastings 2018](https://doi.org/10.1111/geb.12820)）。
3. たとえBombusが存在しても、対象植物へのvisitation、pollen deposition、seed productionへの寄与は別変数である（[King et al. 2013](https://doi.org/10.1111/2041-210X.12074)）。
4. 複数Bombus種の最大適合度は「少なくとも一種に気候が似る」ことを表すだけで、source regionにいた機能群の組成、海を越える到達可能性、島での営巣資源、個体数、長舌機能を保存しない。

したがって、現在の `environmental_compatibility_max` が花サイズ勾配を条件づけなかったことは、channel-gated hypothesis全体の反証ではない。支持されなかったのは、**この気候適合度をrecipient opportunity gateとして使い、隔離勾配が高opportunity側で強まるとする特定の測定モデル**である。

さらにstatus・genus固定後のcheckpointでも、generalized formとSCの予測力は
改善せず、plain colourの小さな改善は単純なdeficit予測と逆方向だった。
したがってこのproxyはChapter 1の直接説明変数ではなく、保存された探索的・
補助的診断である。

## Pollination syndrome 文献をどう使うか

417種の定量レビューでは、複数花形質からなるsyndromeはeffective pollinator
groupを一定程度予測し、とくにpollinator-dependent taxaと熱帯植物で予測性が
高かった（[Rosas-Guerrero et al. 2014](https://doi.org/10.1111/ele.12224)）。
しかし全球テストではtraditional syndromeがdominant functional pollinatorを
しばしば誤分類した（[Ollerton et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2701765/)）。
近年のレビューも、syndrome概念を捨てるのではなく、colour単独ではなくshape、
reward、corolla widthなどを含むmultivariate trait dataとして再検証する必要を
強調する（[Dellinger 2020](https://doi.org/10.1111/nph.16793)）。

したがって本研究で許されるのは次の照合だけである。

```text
supported multivariate residual trait vector
    <-> literature-defined pollinator-functional expectation
```

`red/pink -> bird` や `small flower -> Bombus absence` のような単一形質割当は
行わない。照合はmechanistic evidenceではなく、Chapter 2で区別すべき候補機構
を生成するDiscussion上のconcordanceである。

## 博士論文内での位置づけ

`Channel-gated island assembly` はChapter 1単独の直接仮説ではなく、三章を
統合する博士論文全体の概念仮説とする。Chapter 1が直接識別するのは、
status・lineageを越えて isolation-associated trait signal がどこに残るか
までである。

| 対象 | Chapter 1で識別可能か | 現在の証拠 |
|---|---:|---|
| isolationがfloristic statusを変えるか | **可能** | Northernでregional endemicityが頑健に増える |
| broad syndromeがstatus・genus固定後も残るか | **可能** | native non-endemicのplain/generalized/SCは概ねnull |
| fine categoryが残るか | **可能** | Northern all-nativeの小さなred/pink減少のみ比較的頑健 |
| coherent multivariate syndromeがあるか | **現時点で否定的** | confirmatory form residualはなく、単一colour結果だけではsyndromeと呼べない |
| source側でBombus channelが成立していたか | **不可** | 気候opportunityは有効相互作用を測らない |
| 島でBombus機能が欠損したか | **不可** | effort-corrected occurrence、abundance、functional traits、interactionがない |
| 鳥・蝶・蛾などが機能代替したか | **不可** | visitationとper-visit effectivenessが島別にない |
| assembly filteringか島内進化か | **不可** | population/temporal contrastsがない |

Chapter 1の現在の主張上限は次である。

> Isolation is robustly associated with floristic and endemic-lineage turnover,
> whereas little broad floral or reproductive trait signal remains after
> floristic status and genus composition are represented. The surviving
> Northern colour residual does not identify a pollinator channel.

pollination syndrome文献は、単一形質から送粉者を割り当てるためではなく、
複数の独立した形質が残った場合の**multivariate concordance**にだけ使う。
NorthernではBombus/large long-tongued function、Tropical/Southernでは鳥・蝶・
蛾などを同じ証拠基準でDiscussionの候補機構として照合する。現在の単一の
red/pink残差だけでは、この照合段階へ進まない。

channel-gated仮説を直接検証するChapter 2では、地域帯を説明変数として増やすのではなく、各島について次を構築する必要がある。

1. source regionでのBombus機能的availability（対象植物への有効性を含む）
2. sampling effortを考慮した島のobserved Bombus occurrence / abundance / functional composition
3. `expected - observed` としてのBombus functional deficit
4. 鳥・蝶・蛾・他bee・fly等の `visitation x per-visit effectiveness` によるreplacement
5. source poolに対する島フローラのtrait over-/under-representation

これにより「Northernだからsyndromeが出る」から、**成立可能だったチャネルが失われ、代替されなかった場所でsyndromeが出る**という検証可能な説明へ移行できる。

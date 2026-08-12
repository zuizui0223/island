# 形質取得の戦略転換案（PR #132 の停滞に対して）

## 結論

取得を止める必要はない。**分母を変える**。

現行は `106,295 accepted species × 3 axes = 318,885 cells` を等重みで埋めにいっている。
確証モデルが読むのは種ではなく**島ごとの組成**であり、島組成は
`island × species` の帰属（本書では flora mass と呼ぶ）上で平均される。
この二つは同じものではない。993島に分布する *Trifolium repens* と単一島固有種は、
現行分母では同じ 1 cell、解析寄与では 3 桁違う。

凍結済みの flora universe で測ると、島重み付けの標的化は同一種数予算で
**桁違いに解析可能島を増やす**。数値は本 PR で追加した
`island-v2-acquisition-value` により再現できる（ネットワーク不要、入力は凍結済み artifact）。

| 種予算 | 戦略 | flora mass 被覆 | 島flora被覆の中央値 | 解析可能島 |
|---:|---|---:|---:|---:|
| 5,000 | 島重み付け | 60.7% | 74.2% | **1,921** |
| 5,000 | 種数等重み | 4.5% | 3.6% | 0 |
| 10,000 | 島重み付け | 72.8% | 91.4% | **2,307** |
| 10,000 | 種数等重み | 8.8% | 8.0% | 0 |
| 50,000 | 種数等重み | 42.5% | 42.4% | 168 |

（解析可能 = その島の realised flora のうち 30種以上かつ 50%以上が被覆。
全 4,505 島、flora mass 合計 1,039,757。閾値は `config/island_weighted_acquisition.yml` で宣言。）

種数等重みは 50,000種（5倍予算）を取っても 168島にしか届かない。
島重み付けは 10,000種で 2,307島に届く。**同一予算で 20〜50倍の効率差**である。

## 1. なぜ停滞して見えるのか

停滞は取得能力の問題ではなく、進捗指標が解析目的と一致していないことの帰結である。

### 1.1 直接証拠の限界収量が分母に対して小さすぎる

PR #132 に記録された実測：

| 実行 | 投入 | 直接 species × trait 純増 |
|---|---|---:|
| 2026-07-25 pilot | 300種 / 600ページ | 100 |
| 多ドメイン inventory | 4ドメイン / 600ページ | 158 (cells) |
| Europe PMC batch2 (08-07) | 100 query / 957 全文 | 5 |
| high-leverage batch (08-08) | 個別review 4件 | 2 |

未解決は 222,568 cells 規模である。直近の収量では到達しない。
これは「もっと走らせれば埋まる」種類の差ではない。

### 1.2 報告被覆のほぼ全部が genus 推論の increment である

08-08 の `+735` cells は内訳が `Validated Low +733 / direct species × trait +2`。
2件の直接事実が 2つの genus rule を eligible にし、733 cells を生んでいる。
この倍率は逆方向にも働く。同じ PR で `Poa` 245行、rebuild 一回で 12,306行が
invalidate されている。

### 1.3 その結果、指標が取得と無関係に動く

08-07 コメントの strict 被覆は `96,317 / 318,885`、08-08 コメントの比較基線は
`75,852`。間に取得イベントはない。原因が rebuild なのか contract 変更なのかは
CI artifact が期限切れのため本書では断定しないが、
**約20,000 cells の改訂は PR #132 の全取得バッチの合計より大きい**。
この状態の数値は進捗信号として使えない。突合が先である。

### 1.4 ソース駆動で標的が選ばれている

実際に処理された分類群と、その島フットプリント：

| 分類群 | 種数 | flora mass |
|---|---:|---:|
| Coelogyne | 330 | 0.051% |
| Bulbophyllum | 981 | 0.154% |
| Diospyros | 445 | 0.151% |
| Miconia | 483 | 0.118% |
| Panicum | 159 | 0.108% |

Mikawa / Tsukuba / botanic.jp / PlantNET NSW という地域フロラ起点の設計上、
標的は「そのソースに載っている種」であって「島組成を動かす種」ではない。
狭域固有種に偏るため、島重み付けでは種数等重みより**さらに悪い**可能性がある。

一方で最上位の標的は次のとおり：*Trifolium repens* (993島),
*Achillea millefolium* (954), *Stellaria media* (904), *Plantago major* (894),
*Prunella vulgaris* (822), *Plantago lanceolata* (808)。
これらは広域普通種であり、**既に repo が ingest workflow を持つ bulk source
（USDA PLANTS, AusTraits, GIFT, EOL TraitBank, Ecoflora, BiolFlor 系）が
最もよく収録している種**でもある。最高価値の標的が最も安い標的である。

## 2. 転換案

### 2.1 主指標を差し替える（撤回ではなく併記）

`318,885 cells` の strict 被覆は submission contract の evidence tier 表として残す。
その上に、**解析可能島数**を運用主指標として置く。

- 分子: readiness gate を満たす島
- gate: `min_covered_species=30` かつ `min_covered_fraction=0.5`
- 分母: 4,505 島（realised flora を持つ島）

1バッチごとに「解析可能島が何島増えたか」を読む。0島なら、そのバッチは
cells を増やしていても解析を進めていない。

### 2.2 標的を島重み付け順に固定する

`island-v2-acquisition-value run` が出力する `species_acquisition_value.csv.gz` の
`priority_rank` を、全取得レーンの標的キューにする。
第1期の予算は **上位10,000種**（flora mass 72.8%、解析可能島 2,307）。

### 2.3 取得レーンを scope class で分ける

上位10,000種の内訳：

| scope class | 種数 | flora mass |
|---|---:|---:|
| biotic_floral_evidence_required | 7,363 | 51.9% |
| anemophilous_candidate | 1,799 | 15.9% |
| non_angiosperm_out_of_floral_scope | 838 | 5.1% |

- **biotic（7,363種）** — 実証拠が必要な本命。bulk join を先に尽くし、
  残差だけを Web/全文へ回す。現行の抽出器はそのまま使える。
- **anemophilous（1,799種）** — Poaceae, Cyperaceae, Juncaceae, *Plantago* 等。
  clade 内でほぼ不変な形質について、既存の genus Validated Low と同じ
  dominance + LOO 機構を **family 水準で明示的に検証**することを提案する。
  合格した rule のみ、`family_validated` という**別 tier** として報告する。
  strict 分子には入れない。反例1件で rule 全体を失効させる。
  これは禁止されている暗黙の family inference ではなく、宣言・検証・分離報告
  された規則である。適用是非は人間の判断に委ねる。
- **non-angiosperm（838種）** — 花形質の対象外。
  `not_applicable_taxonomic_scope` として明示し、**未解決から外す**。
  現行 acquisition-rate ledger は既にこの扱いをしている。整合させる。

3つ目だけでも、上位10,000種の 8.4% は「埋める対象ではない」ことが確定する。

### 2.4 Validated Low の役割を縮小する

genus rule は 2件の直接事実から 733 cells を生む。この増幅は
被覆報告の主成分であってはならない。

- Validated Low は感度解析 tier に限定し、主指標（解析可能島）の分子から外す。
- 主指標は direct 証拠のみで計算する。
- genus rule は残す。ただし「被覆が増えた」ではなく
  「rule が n 島の組成推定を安定させた」として報告する。

これにより 1.2 と 1.3 の指標振動が主指標に伝播しなくなる。

## 3. 実行順

1. `island-v2-acquisition-value run` を CI 化し、`priority_rank` artifact を固定する。
2. 既存 bulk source（USDA PLANTS, AusTraits, GIFT, EOL, Ecoflora）を
   **上位10,000種に限定して**再走査し、解析可能島の純増で測る。
   新規取得ゼロで済む可能性が高い区間である。
3. 残差 biotic 種のみを既存の Web / Europe PMC レーンへ回す。
   標的キューは 2.2 の順序に従う。
4. anemophilous family rule を dominance + LOO で検証し、
   別 tier として提案する（採否は別判断）。
5. non-angiosperm を未解決分母から外し、`318,885` の内訳を再掲する。
6. 1.3 の約20,000 cells 差を突合し、strict 被覆の履歴を1本にする。

## 4. 変えないもの

- species/synonym-direct 要件、exact quote と stable URL。
- cultivar / hybrid の fail-closed。
- 暗黙の family inference と global fallback の禁止。
  2.3 の family rule は**明示・検証・別tier**であり、この禁止に抵触しない。
- missing は missing のまま。scope class は取得計画のラベルであって
  形質値でも生物学的不在でもない。
- production promotion は独立した人手 review を要する。

## 5. 縮小を提案するもの

- 地域フロラ inventory の網羅走査（Mikawa / Tsukuba / botanic.jp / PlantNET NSW）。
  標的キュー経由の lookup としてのみ残す。
- 単一島固有種への個別取得。58,660種（全種の50.9%）が単一島で、
  flora mass の 5.6% しか持たない。第1期の標的から外す。
- `cells 純増` を成果として報告する運用。解析可能島の純増に置き換える。

## 実装済みの計測器

本提案の各項目は、主張ではなく**実行可能な計測**として実装してある。
すべてネットワーク不要、入力は凍結済み artifact。

### 標的キューと planning frontier（2.2、2.3）

```bash
island-v2-acquisition-value run \
  --config config/island_weighted_acquisition.yml \
  --output-dir data/v2/staging/traits/acquisition_value
```

出力は `species_acquisition_value.csv.gz`（`priority_rank` 付き標的キュー）、
`acquisition_strategy_frontier.csv`（本書の表）、
`target_scope_breakdown.csv`（レーン内訳）、`acquisition_value_summary.json`。

### 主指標の実測（2.1、2.4）

```bash
island-v2-acquisition-value evaluate \
  --evidence-csv <accepted evidence ledger> \
  --output-dir <out>
```

軸ごと・tier ごとに解析可能島数を出す。`direct_only` が primary、
`direct_plus_validated_low` は sensitivity として**分離報告**される。
family inference と global fallback はどちらの tier でも admit されない
（テストで固定）。

### source の受入判定（実行順 2、3）

```bash
island-v2-acquisition-value source-yield \
  --candidate-csv <source candidates> \
  --baseline-csv <current accepted species> \
  --budget 10000 \
  --output-dir <out>
```

1つのソースが baseline に対して**解析可能島を何島増やしたか**を返す。
`priority_head_hit_rate` は、そのソースの新規種のうち上位予算に入った割合。
種数や cells ではなくこの純増が受入判定になる。

### 風媒 clade rule の検証（2.3）

```bash
island-v2-anemophilous-clade-rule build \
  --evidence-csv <accepted evidence ledger> \
  --unresolved-csv <unresolved species x trait> \
  --output-dir <out>
```

宣言済み clade のみを対象に support / dominance / 源流数 /
species LOO / lineage LOO を採点する。閾値は genus track より厳しく
（`min_species=20`, `min_dominance=0.95`, `min_source_lineages=3`,
両 LOO `0.95`）、1つでも外すと `eligible=False` と失格理由が残る。
出力は `family_validated` tier で、strict 分子には決して入らない。
**この tier を採用するか否かは submission contract の人間判断である。**
本モジュールは「何が通るか」を報告するだけで、採用はしない。

### CI

`.github/workflows/measure-island-weighted-acquisition.yml` の
`validate` job は PR ごとに lint / test と frontier 再構築を行う。
`measure` job は dispatch 限定で、pin した evidence artifact を解決して
上記の実測を実行し、artifact として上げる。

## 最初の一手

`measure` job に現行の integrated coverage artifact を渡して
`evaluate` を走らせる。取得は不要で、いま何島が解析可能かが直接出る。
**これが本提案の最初の実測になる。**

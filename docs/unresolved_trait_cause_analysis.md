# 未解決 222,636 cells の原因分析

## 対象

repo にコミットされている最新の strict coverage ledger：

```
data/v2/staging/traits/direct_llm_pilot/
  20260807_europe_pmc_reproduction_integrated/strict_species_axis_coverage.csv.gz
```

`106,295 species × 3 axes = 318,885` 行を完全に保持しており、
flora master (`island_taxa.csv`) と 100% 突合する（未照合 0 件）。
PR #132 の 08-08 以降のコメントが参照する ledger は artifact 側にあり
本書では扱わない。以下はすべてこの ledger 上の実測である。

## 1. 全体像

| 状態 | cells | 割合 |
|---|---:|---:|
| 未解決 | 222,636 | 69.82% |
| Validated Low | 51,616 | 16.19% |
| Medium | 30,736 | 9.64% |
| High | 13,897 | 4.36% |

解決済み 96,249 のうち **53.6% が Validated Low（genus 推論）**。
直接証拠（High + Medium）は 44,633 cells で、**分母の 14.0%** にとどまる。

## 2. 未解決の 78% は島をほとんど持たない種である

未解決 cells を種の島フットプリントで分解する：

| n_islands | 未解決 cells | 割合 | 種数 |
|---|---:|---:|---:|
| 1 | 126,205 | 56.7% | 53,569 |
| 2–3 | 47,751 | 21.4% | 21,277 |
| 4–10 | 30,603 | 13.7% | 14,534 |
| 11–50 | 14,763 | 6.6% | 7,910 |
| 51–200 | 2,883 | 1.3% | 1,864 |
| 201+ | 431 | 0.2% | 304 |

**未解決の 78.1% は 3島以下にしか分布しない種**である。
科別でも Orchidaceae 8.26%、Rubiaceae 6.21%、Melastomataceae 2.29% と
狭域固有種の多い群が上位を占める。

一方、島フットプリント上位10,000種（priority head）は：

- cells 27,285 のうち **解決済み 16,010（58.7%）**
- 未解決 11,275 は、**全 backlog の 5.1% にすぎない**

つまり「222,636 cells が未解決」という数字が示しているのは、
取得が失敗しているという事実ではなく、
**分母の大半が島組成をほとんど動かさない種で構成されている**という事実である。

## 3. 実際に何島が解析可能か

gate は「その島の realised flora のうち 30種以上かつ 50%以上が被覆」。
全 4,505 島に対する実測：

| 軸 | 直接証拠のみ | Validated Low 込み |
|---|---:|---:|
| flower_colour | **2,085** (46.3%) | 2,109 |
| floral_structural_complexity | **782** (17.4%) | 2,103 |
| reproductive_assurance | **461** (10.2%) | 1,056 |
| **3軸すべて** | **422** | — |

これが現状である。**直接証拠で3軸そろう島は 4,505 島中 422 島**。

## 4. どこが詰まっているか

### 4.1 花色は終わっている

flower_colour は直接証拠だけで 2,085 島。
Validated Low を足しても **+24 島**しか増えない。
この軸に追加取得の余地はほぼない。

### 4.2 律速は花構造、その次が繁殖保証

- 構造軸の直接 782 島が**天井**である。
- 色は ready だが構造が ready でない島が **1,303 島**ある。
- 構造軸は Validated Low が +1,321 島を生んでいる。
  つまり現行報告の構造被覆は**ほぼ genus 推論で出来ている**。
  一括 invalidate で消えるのはここである。

### 4.3 繁殖保証だけで止まっている島が 360 ある

色と構造を直接証拠でクリア済みなのに、繁殖保証だけ足りない島が **360**。
ここを埋めれば3軸 ready は **422 → 782（+85%）** になり、構造の天井に達する。

## 5. 具体的な標的

### 5.1 繁殖保証（360島を解放する）

対象種 11,597。上位200種で 11,671 island-slot、上位500種で 22,050 をカバーする。
最上位：

| 解放される島数 | 種 |
|---:|---|
| 117 | *Taraxacum officinale* |
| 104 | *Chenopodium album* |
| 103 | *Lonicera japonica* |
| 97 | *Cynodon dactylon* |
| 94 | *Euphorbia maculata* |
| 90 | *Pueraria montana* |
| 88 | *Melia azedarach* |
| 86 | *Echinochloa crus-galli* |
| 85 | *Juncus effusus* |

### 5.2 花構造（天井を上げる）

対象種 43,329。上位200種で 39,301 island-slot。最上位：
*Poa pratensis* (389)、*Juncus bufonius* (384)、*Luzula multiflora* (369)、
*Chamaenerion angustifolium* (361)、*Plantago maritima* (358)、
*Sonchus oleraceus* (351)、*Caltha palustris* (347)、*Drosera rotundifolia* (337)。

## 6. 原因

停滞の原因は取得能力ではなく、**最適化している量が律速要因と一致していないこと**である。

1. **分母が狭域固有種で希釈されている。** cells 純増を成果とする限り、
   78% を占める3島以下の種が常に最大の「残り」に見える。
   実際にはそこを埋めても解析可能島はほぼ増えない。

2. **軸を区別していない。** 花色は既に終わっており、律速は構造と繁殖保証である。
   3軸を合算した被覆率は、終わった軸が詰まった軸を隠す。

3. **島の到達状態を見ていない。** 360島が1軸だけで止まっていることは、
   species×axis の集計からは見えない。

4. **標的が literature 密度と逆相関している。** 実際に処理されたのは
   *Coelogyne*、*Bulbophyllum*、*Diospyros* といった狭域・低文献量の群である。
   一方で最大の律速標的は *Taraxacum officinale*、*Chenopodium album*、
   *Poa pratensis*、*Sonchus oleraceus* ——
   繁殖様式も花構造も標準的な文献に記載がある普通種である。

5. **報告被覆が genus 推論に依存している。** 構造軸の Validated Low は
   +1,321 島相当を生んでおり、一括 invalidate のたびに指標が大きく動く。
   直接証拠 44,633 cells（分母の 14.0%）が実体である。

## 7. 直ちに取るべき順序

1. **繁殖保証の上位200種**を取得する。360島が解放され、3軸 ready が 422 → 782 になる。
2. **花構造の上位200種**を取得する。天井そのものを上げる。
3. **花色への追加取得を止める。** 残余は +24 島しかない。
4. 主指標を cells 純増から **3軸 ready 島数**へ差し替える。
   422 という数字を毎バッチ更新すれば、進捗が初めて可視になる。
5. 3島以下の種（未解決の 78.1%）は第1期の標的から外す。

## 再現

```bash
python3 - <<'PY'
# ledger は上記パス、flora master は data/v2/staging/gbif/collected/ 配下。
# 集計は accepted_species x axis の quality 列のみを使う。
PY
```

本書の数値はすべて repo 内の凍結済みファイルから再計算できる。
ネットワークも artifact ダウンロードも不要である。

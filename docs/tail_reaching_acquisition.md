# テールに届く取得（データ増加策 1・2・3）

`docs/coverage_confound_analysis.md` の診断に対する実装。狙いは cells の総量ではなく
**被覆と隔離度の相関（rho = −0.369）を平らにする方向にデータを増やすこと**。

| 手 | 新規取得 | テール到達 | rho への効果 | 状態 |
|---|---|---|---|---|
| 1 捨てた証拠の回収 | 不要 | 中 | 中立 | 実装済・実測済 |
| 2 標本ラベル採掘 | bulk・無料 | **高** | **平坦化** | 実装済・CI 待ち |
| 3 島嶼フロラ／原記載 | IPNI＋BHL | **高** | **平坦化** | 原記載は実装済・CI 待ち／島嶼フロラは許諾待ち |

## 1. 捨てた証拠の回収（実測済み）

棄却台帳を読み直し、正当な棄却と自前の取りこぼしを分ける測定。
実装は本ブランチで一度landしたのち上流で superseded として削除された（`f378497`, `9ecacc3`）。
以下は削除前に取得した実測値で、値自体は有効。正準実装での再測定が要る。

**名前照合が主戦場。** 既存 bulk 2ソースの棄却 1,917行のうち
**1,823行（95.1%）** が GBIF EXACT × ACCEPTED × usage key 付き。
**1,753 個のキー**を1回ずつ引けば accepted name に解決できる。
family conflict、非種階級、FUZZY、HIGHERRANK は決して回収しない。

**抽出は理由ラベルより小さい。** 引用文を読んで分類すると 251行の内訳は：

| 分類 | 件数 | 割合 |
|---|---:|---:|
| 正しく棄却（別器官の計測） | 97 | 38.6% |
| 正しく棄却（記述テンプレートの選択肢列挙） | 51 | 20.3% |
| **回収可能（非英語の明示記載）** | 22 | 8.8% |
| 要人手確認 | 81 | 32.3% |

style は花ではないので大半は正しく捨てられている。回収可能なのは多言語レーンで、
**引用文に語が実在することを要求**する（言語タグだけでは回収しない）。
`language` 列は `en/English/EN/eng/fre/es/Spanish/空` と不統一なので先に正規化する。

## 2. GBIF 標本ラベル採掘

`island-v2-gbif-label-mining`。テールに届く唯一の大規模ソース。

未解決テールは GBIF 記録中央値4件の島嶼固有種で、文献はほぼ存在しない。
だが**その4件はすべて押し葉標本**であり、採集者は乾燥で失われる花色を
現場でラベルに書く。よってラベルは：

- テールに届く
- 同定済み標本に付随するので **species-direct が構造的に保証**される
- 無料・bulk・crawl 不要（DwC-A に標準で含まれる）

**器官紐付けを 1 と同じ規律で行う。** これは理屈ではなく実測に基づく —
WFO 棄却の 38.6% が別器官の値を全体花の値として読んだものだった。よって：

- 色語は花器官語から一定距離内になければ棄却
- 最近傍の器官が fruit / leaf / calyx 等なら**その場で棄却**（同距離なら競合器官を優先＝保守側）
- 「flowers not seen」「sterile」「? white」は無効化
- 「white to pink」は最初の語を採らず `multicolored_variable`
- **重複標本を1採集イベントに畳む**。1回の採集が複数標本館に分配されているだけのものを
  独立証拠として数えない

出力は verbatim 引用付きの unreviewed 候補。候補のない種は
「色がない」のではなく未解決。

## 3. 島嶼フロラと原記載

registry は `config/island_flora_sources.yml`。ここは**2レーンに分かれる**。

### 3a. 原記載（実装済み）

`island-v2-protologue-acquisition`。IPNI（書誌）→ BHL（スキャン）の2段。
両方とも registry の `permitted_routes` にある official API なので、
許諾確認待ちの島嶼フロラと違って有効化してよい経路。

テール 74,846種は定義上ほぼ島嶼固有種で、記録4件の種にとって
**原記載が現存する唯一の花形態記載であることが多い**。しかも
「名前が有効発表されている以上、原記載は必ず存在する」という点で、
文献の有無が種によってばらつく他のソースと性質が違う。

凍結マスタからテールを切ると **83,119種**（n_islands ≤ 3、記録数中央値5）。

実装上、他のテキスト採掘と決定的に違う点が3つある。

**多言語は要件であって好みではない。** 原記載はラテン語・仏語・独語・西語が
普通で、英語前提が失敗することは既に実測済み（WFO 棄却のうち自動回収可能な
唯一の層 8.8% が非英語記載だった）。照合は fold したテキスト
（小文字化・ダイアクリティカル除去・ß→ss）に対して行う。

**ラテン語は語形を展開する。** 形容詞は名詞に一致するので、実際に出る語形は
文が決める — `corolla alba` だが `flores albi` だが `floribus albis`。
原記載は複数形を辞書形と同じくらい使うので、辞書形だけ列挙すると
大半を取りこぼす。よって語幹に閉じた語尾集合を適用する
（161個の literal → 335語形）。実際これはテストで捕まった欠陥で、
最初は `lutei` も `rubri` も引けていなかった。

**文境界を越えない。** 器官探索を文（`.` `;` 改行区切り）の内側に閉じる。
`Folia coriacea, subtus fusca. Flores albi.` で「fusca」から見ると
自分の `Folia` より次の文の `Flores` の方が近く、境界がないと
**葉の褐色が花色に混入する**。これも実装中に実データ形のサンプルで発覚した。
WFO 棄却の 38.6% がまさにこの器官取り違えである。

保守側に倒す設計なので、既知の取りこぼしがある。
`Corolla alba demum rosea, calyce viridi` は `rosea` が `corolla` より
`calyce` に近いため落ち、`white` だけが残る。誤った帰属を作るより
ニュアンスを失う方を選んでおり、テストで意図として固定してある。

著作権は二重ゲート。BHL が public domain と宣言し、かつ刊行年が 1928 年以前。
**保護期間内のものは skip であって scrape ではない。**
lineage 5項目（source_url / source_citation / license / retrieval_date /
exact_supporting_quote）を欠く行は降格ではなく棄却。

### 3b. 島嶼フロラ（未着手・許諾待ち）

Flore de Madagascar、Flora Malesiana、Flora of Socotra など。
registry に access 種別を付けてあるが、多くが `needs_review`。
**各ソースの許諾確認が先**で、こちらでは一切取得していない。
robots.txt に反する crawl と無許可の portal scrape は明示的に禁止。
WFO portal は禁止済みとして記録し、再発見されないようにしてある。

## 実行

```bash
# 2（GBIF download 認証が要る。workflow から dispatch）
island-v2-gbif-label-mining run --occurrences-csv <occ.csv> --output-dir <out>

# 3a（BHL_API_KEY が要る。workflow から dispatch）
island-v2-protologue-acquisition queue --output-csv <queue.csv>
island-v2-protologue-acquisition run --queue-csv <queue.csv> --output-dir <out>
```

CI: `mine-gbif-specimen-labels.yml`、`acquire-protologue-traits.yml`
（いずれも validate は PR ごと、取得は dispatch）。

1 の測定を再現するには正準の bulk-recovery 実装が要る。

## 未確定

- 1 の解決収量（1,753キーのうち何件が island master に着地するか）は
  GBIF API が要るため CI 実行待ち。作成環境からは到達できない。
- 2 のラベル記載ヒット率は未知。パイロットで測る。
- 3a のヒット率は未知。`BHL_API_KEY` が未設定なので一度も走っていない。
  効かない側の見込みも書いておく：BHL に該当ページのスキャンが無い、
  OCR が古い植物ラテン語で崩れる、citation の page が BHL の
  `PageNumbers` と一致しない、のいずれかで大半が落ちる可能性がある。
  まず `--limit 500` のパイロットで棄却理由の内訳を見るべき。
- 3b は access 欄が `needs_review` のものが多い。**各ソースの許諾確認が先**で、
  こちらでは一切取得していない。

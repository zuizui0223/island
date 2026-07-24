# Trait acquisition: next implementation steps

## Decision

既存artifactの統合後は、次の順で追加取得を進める。順位は「実装を自動化できるか」「現在の3軸（花色・花構造・繁殖保証）へ直接寄与するか」「source lineageとライセンスを再現可能に保存できるか」を合わせたもの。期待値は事前評価であり、採用前に固定した未解決種で純増を測る。

| 順位 | 改善案 | 実装可能性 | 3軸への期待 | 判断 |
|---:|---|---|---|---|
| 1 | synonym検索の強化 | 高 | 高（全軸） | 静的な版付きalias索引を一度作れば、既存providerを追加走査せず再利用できる。 |
| 2 | tube depth・flower size・reward・inflorescence display抽出 | 高 | 高（主に花構造） | 既取得本文へ決定的抽出を追加でき、新規Web取得が不要。 |
| 3 | Europe PMC full-text XML・supplementary table | 高〜中 | 中〜高（繁殖保証、花構造） | 現行のtitle/abstract discoveryをOA本文へ自然に拡張できるが、表形式の多様性を扱う必要がある。 |
| 4 | botanical Latin・多言語抽出 | 中 | 中（花色、花構造） | 取得路は既にある。言語別辞書と検証用gold setがボトルネック。 |
| 5 | TRY対象形質export | 外部申請後は高 | 非常に高（全軸） | export取込は容易だが、登録済みPIによる手動requestと利用条件確認をCIから代替できない。 |
| 6 | WFO description全件取得 | 低 | 高い可能性（花色、花構造） | 公開portalの自動crawlは禁止。公式content exportか元providerのDwC-A確保が先。 |

## 1. Synonym検索の強化

WFO Plant Listの版付きsnapshotはJSON、Darwin Core、SQLとaccepted/synonym関係をまとめて取得できるため、オンライン照会より先にローカルalias索引へ変換する。[WFO Plant List June 2026 snapshot](https://zenodo.org/records/20782718) をpinし、WFO ID、snapshot DOI、accepted name、synonym、authorshipを保存する。snapshotで解決しない名前だけをGBIFの `species/match` と `species/{usageKey}/synonyms` で補完する。[GBIF Species API](https://techdocs.gbif.org/en/openapi/v1/species) は多数照会時に429を返し得るため、静的snapshot優先とbounded backoffを守る。[GBIF API reference](https://techdocs.gbif.org/en/openapi/)

実装 seam は `src/island_v2/alias_wikimedia_discovery.py` の入力alias表である。`accepted_species, search_name, name_role, taxon_id, source_snapshot` を持つ版付きartifactを生成し、既存のWikidata P225検証を維持する。fuzzy、高次分類群、homonym未解決は自動採用せずreview候補へ残す。providerへはaliasで検索しても、結果は常に元の `accepted_species` へ戻し、`discovery_name` とtaxonomic sourceを保持する。

## 2. 既取得本文から4形質を追加抽出

対象形質はTRYの公式trait catalogにも独立項目として存在する（flower size 3568、floral depth 3574、nectar tube depth 3579、floral reward 205、inflorescence type 2934）。[TRY trait selection](https://www.try-db.org/TryWeb/Prop023.php)

`config/trait_ontology.yml` は既に `tube_depth_class`、`flower_size_class`、`inflorescence_display`、`reward_type` を定義している。最小実装は次のとおり。

1. `config/global_floral_access_keywords.yml` にinflorescence displayとrewardの保守的な語彙を追加し、既存のtube/size語彙には数値・単位・範囲規則を足す。
2. `src/island_v2/trait_web_reported_scout.py` で4形質すべてに花器官近傍、否定、曖昧表現のfail-closed判定を適用する。
3. `source_texts.csv` の固定済み本文を再解析し、legacy nine-columnを広げず `v2_reported_candidates.csv` のlong形式へ出す。
4. source URL、本文hash、excerpt、matched ruleを維持し、unreviewed候補をacceptedへ直結しない。

`tube_depth_class`、`flower_size_class`、`inflorescence_display` は花構造軸へ入る。`reward_type` は有用な補助形質だが、3軸への対応を別途明示しない限りcoverage分子へ加えない。

## 3. Europe PMC OA全文・supplementary

現行の `src/island_v2/europe_pmc_discovery.py` はexact speciesを含むtitle/abstractとPMCID、DOI、OA状態、licenseを既に保存している。OA候補だけを `/{PMCID}/fullTextXML` へ進め、表・caption・methods/resultsを抽出する。補遺は `/{PMCID}/supplementaryFiles` のZIPを取得し、NXML、CSV、TSV、XLSXを優先する。[Europe PMC REST API](https://europepmc.org/RestfulWebService)

Europe PMCのOA subsetはAPI、OAI、FTPによる自動取得を認める一方、記事ごとのライセンスは同一ではない。[Open access subset](https://europepmc.org/downloads/openaccess) main websiteのbulk crawlは禁止され、公式サービス以外を使ってはならない。[Europe PMC developer resources](https://europepmc.org/developers)

実装は `europe_pmc_discovery.py` にOA本文取得adapterを加え、`src/island_v2/v1_category_search.py::_europe_pmc_sources` へ接続する。lineageにはPMCID、DOI、article license、supplement filename/hash、table/row座標を保存する。accepted nameまたは検証済みsynonymが対象表・caption・近傍本文に現れない行はspecies-directにしない。HTML/PDF fallbackは設けず、非OAとunsupported supplementはaccounted missにする。

## 4. Botanical Latin・多言語抽出

WFOのcontent schemaはdescriptionごとに `taxonID`、`type`、`language`、`description`、`source`、`rights`、`license` を分離するため、言語別抽出とlineage保存に適している。[WFO Content Data Contributor Guidelines](https://about.worldfloraonline.org/images/uploads/documents/WFOGuidelinesforContentDataContributorsV._2.04.pdf)

現状は `config/global_trait_campaign.yml` が英語に続けてes/fr/de/pt/itを検索し、`src/island_v2/multilingual_wikimedia_recovery.py` がja/zh/ruの繁殖語彙を持つ。次は取得providerを増やす前に、同じ固定本文へ次を行う。

- language tagをsource recordへ必須化し、言語別の花器官・色・形・繁殖語彙を分離する。
- botanical Latinは `flos/floris`、`corolla`、色・形容詞の性数格変化をlemma相当の規則で扱う。
- 否定、比較、複数器官、果実・葉色との混同をgold setで評価する。
- 言語未判定、Latin/現代語混在、低精度ruleはunreviewedのままにする。

一般言語モデルだけでbotanical Latinの正規化精度は保証されない。必要ならStanzaのlanguage identificationとLatinを含むモデルを前処理に使えるが、最終値はrepository固有の辞書・テストで検証する。[Stanza language identification](https://stanfordnlp.github.io/stanza/langid.html) [Stanza model performance](https://stanfordnlp.github.io/stanza/performance.html)

## 5. TRY対象形質export

TRYは登録済みPIによるdata requestでtraitsと任意のspeciesを指定する。public-onlyは処理後に提供され、restricted recordを含むrequestはdata ownerの許可を待つ。[TRY request portal](https://www.try-db.org/TryWeb/Prop0.php) 利用はTRYのIntellectual Property Guidelinesに従い、raw dataの再配布可否をrequest単位で確認する。[TRY IP guidelines](https://www.try-db.org/TryWeb/TRY_Intellectual_Property_Guidelines.pdf) [TRY Data Use Policy](https://www.try-db.org/TryWeb/DataUsePolicy.php)

申請対象はflower colour、corolla type/symmetry/size/depth、self-incompatibility、autogamy/xenogamy、sexual syndrome、reward、inflorescence typeに限定する。受領後のimporterはTRY release formatの `AccSpeciesName`、`OrigObsDataID`、`Reference`、trait ID、unit、dataset metadataを保持する。[TRY release notes](https://www.try-db.org/TryWeb/TRY_Data_Release_Notes.pdf) `OrigObsDataID` でTRY内重複を除き、`Reference` をsource lineageの基本単位とする。Actionsには認証情報やraw restricted exportを置かず、許可済み入力のhashとderived candidate artifactだけをmanifestへ記録する。

## 6. WFO description全件取得

WFO portalは多数のdescriptionを保持するが、個々のcontent elementに別ライセンスがあり得るため、record単位でsource、rights、licenseを保持する必要がある。[WFO portal](https://www.worldfloraonline.org/) [WFO Terms of Use](https://www.worldfloraonline.org/termsOfUse)

現行 `src/island_v2/v1_category_search.py::_world_flora_sources` はsearch/taxon HTMLを読む実装だが、`https://www.worldfloraonline.org/robots.txt` は一般botを `Disallow: /` としており、`docs/public_web_trait_shard_campaign.md` でもproviderは無効化されている。このcollectorをActionsで有効化してはならない。

実装開始条件は、WFOからdescription extensionを含む公式DwC-A exportを得るか、WFOが明示する元content providerのDwC-Aを取得できることである。受領後はWFO Content Guidelinesの `taxonID`、`type`、`language`、`description`、`source`、`rights`、`license` をそのままsnapshotし、WFO IDからaccepted speciesへ結ぶ。公式exportが得られない間は、既存eFlorasと許可済みprovider artifactの拡張を優先する。

## Delivery order and gates

1. WFO Plant List snapshotからalias artifactを作り、未解決種だけで既存providerの回収率を再測定する。
2. 固定済みsource textへ4形質抽出を追加し、hand-reviewed gold setでprecision gateを通す。
3. Europe PMC OA全文・supplementをtargeted queueとして実装する。
4. 多言語・botanical Latinを言語別gold setで段階的に有効化する。
5. TRY requestは並行して人手申請し、受領・利用条件確認後にimportする。
6. WFO全件は公式content exportを得るまで着手しない。

どのlaneもfamily inference、global fallback、pending候補の自動accepted化を行わない。評価は固定した未解決 `accepted_species × axis` を使い、source lineage重複除去後の純増で比較する。

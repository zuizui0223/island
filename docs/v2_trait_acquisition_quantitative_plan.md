# v2 trait acquisition and quantitative outcomes

## Decision

v2 does not repeat the v1 workflow of immediately coding every species as 0/1.
Acquisition stores source-backed categorical observations in long format.
Categorical composition is the broad-coverage primary analysis; binary
contrasts remain secondary summaries for continuity. Numeric measurements are
optional enrichment and never determine whether a species enters the core set.

The broad-coverage minimum is flower colour plus flower form (symmetry is the
shape fallback). Reproductive assurance is the next priority. The exact
15-trait order is frozen in `config/v2_trait_acquisition.yml`. For each trait,
species-direct evidence is attempted first, then only unresolved cells may use
a qualified genus inference, then a qualified family inference; unsupported
cells remain explicit unresolved values. Direct pollen-vector evidence is
collected whenever available, but effective pollinator guild is not an
exhaustive completion target and cannot block acquisition of intrinsic floral
or reproductive traits.

The machine-readable contract is
[`config/v2_trait_acquisition.yml`](../config/v2_trait_acquisition.yml).

## Acquisition layers

1. **Raw evidence**: scientific name, original description or measurement,
   source URL/citation, evidence scope, wild/cultivated status and review state.
2. **Harmonized traits**: v2 ontology values plus normalized continuous
   measurements. Multiple observations remain multiple rows.
3. **Likely proxies**: morphology-based or taxonomic inference, always marked
   `inference_status=likely`, kept outside direct evidence and outside the main
   analysis.
4. **Analysis features**: generated reproducibly from harmonized traits. These
   can change without recollecting the source evidence.

## Quantitative axes

### Visual signal

Human colour labels provide a broad categorical composition outcome; they do
not provide a pollinator's colour contrast. Bee, fly or bird visual models need
flower reflectance by wavelength, an illuminant, a background spectrum and a
declared receptor model. Receptor contrast is therefore a continuous outcome
for the spectral subset only. It is never back-filled from `red`, `yellow` or
other human labels, and spectral collection is not a core acquisition target.

### Floral access and specialization

Retain floral form and symmetry for every possible species. The main analysis
uses the form/symmetry categories and may also show a simple open/generalized
versus restrictive/specialized sensitivity contrast. Flower diameter, tube
depth, aperture width and display diameter are retained when found, but no
extra web effort is spent solely to obtain them. A continuous access descriptor
is attempted only if measurement coverage proves adequate.

### Reproductive assurance

Keep self-incompatibility, autonomous selfing, mating system and cleistogamy as
separate observations. Quantitative autofertility, outcrossing rate and seed-set
ratios should be retained when sources report them. Self-compatibility alone is
not autonomous selfing and cannot fill that field.

### Functional replacement

Visitor, pollen-contact and seed-set evidence stays long by pollinator guild.
A red tubular phenotype may generate a `likely` bird-or-butterfly phenotype
proxy, but not a confirmed guild. Actual replacement requires accepted evidence
for effective guilds. Island-level Bombus state is never an input to the floral
proxy, so the later Bombus-by-alternative-channel test is not circular.

## Scaling to 100,000 species

Run free bulk downloads and scientific-name joins first. Extract cached flora,
monograph, encyclopedia and institutional horticulture text second. Query public
web pages only for remaining colour/form gaps. Proxy generation and measurement
normalization are local post-processing passes and make no network requests.

Coverage reports must show separate denominators for colour, form, continuous
size/access measurements, reproductive assurance, direct pollinator evidence
and likely proxies. Missing direct ecology is unresolved, not biological
absence.

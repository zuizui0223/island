# Chapter 1 response-compression audit — 2026-09-15

Status: **completed diagnostic on the all-observed beta-binomial route**.

## Provenance

- workflow run: `34964112334`
- artifact: `chapter1-response-compression-34964112334`
- artifact ID: `10394930700`
- digest: `sha256:960802ffc47a2287849693fc37788466c6778beae6388823453bf39eb51782fc`

## Question

The broad six-atomic North--Tropical response is supported, whereas the historical species-level two-axis syndrome-concordance route is weak on the all-observed flora. Is that discrepancy caused simply by reducing six outcomes to two dimensions?

To answer this without changing the fitted beta-binomial models, the frozen six-component North--Tropical interaction vector and its full spatial-block sandwich covariance were projected onto:

1. the full six-axis vector;
2. two equal-weight family means: accessibility/generalisation and reproductive assurance;
3. four within-family contrasts discarded by a two-family mean compression.

## Result

### All-analysis evidence

- full six-axis vector: chi-square `23.1011`, df=6, **p=0.0007633**;
- two family means: chi-square `12.8179`, df=2, **p=0.0016467**;
- within-family contrasts discarded by the two-family projection: chi-square `3.3028`, df=4, **p=0.5085**;
- Euclidean signal retained by the two-family subspace: **78.22%**.

### High/Medium direct evidence

- full six-axis vector: chi-square `18.8072`, df=6, **p=0.0045020**;
- two family means: chi-square `8.8304`, df=2, **p=0.0120921**;
- within-family contrasts discarded by the two-family projection: chi-square `6.8548`, df=4, **p=0.1438**;
- Euclidean signal retained by the two-family subspace: **64.82%**.

## Interpretation

The previous explanation "the two-axis compression destroys the signal" is too strong and is rejected by this audit.

A simple linear two-family projection of the six atomic interaction slopes remains supported in both evidence scopes. The structure orthogonal to those two family means is not itself supported. Therefore the discrepancy with the historical `generalized_accessible + selfing_core` route is not caused by dimensionality alone.

The remaining difference lies in **response construction and support weighting**. The historical syndrome-concordance route:

- first constructs species-level multivariate concordance scores;
- requires at least two informative traits per species;
- ignores missing traits in the within-species denominator;
- then averages those species scores within islands;
- finally fits a continuous island-score model.

The atomic probability route instead uses every species resolved for each individual outcome and retains each outcome-specific denominator explicitly.

## Next gate

Diagnose the syndrome-score construction directly by rerunning the same two historical axes under controlled support variants, especially:

1. `minimum_informative_traits = 1` versus the frozen value 2;
2. same original trait weights;
3. no change to evidence scope, island covariates, contexts or multiple-testing family.

This will determine whether the weak historical two-axis result is primarily a multi-trait support/intersection effect rather than a genuine absence of a two-family response.

<!-- README.md is generated from README.Rmd. Please edit that file -->

# exchanger: Bayesian Entity Resolution with Exchangeable Random Partition Priors

## Overview

This R package implements the Bayesian entity resolution model described
in (Marchant, Rubinstein, and Steorts 2020). The model assumes records
are generated from a latent population of entities via a distortion
process. Entity resolution is performed by inferring the latent *linkage
structure* using a Markov chain Monte Carlo algorithm.

## Features

  - Supports user-defined attribute distance functions for modeling the
    distortion process.
  - Supports missing values.
  - Implements the Ewens-Pitman (EP) family of exchangeable clustering
    priors (Pitman 2006), with support for hyperpriors on the EP
    parameters. This includes the Ewens partition (related to the
    Dirichlet process), Pitman-Yor partition (related to the Pitman-Yor
    process) and Generalized Coupon partition (related to finite mixture
    models).
  - Inference is implemented in C++ using partially-collapsed Gibbs
    sampling (Dyk and Park 2008). Several computational tricks are used
    to improve scalability.

## Installation

`exchanger` is not currently available on CRAN. The latest development
version can be installed from source using `devtools` as follows:

``` r
library(devtools)
install_github("cleanzr/exchanger")
```

## Example

We demonstrate how to perform entity resolution using the `RLdata500`
test data set included with the package. `RLdata500` contains 500
synthetic personal records, 50 of which are duplicates with
randomly-generated errors. In the code block, below we load the package
and examine the first few rows of `RLdata500`.

``` r
library(exchanger)
library(comparator)     # provides string comparison functions
library(clevr)          # provides functions for supervised evaluation
head(RLdata500)
#>   fname_c1 fname_c2 lname_c1 lname_c2   by bm bd
#> 1  CARSTEN     <NA>    MEIER     <NA> 1949  7 22
#> 2     GERD     <NA>    BAUER     <NA> 1968  7 27
#> 3   ROBERT     <NA> HARTMANN     <NA> 1930  4 30
#> 4   STEFAN     <NA>    WOLFF     <NA> 1957  9  2
#> 5     RALF     <NA>  KRUEGER     <NA> 1966  1 13
#> 6  JUERGEN     <NA>   FRANKE     <NA> 1929  7  4
```

Next we specify the model parameters for the entity attributes. For
simplicity, we use a uniform prior for the distortion probability
associated with each attribute.

``` r
unif_prior <- BetaRV(1, 1)
```

We model the distortion for the name attributes (`fname_c1` and
`lname_c1`) using Levenshtein (edit) distance. The `transform_dist_fn`
applies a threshold at a distance of 3, so that any pair of names `y`,
`x` with an edit distance strictly greater than 3 are assumed to be
“infinitely” apart. This improves computational efficiency, as the
inference algorithm exploits the fact that `x` can never be a distortion
of `y` if their distance is strictly greater than 3. For the date
attributes (`by`, `bm` and `bd`) we assume a constant distance function,
which means all values in the domain are equally likely a priori as a
distorted alternative. The `CategoricalAttribute` constructor provides a
shortcut for an attribute with a constant distance function.

``` r
attr_params <- list(
  fname_c1 = Attribute(transform_dist_fn(Levenshtein(), threshold=3), 
                         distort_prob_prior = unif_prior),
  lname_c1 = Attribute(transform_dist_fn(Levenshtein(), threshold=3), 
                         distort_prob_prior = unif_prior),
  by = CategoricalAttribute(distort_prob_prior = unif_prior),
  bm = CategoricalAttribute(distort_prob_prior = unif_prior),
  bd = CategoricalAttribute(distort_prob_prior = unif_prior)
)
```

Finally we specify the prior over the linkage structure (clustering).
Here we use a Pitman-Yor random partition for the prior, with
hyperpriors on the concentration and discount parameters.

``` r
clust_prior <- PitmanYorRP(alpha = GammaRV(1, 1), d = BetaRV(1, 1))
```

All that remains is to initialize the model and run inference. For
demonstration purposes, we only generate 100 posterior samples of the
linkage structure.

``` r
model <- exchanger(RLdata500, attr_params, clust_prior)
result <- run_inference(model, n_samples=100, thin_interval=10, burnin_interval=1000)
```

We can obtain a point estimate of the most likely linkage structure
using the *shared most probable maximal matching sets method* (Steorts,
Hall, and Fienberg 2016).

``` r
pred_clust <- smp_clusters(result)
```

Then we evaluate with respect to the true matching.

``` r
n_records <- nrow(RLdata500)
true_pairs <- membership_to_pairs(identity.RLdata500)
pred_pairs <- clusters_to_pairs(pred_clust)
measures <- eval_report_pairs(true_pairs, pred_pairs, num_pairs=n_records*(n_records-1)/2)
print(measures)
#> $precision
#> [1] 1
#> 
#> $recall
#> [1] 0.98
#> 
#> $specificity
#> [1] 1
#> 
#> $sensitivity
#> [1] 0.98
#> 
#> $f1score
#> [1] 0.989899
#> 
#> $accuracy
#> [1] 0.999992
#> 
#> $balanced_accuracy
#> [1] 0.99
```

## License

GPL-3

## References

<div id="refs" class="references hanging-indent">

<div id="ref-dyk_partially_2008">

Dyk, David A. van, and Taeyoung Park. 2008. “Partially Collapsed Gibbs
Samplers.” *Journal of the American Statistical Association* 103 (482):
790–96. <https://doi.org/10.1198/016214508000000409>.

</div>

<div id="ref-marchant_exchangeable_2020">

Marchant, Neil G., Benjamin I. P. Rubinstein, and Rebecca C. Steorts.
2020. “Exchangeable Clustering Priors for Bayesian Entity Resolution.”

</div>

<div id="ref-pitman_exchangeable_2006">

Pitman, Jim. 2006. “Exchangeable Random Partitions.” In *Combinatorial
Stochastic Processes: Ecole d’Eté de Probabilités de Saint-Flour XXXII –
2002*, edited by Jean Picard, 37–53. Berlin, Heidelberg: Springer Berlin
Heidelberg. <https://doi.org/10.1007/3-540-34266-4_3>.

</div>

<div id="ref-steorts2016">

Steorts, Rebecca C., Rob Hall, and Stephen E. Fienberg. 2016. “A
Bayesian Approach to Graphical Record Linkage and Deduplication.”
*Journal of the American Statistical Association* 111 (516): 1660–72.
<https://doi.org/10.1080/01621459.2015.1105807>.

</div>

</div>

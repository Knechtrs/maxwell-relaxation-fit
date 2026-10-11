# maxwellfit

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23298133.svg)](https://doi.org/10.5281/zenodo.23298133)

Robust fitting of generalized Maxwell models to stress-relaxation curves, with checks that tell you when a fit should not be trusted.

I wrote this for the analysis of viscoelastic stress relaxation in human fracture hematoma and engineered alginate hydrogels in Knecht et al., *Distinct matrix viscoelasticity in bone fracture hematoma steers macrophage phenotype* (*Nature Communications* under revision). The normalized stress is modelled as

```
sigma(t) = 1 - sum_i A_i * (1 - exp(-t / tau_i))
```

Multi-exponential fits are ill-conditioned: different starting values can converge to different solutions, time constants can swap labels between samples, and information criteria can reward elements that the data do not support. The package addresses each of these explicitly.

## What it does

| Problem | How it is handled |
| --- | --- |
| Local minima | 20 random starts per fit (log-uniform taus across the recorded window); the lowest-RSS solution is reported, and `starts_at_best` shows how many starts reached it |
| Taus spanning orders of magnitude | taus are fitted on the log scale with bounds; amplitudes bounded to [0, 1] |
| Label switching between samples | taus are always returned sorted, so `tau1` is the fastest process in every sample |
| Over-complex models | `compare_maxwell()` reports AIC, BIC and a nested F-test, plus residual lag-1 autocorrelation as a warning sign |
| Unphysical or unidentifiable fits | flags when `sum(A) > 1`, a tau sits at its bound, or a tau exceeds the recording time |

## Usage

```r
# install.packages("remotes")
remotes::install_github("Knechtrs/maxwell-relaxation-fit@v0.1.0")  # version used in the paper
library(maxwellfit)

fit <- fit_maxwell(time, stress_norm, n_elements = 2, seed = 1)
fit
#> 2-element Maxwell fit (n = 6000 points)
#>   A1 = 0.3002   tau1 = 20.09 s
#>   A2 = 0.4999   tau2 = 800.56 s
#>   RMSE = 0.003056, R2 = 0.99948, AIC = -52451.5, BIC = -52418.0
#>   starts reaching the best solution: 18 of 19 converged

compare_maxwell(time, stress_norm, max_elements = 3, seed = 1)

# Many samples in long format (one row per time point)
fit_maxwell_by_sample(df, id = "sample", time = "time", stress = "load_norm",
                      t_max = 3000, n_elements = 2, seed = 1)
```

Preprocessing (Butterworth filtering, normalization to peak load, defining t = 0 at peak load) is done upstream and is not part of this package; the helpers used for the paper are in [ktools](https://github.com/Knechtrs/ktools).

## Use in the manuscript

This package (version 0.1.0) is the fitting code for the article. For the one- and two-element models it reproduces the time constants obtained with the original analysis scripts to within 0.01% for every hematoma and hydrogel sample, with identical residual sums of squares. Models with three or more elements were frequently not identifiable (total relaxation > 1, time constants beyond the recording window, different starts converging to different solutions); `compare_maxwell()` reports these cases as flags.

## Tests

```r
testthat::test_dir("tests/testthat")
```

The tests check parameter recovery for one- and two-element models across noise realisations, consistent tau ordering, that the reported parameters belong to the lowest-RSS start, and that BIC neither under- nor over-fits when noise is independent.

## Citation

If you use this package, please cite the article above and the archived software release:

> Knecht, R. S. (2026). maxwellfit: Robust Generalized Maxwell Fitting of Stress-Relaxation Curves (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.23298133

## License

MIT

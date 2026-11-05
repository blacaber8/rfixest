# rfixest: fast fixed-effects regressions in Stata, powered by R's fixest

`rfixest` is a Stata command that estimates linear, Poisson and negative binomial regressions with many
high-dimensional fixed effects. It hands the data to the R package
[fixest](https://lrberge.github.io/fixest/), one of the fastest implementations available, and brings the
results back to Stata as a regular estimation command. You never have to open R.

- Same syntax as `reghdfe`: `rfixest y x1 x2, absorb(id year) vce(cluster id)`
- Same point estimates and standard errors as `reghdfe` (OLS) and `ppmlhdfe` (Poisson)
- Several outcomes in one call, `rfixest (y1 y2 y3) x, absorb(id year)`, with the fixed effects computed once
- Interacted fixed effects (`absorb(id industry#year)`), factor variables, `aweight`/`pweight`,
  robust and multiway clustered standard errors
- Saves the estimated fixed effects and the residuals, as in `reghdfe`: `absorb(FE_id=id year) residuals(res)`
- Works with `estimates store`, `esttab`, `test`, `lincom` and `e(sample)`
- R runs in the background: no R window, no R code to write

## Installation

1. Install R (any recent version) from [cran.r-project.org](https://cran.r-project.org) if you do not have it.
2. In Stata, type:

```stata
net install rfixest, from("https://raw.githubusercontent.com/blacaber8/rfixest/main/") replace
rfixest setup
```

`rfixest setup` finds R on your computer, installs the R package fixest if needed, and remembers where R is.
You only do this once. If R is in an unusual folder, tell it where `Rscript` is:

```stata
rfixest setup, rscript("D:/Programs/R/R-4.4.1/bin/x64/Rscript.exe")
```

## Examples

```stata
webuse nlswork, clear

* individual and year fixed effects, clustered by individual
rfixest ln_wage tenure ttl_exp, absorb(idcode year) vce(cluster idcode)

* industry-by-year fixed effects and a factor variable
rfixest ln_wage tenure ttl_exp i.union, absorb(idcode ind_code#year) cluster(idcode)

* two outcomes at once, stored as m1 and m2
rfixest (ln_wage hours) tenure ttl_exp, absorb(idcode year) cluster(idcode) estore(m)
esttab m1 m2, se

* save the fixed effects and the residuals
rfixest ln_wage tenure ttl_exp, absorb(FE_id=idcode FE_year=year) residuals(res) cluster(idcode)

* Poisson
rfixest hours tenure ttl_exp, absorb(idcode year) cluster(idcode) family(poisson)
```

See `help rfixest` and the guide [rfixest_guide.pdf](rfixest_guide.pdf) for all options and worked examples.

## Speed

Three outcomes, 2 million observations, three fixed effects (100,000 individuals, 20,000 firms, 20 years),
standard errors clustered by individual, on a Windows desktop:

| | Time |
|---|---|
| `reghdfe`, one call per outcome | 16.3 s |
| `rfixest (y1 y2 y3) ...`, one call | 3.9 s |

## Requirements

- Stata 16 or later (the command uses Stata data frames, introduced in version 16)
- R 4.0 or later, with the R package `fixest` (installed by `rfixest setup`)
- On Windows, Stata must run with its graphical interface (as it normally does). On macOS and Linux any mode works.

## Credits

All the estimation is done by [fixest](https://lrberge.github.io/fixest/) (Laurent Bergé). Please cite it:
Bergé, L. (2018), "Efficient estimation of maximum likelihood models with multiple fixed-effects: the R package FENmlm", CREA Discussion Papers 13.

The R code that `rfixest` runs is embedded in `rfixest.ado`, so installation needs nothing else. `rfixest.R` in this repository is a readable copy of that code.

License: MIT.

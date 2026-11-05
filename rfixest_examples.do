* Examples of rfixest with the NLS Women dataset that ships with Stata
webuse nlswork, clear

* 1. Individual and year fixed effects, standard errors clustered by individual
rfixest ln_wage tenure ttl_exp, absorb(idcode year) vce(cluster idcode)

* 2. Same model with reghdfe, for comparison
reghdfe ln_wage tenure ttl_exp, absorb(idcode year) vce(cluster idcode)

* 3. Interacted fixed effects (industry by year) and a factor variable (union membership)
rfixest ln_wage tenure ttl_exp i.union, absorb(idcode ind_code#year) cluster(idcode)

* 4. Two-way clustering
rfixest ln_wage tenure ttl_exp, absorb(idcode year) cluster(idcode ind_code)

* 5. Several outcomes at once, stored as m1 and m2, and a table
rfixest (ln_wage hours) tenure ttl_exp, absorb(idcode year) cluster(idcode) estore(m) quietly
esttab m1 m2, se star(* 0.10 ** 0.05 *** 0.01) b(%9.4f) se(%9.4f) stats(N r2_within, fmt(%9.0fc %9.4f) labels("Observations" "Within R-squared"))

* 6. Save the fixed effects and the residuals
rfixest ln_wage tenure ttl_exp, absorb(FE_id=idcode FE_year=year) residuals(res) cluster(idcode) quietly
summarize FE_id FE_year res

* 7. Poisson regression
rfixest hours tenure ttl_exp, absorb(idcode year) cluster(idcode) family(poisson)

* 8. Postestimation
rfixest ln_wage tenure ttl_exp, absorb(idcode year) cluster(idcode) quietly
test tenure = ttl_exp
lincom tenure + ttl_exp
count if e(sample)

{smcl}
{* *! version 1.1.0  04nov2026}{...}
{title:Title}

{phang}
{bf:rfixest} {hline 2} Fast fixed-effects regressions in Stata, estimated by the R package {bf:fixest}


{title:Syntax}

{p 8 15 2}
{cmd:rfixest} {depvar} [{indepvars}] {ifin} {weight} [{cmd:,} {it:options}]

{p 8 15 2}
{cmd:rfixest} {cmd:(}{it:depvar1 depvar2 ...}{cmd:)} [{indepvars}] {ifin} {weight} [{cmd:,} {it:options}]

{p 8 15 2}
{cmd:rfixest setup} [{cmd:,} {opt rscript(path)}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{opt a:bsorb(terms)}}fixed effects to absorb, for example {cmd:absorb(id year)} or {cmd:absorb(id year#industry)}; {cmd:absorb(FE_city=city)} also saves the estimated fixed effect in the new variable {cmd:FE_city}{p_end}
{synopt:{opt f:amily(name)}}{cmd:ols} (default), {cmd:poisson} or {cmd:negbin}{p_end}
{synopt:{opt keepsing:letons}}keep singleton groups (dropped by default, as in {cmd:reghdfe}){p_end}

{syntab:Standard errors}
{synopt:{opt vce(vcetype)}}{cmd:unadjusted}, {cmd:robust} or {cmd:cluster} {it:varlist} (two or more variables give multiway clustering){p_end}
{synopt:{opt cl:uster(varlist)}}same as {cmd:vce(cluster} {it:varlist}{cmd:)}{p_end}
{synopt:{opt dfres:id}}t and p values with N - K residual degrees of freedom instead of G - 1 when clustering{p_end}
{synopt:{opt ppmlse}}no small-sample adjustment of the variance, as {cmd:ppmlhdfe}{p_end}

{syntab:Reporting and storing}
{synopt:{opt res:iduals(newvar)}}save the residuals in {it:newvar} (missing outside the estimation sample){p_end}
{synopt:{opt est:ore(name)}}store the estimates; with several outcomes, {it:name}1, {it:name}2, ...{p_end}
{synopt:{opt esta:rt(#)}}first number used to name stored estimates; default 1{p_end}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}
{synopt:{opt qui:etly}}do not display the table{p_end}
{synopt:{opt full:sample}}{cmd:e(sample)} marks every observation passed to R, including dropped singletons{p_end}
{synopt:{opt time:it}}report the time spent in Stata and in R{p_end}

{syntab:Technical}
{synopt:{opt thr:eads(#)}}number of threads used by fixest; default is fixest's own default{p_end}
{synopt:{opt tol(#)}}convergence tolerance of the fixed effects; default {cmd:1e-6}{p_end}
{synopt:{opt rscript(path)}}path to {cmd:Rscript}, if it is not found automatically{p_end}
{synopt:{opt tmp:dir(path)}}folder for the temporary files; default is Stata's temporary folder{p_end}
{synopt:{opt keep:files}}keep the temporary files (for debugging){p_end}
{synoptline}
{p 4 6 2}{it:indepvars} may contain factor variables and time-series operators; see {help fvvarlist}.{p_end}
{p 4 6 2}{cmd:aweight}s and {cmd:pweight}s are allowed; {cmd:pweight}s imply {cmd:vce(robust)} unless {cmd:vce()} is given.{p_end}


{title:Description}

{pstd}
{cmd:rfixest} estimates linear, Poisson and negative binomial regressions with any number of high-dimensional fixed effects.
It sends the data to the R package {browse "https://lrberge.github.io/fixest/":fixest} (Berge 2018), which is among the fastest
implementations available, and brings the results back to Stata as a regular estimation command: {cmd:e(b)}, {cmd:e(V)},
{cmd:e(sample)}, so {cmd:estimates store}, {cmd:esttab}, {cmd:test} and {cmd:lincom} work as usual.
{cmd:predict} and {cmd:margins} are not supported.

{pstd}
Point estimates and standard errors match {cmd:reghdfe} (OLS) and {cmd:ppmlhdfe} (Poisson) with the same options.

{pstd}
Several outcomes can be estimated at once with the same regressors, fixed effects and sample:
{cmd:rfixest (y1 y2 y3) x, absorb(id year)}. R reads the data once and, when the outcomes are missing in the same
observations, computes the fixed effects once, which is much faster than separate regressions. Each outcome uses its own
nonmissing observations.


{title:Requirements and first use}

{pstd}
{cmd:rfixest} needs Stata 16 or later (it uses data frames) and R ({browse "https://cran.r-project.org":cran.r-project.org}), version 4.0 or later. The first time, type

{phang2}{cmd:. rfixest setup}{p_end}

{pstd}
This finds R on your computer, installs the R package fixest if needed, and remembers the location of R for future sessions.
If R is installed in an unusual place, give its location:

{phang2}{cmd:. rfixest setup, rscript("D:/Programs/R/R-4.4.1/bin/x64/Rscript.exe")}{p_end}

{pstd}
R runs in the background without opening any window. Stata must run with its graphical interface on Windows
(not in batch mode); on macOS and Linux both modes work.


{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}

{pstd}Individual and year fixed effects, standard errors clustered by individual{p_end}
{phang2}{cmd:. rfixest ln_wage tenure ttl_exp, absorb(idcode year) vce(cluster idcode)}{p_end}

{pstd}Interacted fixed effects (industry-by-year) and a factor variable{p_end}
{phang2}{cmd:. rfixest ln_wage tenure ttl_exp i.union, absorb(idcode ind_code#year) cluster(idcode)}{p_end}

{pstd}Two-way clustering{p_end}
{phang2}{cmd:. rfixest ln_wage tenure ttl_exp, absorb(idcode year) cluster(idcode ind_code)}{p_end}

{pstd}Several outcomes at once, then a table with {cmd:esttab}{p_end}
{phang2}{cmd:. rfixest (ln_wage hours) tenure ttl_exp, absorb(idcode year) cluster(idcode) estore(m)}{p_end}
{phang2}{cmd:. esttab m1 m2, se star(* 0.10 ** 0.05 *** 0.01)}{p_end}

{pstd}Save the fixed effects and the residuals{p_end}
{phang2}{cmd:. rfixest ln_wage tenure ttl_exp, absorb(FE_id=idcode FE_year=year) residuals(res) cluster(idcode)}{p_end}

{pstd}Poisson regression{p_end}
{phang2}{cmd:. rfixest hours tenure ttl_exp, absorb(idcode year) cluster(idcode) family(poisson)}{p_end}


{title:Saved fixed effects and residuals}

{pstd}
{cmd:absorb(}{it:newvar}{cmd:=}{it:term}{cmd:)} saves the estimated fixed effect of {it:term} for each observation, and
{cmd:residuals(}{it:newvar}{cmd:)} saves the residuals, both only for observations in the estimation sample. They need a single
dependent variable. With two or more fixed effects, each one is identified only up to a constant, so the level of a saved fixed
effect can differ from {cmd:reghdfe} by a constant; differences between groups, the sum of all fixed effects plus the constant
and the residuals are the same. For Poisson and negative binomial models the fixed effects are on the scale of the linear index
and the residuals are {it:y} minus the predicted mean.


{title:Stored results}

{pstd}{cmd:rfixest} stores the following in {cmd:e()}:{p_end}
{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations used{p_end}
{synopt:{cmd:e(N_full)}}observations before dropping singletons{p_end}
{synopt:{cmd:e(num_singletons)}}observations dropped as singletons or perfectly predicted{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters (smallest dimension with multiway clustering){p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom used for t tests{p_end}
{synopt:{cmd:e(df_a)}}degrees of freedom absorbed by the fixed effects{p_end}
{synopt:{cmd:e(r2)}, {cmd:e(r2_a)}, {cmd:e(r2_within)}, {cmd:e(r2_a_within)}, {cmd:e(rmse)}}fit statistics (OLS){p_end}
{synopt:{cmd:e(ll)}, {cmd:e(r2_p)}}log likelihood and pseudo R-squared (Poisson, negative binomial){p_end}
{synopt:{cmd:e(time_r)}}seconds spent in R{p_end}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:rfixest}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}dependent variable{p_end}
{synopt:{cmd:e(absvars)}}absorbed fixed effects{p_end}
{synopt:{cmd:e(family)}}{cmd:ols}, {cmd:poisson} or {cmd:negbin}{p_end}
{synopt:{cmd:e(vce)}, {cmd:e(clustvar)}}type of standard errors and cluster variables{p_end}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}coefficients and variance matrix{p_end}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}estimation sample{p_end}


{title:Reference}

{phang}
Berge, L. 2018. Efficient estimation of maximum likelihood models with multiple fixed-effects: the R package FENmlm.
CREA Discussion Papers 13.


{title:Author}

{pstd}
{browse "https://github.com/blacaber8/rfixest":github.com/blacaber8/rfixest}. Please report problems there.
{p_end}

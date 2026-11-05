*! rfixest_p 1.2.0  05nov2026
*! predict after rfixest
program define rfixest_p
	version 16
	syntax newvarname [if] [in] [, xb xbd d Residuals mu]
	local opt "`xb'`xbd'`d'`residuals'`mu'"
	if `: word count `xb' `xbd' `d' `residuals' `mu'' > 1 {
		di as error "only one of xb, xbd, d, residuals and mu is allowed"
		exit 198
	}
	if "`opt'" == "" {
		local opt "xb"
		di as text "(option xb assumed; linear prediction without the fixed effects)"
	}
	if "`opt'" == "mu" & "`e(family)'" == "ols" {
		di as error "mu is available after family(poisson) and family(negbin); use xbd after OLS"
		exit 198
	}
	marksample touse, novarlist

	if "`opt'" == "xb" {
		_predict double `varlist' if `touse', xb
		exit
	}

	local need = e(n_fe) > 0
	if `need' & e(fe_all) != 1 {
		di as error "predict, `opt' needs every fixed effect saved at estimation:"
		di as error "add the option {bf:savefe}, or name each one, for example {bf:absorb(FE1=id FE2=year)}"
		exit 198
	}
	tempvar xbv dv
	qui _predict double `xbv' if `touse', xb
	qui gen double `dv' = 0 if `touse'
	foreach v in `e(fe_saved)' {
		cap confirm numeric variable `v'
		if _rc {
			di as error "saved fixed effect `v' not found; estimate the model again"
			exit 111
		}
		qui replace `dv' = `dv' + `v' if `touse'
	}

	if "`opt'" == "d" {
		qui gen double `varlist' = `dv' if `touse'
		label variable `varlist' "Sum of the fixed effects"
	}
	else if "`opt'" == "xbd" {
		qui gen double `varlist' = `xbv' + `dv' if `touse'
		label variable `varlist' "Linear prediction including the fixed effects"
	}
	else if "`opt'" == "mu" {
		qui gen double `varlist' = exp(`xbv' + `dv') if `touse'
		label variable `varlist' "Predicted mean"
	}
	else {
		if "`e(family)'" == "ols" qui gen double `varlist' = `e(depvar)' - (`xbv' + `dv') if `touse'
		else qui gen double `varlist' = `e(depvar)' - exp(`xbv' + `dv') if `touse'
		label variable `varlist' "Residuals"
	}
end

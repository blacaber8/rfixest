*! rfixest 1.2.0  05nov2026
*! Fixed-effects regressions in Stata estimated by the R package fixest
program define rfixest, eclass
	version 16
	if replay() {
		if "`e(cmd)'" != "rfixest" error 301
		_rfx_show `0'
		exit
	}

	gettoken sub rest : 0, parse(" ,")
	if "`sub'" == "setup" {
		_rfx_setup `rest'
		exit
	}

	local cmdline `"rfixest `0'"'
	local ys ""
	local multi = 0
	gettoken first rest : 0, match(par) parse(" ,")
	if "`par'" == "(" {
		local ys `first'
		local 0 `"`rest'"'
		local multi = 1
	}

	syntax varlist(numeric fv ts) [if] [in] [aw pw] , ///
		[ Absorb(string) NOAbsorb VCE(string) CLuster(varlist) Family(string) ///
		  KEEPSINGletons DFRESid PPMLse FULLsample THReads(integer 0) TOL(real 1e-6) ///
		  ESTore(name) ESTArt(integer 1) TMPdir(string) RSCRIPT(string) Level(cilevel) QUIetly TIMEit KEEPfiles ///
		  RESiduals(name) SAVEFE ]

	cap timer off 100
	timer clear 100
	timer on 100

	if !`multi' {
		gettoken ys varlist : varlist
	}
	else unab ys : `ys'
	local xlist `varlist'

	if "`absorb'" != "" & "`noabsorb'" != "" {
		di as error "absorb() and noabsorb may not be combined"
		exit 198
	}
	if "`family'" == "" local family "ols"
	if !inlist("`family'", "ols", "poisson", "negbin") {
		di as error "family() must be ols, poisson or negbin"
		exit 198
	}

	local fevars ""
	local feterms ""
	local fesave ""
	local fenames ""
	local absorb = ustrregexra(`"`absorb'"', "\s*=\s*", "=")
	foreach t of local absorb {
		local nm "."
		if strpos("`t'", "=") {
			local nm = substr("`t'", 1, strpos("`t'", "=") - 1)
			local t = substr("`t'", strpos("`t'", "=") + 1, .)
			confirm new variable `nm'
			local fenames `fenames' `nm'
		}
		local fesave `fesave' `nm'
		local t = subinstr("`t'", "##", "#", .)
		local parts = subinstr("`t'", "#", " ", .)
		local term ""
		foreach p of local parts {
			if substr("`p'", 1, 2) == "i." local p = substr("`p'", 3, .)
			cap unab p : `p'
			if _rc | `: word count `p'' != 1 {
				di as error "absorb(): `t' is not a variable or an interaction of variables (a#b)"
				exit 198
			}
			local fevars `fevars' `p'
			local term = cond("`term'" == "", "`p'", "`term'#`p'")
		}
		local feterms `feterms' `term'
	}
	local fevars : list uniq fevars
	if "`savefe'" != "" {
		cap drop __hdfe*__
		local i = 0
		local fs2 ""
		foreach nm of local fesave {
			local ++i
			if "`nm'" == "." {
				local nm "__hdfe`i'__"
				local fenames `fenames' `nm'
			}
			local fs2 `fs2' `nm'
		}
		local fesave `fs2'
	}
	if "`residuals'" != "" confirm new variable `residuals'
	if "`fenames'`residuals'" != "" & `multi' {
		di as error "residuals() and saved fixed effects (absorb(newvar=var)) need a single dependent variable"
		exit 198
	}
	local dup : list dups fenames
	if "`dup'" != "" | `: list residuals in fenames' & "`residuals'" != "" {
		di as error "the names of saved fixed effects and residuals must be different"
		exit 198
	}

	local vtype "unadjusted"
	if "`weight'" == "pweight" local vtype "robust"
	if "`vce'" != "" {
		gettoken v1 vrest : vce
		if inlist("`v1'", "cluster", "cl", "clu") {
			local vtype "cluster"
			unab cluster : `vrest'
		}
		else if inlist("`v1'", "robust", "r", "hetero") local vtype "robust"
		else if inlist("`v1'", "unadjusted", "ols", "iid") local vtype "unadjusted"
		else {
			di as error "vce() must be unadjusted, robust or cluster varlist"
			exit 198
		}
	}
	else if "`cluster'" != "" local vtype "cluster"

	marksample touse, novarlist

	local xnames ""
	local xvars ""
	if "`xlist'" != "" {
		fvexpand `xlist' if `touse'
		local xnames `r(varlist)'
		fvrevar `xlist' if `touse'
		local xvars `r(varlist)'
	}
	if `: word count `xnames'' != `: word count `xvars'' {
		di as error "could not expand the factor variables"
		exit 198
	}
	if "`xvars'" != "" markout `touse' `xvars'
	if "`fevars'`cluster'" != "" markout `touse' `fevars' `cluster', strok

	local wv ""
	if "`weight'" != "" {
		tempvar wv
		qui gen double `wv' `exp' if `touse'
		qui replace `touse' = 0 if !(`wv' > 0) | missing(`wv')
	}

	_rfx_findR, rscript(`"`rscript'"')
	local rscript `"`r(rscript)'"'

	local tdir `"`tmpdir'"'
	if `"`tdir'"' == "" local tdir `"`c(tmpdir)'"'
	local tdir = subinstr(`"`tdir'"', "\", "/", .)
	if substr(`"`tdir'"', -1, 1) == "/" local tdir = substr(`"`tdir'"', 1, strlen(`"`tdir'"') - 1)
	cap mkdir `"`tdir'"'
	tempfile stamp
	mata: st_local("tok", subinstr(pathbasename(st_local("stamp")), ".", "_"))
	local dir `"`tdir'/rfixest_`tok'"'
	mkdir `"`dir'"'

	cap noisily _rfx_run, ys(`ys') xvars(`xvars') xnames(`xnames') fevars(`fevars') feterms(`feterms') absorb(`absorb') ///
		cl(`cluster') wv(`wv') weight(`weight') vtype(`vtype') family(`family') `keepsingletons' `dfresid' `ppmlse' `fullsample' ///
		estart(`estart') threads(`threads') tol(`tol') touse(`touse') dir(`dir') estore(`estore') level(`level') ///
		`quietly' `timeit' multi(`multi') rscript(`rscript') cmdline(`cmdline') fesave(`fesave') residuals(`residuals')
	local rc = _rc
	if "`keepfiles'" == "" {
		local fs : dir `"`dir'"' files "*"
		foreach f of local fs {
			cap erase `"`dir'/`f'"'
		}
		cap rmdir `"`dir'"'
	}
	else di as text "Files kept in " as result `"`dir'"'
	cap scalar drop _rfx_*
	cap timer off 100
	if `rc' exit `rc'
end

program define _rfx_run, eclass
	syntax , ys(varlist) touse(name) dir(string) vtype(string) family(string) threads(integer) tol(real) multi(integer) ///
		rscript(string) [ xvars(string) xnames(string) fevars(string) feterms(string) absorb(string) cl(string) wv(string) ///
		weight(string) keepsingletons dfresid ppmlse fullsample estart(integer 1) estore(string) level(cilevel) quietly timeit ///
		cmdline(string asis) fesave(string) residuals(string) ]

	tempvar rid
	qui gen long `rid' = _n

	local allv `ys' `xvars' `fevars' `cl' `wv' `rid'
	local allv : list uniq allv

	tempname FR
	frame put `allv' if `touse', into(`FR')
	frame `FR' {
		if !_N {
			di as error "no observations"
			exit 2000
		}
		local i = 0
		local newn ""
		foreach v of local allv {
			local ++i
			local newn `newn' rfxv`i'
		}
		rename (`allv') (`newn')
		cap label values `newn' .
		qui saveold `"`dir'/d.dta"', version(12) replace
	}
	frame drop `FR'

	foreach L in ys xvars fevars cl wv rid {
		local m_`L' ""
		foreach v of local `L' {
			local p : list posof "`v'" in allv
			local m_`L' `m_`L'' rfxv`p'
		}
	}
	local m_fe ""
	foreach t of local feterms {
		local parts = subinstr("`t'", "#", " ", .)
		local mt ""
		foreach v of local parts {
			local p : list posof "`v'" in allv
			local mt = cond("`mt'" == "", "rfxv`p'", "`mt'^rfxv`p'")
		}
		local m_fe `m_fe' `mt'
	}

	tempname fh
	qui file open `fh' using `"`dir'/meta.txt"', write text replace
	file write `fh' "y `m_ys'" _n "x `m_xvars'" _n "fe `m_fe'" _n "cl `m_cl'" _n "w `m_wv'" _n
	file write `fh' "rid `m_rid'" _n "vce `vtype'" _n "family `family'" _n "threads `threads'" _n "tol `tol'" _n
	file write `fh' "keepsingletons " (cond("`keepsingletons'" != "", 1, 0)) _n
	file write `fh' "ppmlse " (cond("`ppmlse'" != "", 1, 0)) _n
	local flags ""
	local savev ""
	foreach nm of local fesave {
		local flags `flags' `=cond("`nm'" != ".", 1, 0)'
		if "`nm'" != "." local savev `savev' `nm'
	}
	file write `fh' "fesave `flags'" _n "resid " (cond("`residuals'" != "", 1, 0)) _n
	file close `fh'

	local rfile `"`dir'/rfixest.R"'
	mata: _rfx_writeR(st_local("rfile"))

	timer off 100
	qui timer list 100
	local t_stata0 = r(t100)
	timer clear 100
	timer on 100
	_rfx_callR, rscript(`rscript') rfile(`rfile') dir(`dir')
	timer off 100
	qui timer list 100
	local t_r = r(t100)
	timer clear 100
	timer on 100

	cap confirm file `"`dir'/ok.txt"'
	if _rc {
		_rfx_rerror, dir(`dir') rscript(`rscript')
		exit 459
	}

	local stored ""
	local nys : word count `ys'
	forvalues k = 1/`nys' {
		local y : word `k' of `ys'
		tempname B V
		mata: _rfx_read(`"`dir'/r`k'.txt"', "`B'", "`V'")

		local cn `xnames'
		if "`feterms'" == "" local cn `cn' _cons
		matrix colnames `B' = `cn'
		matrix colnames `V' = `cn'
		matrix rownames `V' = `cn'
		matrix rownames `B' = `y'

		local i = 0
		local cn2 ""
		local ncol = 0
		foreach c of local cn {
			local ++i
			if `V'[`i',`i'] == 0 & !strpos("`c'", "b.") {
				local ++ncol
				if !strpos("`c'", "o.") local c "o.`c'"
			}
			local cn2 `cn2' `c'
		}
		matrix colnames `B' = `cn2'
		matrix colnames `V' = `cn2'
		matrix rownames `V' = `cn2'

		tempvar tk
		qui gen byte `tk' = `touse' & !missing(`y')
		cap confirm file `"`dir'/x`k'.txt"'
		if !_rc & "`fullsample'" == "" mata: _rfx_excl(`"`dir'/x`k'.txt"', "`tk'")
		qui count if `touse' & !missing(`y')
		local n_cand = r(N)

		local N = scalar(_rfx_N)
		local K = scalar(_rfx_K)
		local ming = .
		local g = 0
		foreach c of local cl {
			local ++g
			local ming = min(`ming', scalar(_rfx_G`g'))
		}
		local dfr = `N' - `K'
		if "`vtype'" == "cluster" & "`dfresid'" == "" local dfr = `ming' - 1

		ereturn post `B' `V', obs(`N') dof(`dfr') depname(`y') esample(`tk')
		ereturn local cmd "rfixest"
		ereturn local cmdline `"`cmdline'"'
		ereturn local depvar "`y'"
		ereturn local absvars "`absorb'"
		ereturn local family "`family'"
		ereturn local vce "`vtype'"
		ereturn local vcetype = cond("`vtype'" == "cluster", "Clustered", cond("`vtype'" == "robust", "Robust", ""))
		ereturn local clustvar "`cl'"
		ereturn local title = cond("`family'" == "ols", "Linear regression with fixed effects (fixest)", ///
			cond("`family'" == "poisson", "Poisson regression with fixed effects (fixest)", "Negative binomial regression with fixed effects (fixest)"))
		ereturn local properties "b V"
		ereturn local predict "rfixest_p"
		ereturn local fe_saved "`savev'"
		ereturn scalar n_fe = `: word count `feterms''
		ereturn scalar fe_all = `: word count `savev'' == `: word count `feterms''
		foreach s in df_a ssr r2 r2_a r2_within r2_a_within rmse ll r2_p {
			ereturn scalar `s' = scalar(_rfx_`s')
		}
		ereturn scalar n_collin = `ncol'
		ereturn scalar N_full = `n_cand'
		ereturn scalar num_singletons = `n_cand' - `N'
		ereturn scalar df_r = `dfr'
		ereturn scalar K = `K'
		ereturn scalar rank = `K' - scalar(_rfx_df_a)
		if "`vtype'" == "cluster" {
			ereturn scalar N_clust = `ming'
			local g = 0
			foreach c of local cl {
				local ++g
				ereturn scalar N_clust`g' = scalar(_rfx_G`g')
			}
		}
		ereturn scalar time_r = `t_r'

		if "`quietly'" == "" _rfx_show, level(`level')

		if "`estore'" != "" | `multi' {
			if `multi' local nm = cond("`estore'" != "", "`estore'`=`k'+`estart'-1'", "rfx`=`k'+`estart'-1'")
			else local nm "`estore'"
			estimates store `nm'
			local stored `stored' `nm'
		}
	}

	if "`residuals'`savev'" != "" {
		cap confirm file `"`dir'/save1.txt"'
		if !_rc {
			tempname sh
			file open `sh' using `"`dir'/save1.txt"', read text
			file read `sh' nr
			file read `sh' nc
			file close `sh'
			local newv `residuals' `savev'
			foreach v of local newv {
				qui gen double `v' = .
			}
			mata: _rfx_readsave(`"`dir'/save1.bin"', `nr', `nc', "`newv'")
			if "`residuals'" != "" label variable `residuals' "Residuals of rfixest `ys'"
			local i = 0
			foreach nm of local fesave {
				local ++i
				if "`nm'" == "." continue
				local term : word `i' of `feterms'
				label variable `nm' "Fixed effect of `term'"
			}
		}
	}

	timer off 100
	qui timer list 100
	local t_stata1 = r(t100)
	timer clear 100

	if "`stored'" != "" di as text "Estimates stored: " as result "`stored'"

	if "`timeit'" != "" {
		tempname t
		file open `t' using `"`dir'/time.txt"', read text
		file read `t' line
		while r(eof) == 0 {
			tokenize `line'
			local t_`1' = `2'
			file read `t' line
		}
		file close `t'
		di as text "Time (s): Stata before R " as result %6.2f `t_stata0' ///
			as text "  | R " as result %6.2f `t_r' ///
			as text " [start " as result %5.2f `t_lib' as text ", read " as result %5.2f `t_read' ///
			as text ", fit " as result %6.2f `t_fit' as text ", write " as result %5.2f `t_write' as text "]" ///
			as text "  | Stata after R " as result %6.2f `t_stata1'
	}
end

program define _rfx_callR
	syntax , rscript(string) rfile(string) dir(string) [ args(string) ]
	if "`c(os)'" == "Windows" {
		local wr = subinstr(`"`rscript'"', "/", "\", .)
		local wf = subinstr(`"`rfile'"', "/", "\", .)
		local wd = subinstr(`"`dir'"', "/", "\", .)
		tempname fh
		qui file open `fh' using `"`dir'/run.vbs"', write text replace
		file write `fh' "On Error Resume Next" _n
		file write `fh' `"Set sh = CreateObject("WScript.Shell")"' _n
		file write `fh' `"rc = sh.Run("""`wr'"" --vanilla ""`wf'"" ""`wd'"" `args'", 0, True)"' _n
		file write `fh' "If Err.Number <> 0 Then rc = -1" _n
		file write `fh' `"Set fs = CreateObject("Scripting.FileSystemObject")"' _n
		file write `fh' `"Set f = fs.CreateTextFile("`wd'\done.tmp", True)"' _n
		file write `fh' "f.WriteLine rc" _n
		file write `fh' "f.Close" _n
		file write `fh' `"fs.MoveFile "`wd'\done.tmp", "`wd'\done.txt""' _n
		file close `fh'
		winexec wscript.exe "`wd'\run.vbs"
		local n = 0
		while 1 {
			cap confirm file `"`dir'/done.txt"'
			if !_rc continue, break
			local ++n
			sleep `=cond(`n' < 40, 25, 100)'
		}
	}
	else {
		qui shell "`rscript'" --vanilla "`rfile'" "`dir'" `args' > "`dir'/stdout.txt" 2>&1
	}
end

program define _rfx_rerror
	syntax , dir(string) rscript(string)
	local msg ""
	cap confirm file `"`dir'/err.txt"'
	if !_rc {
		tempname fh
		file open `fh' using `"`dir'/err.txt"', read text
		file read `fh' msg
		file close `fh'
	}
	if strpos(`"`msg'"', "RFIXEST_NO_FIXEST") {
		di as error "The R package fixest is not installed. Type {stata rfixest setup} once to install it."
		exit
	}
	cap confirm file `"`dir'/log.txt"'
	if _rc {
		di as error `"rfixest could not run R (`rscript')."'
		di as error `"Install R from {browse "https://cran.r-project.org":https://cran.r-project.org}, then type {stata rfixest setup},"'
		di as error `"or give the path to Rscript, for example: {cmd:rfixest setup, rscript("C:/Program Files/R/R-4.4.1/bin/x64/Rscript.exe")}"'
		exit
	}
	local code ""
	cap confirm file `"`dir'/done.txt"'
	if !_rc {
		tempname fh
		file open `fh' using `"`dir'/done.txt"', read text
		file read `fh' code
		file close `fh'
		local code = strtrim(`"`code'"')
	}
	di as error "R stopped with an error" cond("`code'" != "", " (exit code `code')", "") ". R output:"
	type `"`dir'/log.txt"'
end

program define _rfx_findR, rclass
	syntax [, rscript(string)]
	local r `"`rscript'"'
	if `"`r'"' == "" local r `"${RFIXEST_RSCRIPT}"'
	if `"`r'"' == "" {
		cap confirm file `"`c(sysdir_personal)'rfixest_rscript.txt"'
		if !_rc {
			tempname fh
			file open `fh' using `"`c(sysdir_personal)'rfixest_rscript.txt"', read text
			file read `fh' r
			file close `fh'
			local r = strtrim(`"`r'"')
			cap confirm file `"`r'"'
			if _rc local r ""
		}
	}
	if `"`r'"' == "" {
		local best ""
		if "`c(os)'" == "Windows" {
			local la : env LOCALAPPDATA
			local pf : env ProgramFiles
			local bases `""`pf'/R" "C:/Program Files/R" "`la'/Programs/R""'
			foreach b of local bases {
				local b = subinstr(`"`b'"', "\", "/", .)
				local ds ""
				cap local ds : dir `"`b'"' dirs "R-*"
				local ds : list sort ds
				foreach d of local ds {
					foreach s in "bin/x64/Rscript.exe" "bin/Rscript.exe" {
						cap confirm file `"`b'/`d'/`s'"'
						if !_rc {
							local best `"`b'/`d'/`s'"'
							continue, break
						}
					}
				}
			}
		}
		else {
			foreach p in /opt/homebrew/bin/Rscript /usr/local/bin/Rscript /usr/bin/Rscript /Library/Frameworks/R.framework/Resources/bin/Rscript {
				cap confirm file "`p'"
				if !_rc & `"`best'"' == "" local best "`p'"
			}
		}
		local r `"`best'"'
	}
	if `"`r'"' == "" local r "Rscript"
	local r = subinstr(`"`r'"', "\", "/", .)
	if `"`rscript'"' == "" global RFIXEST_RSCRIPT `"`r'"'
	return local rscript `"`r'"'
end

program define _rfx_setup
	syntax [, RSCRIPT(string) TMPdir(string)]
	if `"`rscript'"' != "" global RFIXEST_RSCRIPT ""
	_rfx_findR, rscript(`"`rscript'"')
	local r `"`r(rscript)'"'
	local tdir `"`tmpdir'"'
	if `"`tdir'"' == "" local tdir `"`c(tmpdir)'"'
	local tdir = subinstr(`"`tdir'"', "\", "/", .)
	if substr(`"`tdir'"', -1, 1) == "/" local tdir = substr(`"`tdir'"', 1, strlen(`"`tdir'"') - 1)
	tempfile stamp
	mata: st_local("tok", subinstr(pathbasename(st_local("stamp")), ".", "_"))
	local dir `"`tdir'/rfixest_setup_`tok'"'
	mkdir `"`dir'"'
	local rfile `"`dir'/rfixest.R"'
	mata: _rfx_writeR(st_local("rfile"))
	di as text "Using R at " as result `"`r'"'
	di as text "Checking the R package fixest (installing it from CRAN if needed; this can take a few minutes the first time)..."
	_rfx_callR, rscript(`r') rfile(`rfile') dir(`dir') args(setup)
	local ok = 0
	cap confirm file `"`dir'/ok.txt"'
	if !_rc {
		tempname fh
		file open `fh' using `"`dir'/ok.txt"', read text
		file read `fh' ver
		file close `fh'
		local ok = 1
	}
	if `ok' {
		global RFIXEST_RSCRIPT `"`r'"'
		cap mkdir `"`c(sysdir_personal)'"'
		tempname fh
		cap file open `fh' using `"`c(sysdir_personal)'rfixest_rscript.txt"', write text replace
		if !_rc {
			file write `fh' `"`r'"' _n
			file close `fh'
		}
		di as result "rfixest is ready." as text " R: " as result `"`r'"' as text "; fixest " as result "`ver'" as text "."
	}
	else _rfx_rerror, dir(`dir') rscript(`r')
	local fs : dir `"`dir'"' files "*"
	foreach f of local fs {
		cap erase `"`dir'/`f'"'
	}
	cap rmdir `"`dir'"'
	if !`ok' exit 459
end

program define _rfx_show
	syntax [, Level(cilevel)]
	tempname b V
	matrix `b' = e(b)
	matrix `V' = e(V)
	local names : colnames `b'
	local k = colsof(`b')
	local df = e(df_r)
	local z = invttail(`df', (100 - `level') / 200)

	di
	di as text "`e(title)'" _col(51) "Number of obs" _col(68) "=" as result %11.0fc e(N)
	di as text "Dependent variable: " as result abbrev("`e(depvar)'", 28) _col(51) as text "Residual df" _col(68) "=" as result %11.0fc e(df_r)
	local abs = cond("`e(absvars)'" == "", "(none)", "`e(absvars)'")
	local se = "`e(vce)'" + cond("`e(clustvar)'" != "", " (`e(clustvar)')", "")
	if "`e(family)'" == "ols" {
		di as text "Fixed effects: " as result abbrev("`abs'", 34) _col(51) as text "R-squared" _col(68) "=" as result %11.4f e(r2)
		di as text "Std. errors: " as result abbrev("`se'", 36) _col(51) as text "Adj R-squared" _col(68) "=" as result %11.4f e(r2_a)
		if "`e(absvars)'" != "" di _col(51) as text "Within R-sq." _col(68) "=" as result %11.4f e(r2_within)
		di _col(51) as text "Root MSE" _col(68) "=" as result %11.4f e(rmse)
	}
	else {
		di as text "Fixed effects: " as result abbrev("`abs'", 34) _col(51) as text "Pseudo R-sq." _col(68) "=" as result %11.4f e(r2_p)
		di as text "Std. errors: " as result abbrev("`se'", 36) _col(51) as text "Log likelihood" _col(68) "=" as result %11.1f e(ll)
	}
	if "`e(vce)'" == "cluster" di _col(51) as text "Clusters" _col(68) "=" as result %11.0fc e(N_clust)
	di as text "{hline 13}{c TT}{hline 66}"
	di as text %12s abbrev("`e(depvar)'", 12) " {c |}" _col(19) "Coef." _col(29) "Std. err." _col(42) "t" _col(48) "P>|t|" _col(57) "[`level'% conf. interval]"
	di as text "{hline 13}{c +}{hline 66}"
	forvalues i = 1/`k' {
		local n : word `i' of `names'
		local bi = `b'[1, `i']
		local si = sqrt(`V'[`i', `i'])
		if `si' == 0 | missing(`si') {
			di as text %12s abbrev("`n'", 12) " {c |}" _col(19) as text cond(strpos("`n'", "b."), "(base)", "(omitted)")
			continue
		}
		local t = `bi' / `si'
		local p = 2 * ttail(`df', abs(`t'))
		local st = cond(`p' < .01, "***", cond(`p' < .05, "**", cond(`p' < .1, "*", "")))
		di as text %12s abbrev("`n'", 12) " {c |}" as result ///
			_col(16) %9.0g `bi' _col(28) %9.0g `si' _col(38) %7.2f `t' _col(47) %6.3f `p' ///
			_col(57) %9.0g `bi' - `z' * `si' _col(68) %9.0g `bi' + `z' * `si' " " as text "`st'"
	}
	di as text "{hline 13}{c BT}{hline 66}"
	di as text "* p<0.10, ** p<0.05, *** p<0.01"
	if e(n_collin) > 0 & e(n_collin) < . di as text "Omitted because of collinearity: " as result e(n_collin)
	if e(num_singletons) > 0 & e(num_singletons) < . di as text "Singletons or perfectly predicted observations dropped: " as result strtrim(string(e(num_singletons), "%12.0fc"))
	if "`e(absvars)'" != "" di as text "Degrees of freedom absorbed by fixed effects: " as result strtrim(string(e(df_a), "%12.0fc"))
end

mata:
void _rfx_read(string scalar fn, string scalar bname, string scalar Vname)
{
	real scalar fh, nb
	string matrix l
	string rowvector t
	real colvector b, v
	b = J(0, 1, .)
	v = J(0, 1, .)
	fh = fopen(fn, "r")
	while ((l = fget(fh)) != J(0, 0, "")) {
		t = tokens(l)
		if (t[1] == "b") b = b \ strtoreal(t[2])
		else if (t[1] == "v") v = v \ strtoreal(t[2])
		else st_numscalar("_rfx_" + t[2], strtoreal(t[3]))
	}
	fclose(fh)
	nb = rows(b)
	st_matrix(bname, b')
	st_matrix(Vname, colshape(v, nb)')
}

void _rfx_readsave(string scalar fn, real scalar n, real scalar nc, string scalar vars)
{
	real scalar fh, j
	real matrix x
	string rowvector v
	transmorphic C
	fh = fopen(fn, "r")
	C = bufio()
	bufbyteorder(C, 2)
	x = bufget(C, fread(fh, n * nc * 8), 0, "%8z", nc, n)'
	fclose(fh)
	v = tokens(vars)
	for (j = 1; j <= cols(v); j++) st_store(x[., 1], v[j], x[., j + 1])
}

void _rfx_excl(string scalar fn, string scalar var)
{
	real colvector x
	x = strtoreal(cat(fn))
	st_store(x, var, J(rows(x), 1, 0))
}
end

* begin embedded R code
mata:
string colvector _rfx_rcode()
{
	string colvector s
	s = J(184, 1, "")
	s[1] = "# rfixest 1.2.0: R side of the Stata command rfixest. Called by rfixest.ado; not meant to be run by hand."
	s[2] = "t_start <- proc.time()[[" + char(34) + "elapsed" + char(34) + "]]"
	s[3] = "argv <- commandArgs(trailingOnly = TRUE)"
	s[4] = "dir <- argv[1]"
	s[5] = "mode <- if (length(argv) > 1) argv[2] else " + char(34) + "fit" + char(34)
	s[6] = ""
	s[7] = "logcon <- file(file.path(dir, " + char(34) + "log.txt" + char(34) + "), open = " + char(34) + "wt" + char(34) + ")"
	s[8] = "sink(logcon)"
	s[9] = "sink(logcon, type = " + char(34) + "message" + char(34) + ")"
	s[10] = "options(warn = 1)"
	s[11] = ""
	s[12] = "fail <- function(msg) {"
	s[13] = "  writeLines(msg, file.path(dir, " + char(34) + "err.txt" + char(34) + "))"
	s[14] = "  cat(msg, " + char(34) + "\n" + char(34) + ")"
	s[15] = "  sink(type = " + char(34) + "message" + char(34) + ")"
	s[16] = "  sink()"
	s[17] = "  close(logcon)"
	s[18] = "  quit(save = " + char(34) + "no" + char(34) + ", status = 1)"
	s[19] = "}"
	s[20] = ""
	s[21] = "user_lib <- Sys.getenv(" + char(34) + "R_LIBS_USER" + char(34) + ")"
	s[22] = "if (nzchar(user_lib) && dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))"
	s[23] = ""
	s[24] = "if (mode == " + char(34) + "setup" + char(34) + ") {"
	s[25] = "  if (!requireNamespace(" + char(34) + "fixest" + char(34) + ", quietly = TRUE)) {"
	s[26] = "    lib <- .libPaths()[1]"
	s[27] = "    if (file.access(lib, 2) != 0) {"
	s[28] = "      lib <- user_lib"
	s[29] = "      dir.create(lib, recursive = TRUE, showWarnings = FALSE)"
	s[30] = "      .libPaths(c(lib, .libPaths()))"
	s[31] = "    }"
	s[32] = "    res <- try(install.packages(" + char(34) + "fixest" + char(34) + ", lib = lib, repos = " + char(34) + "https://cloud.r-project.org" + char(34) + "), silent = TRUE)"
	s[33] = "    if (!requireNamespace(" + char(34) + "fixest" + char(34) + ", quietly = TRUE)) fail(" + char(34) + "Could not install the R package fixest. See the R output above." + char(34) + ")"
	s[34] = "  }"
	s[35] = "  writeLines(as.character(utils::packageVersion(" + char(34) + "fixest" + char(34) + ")), file.path(dir, " + char(34) + "ok.txt" + char(34) + "))"
	s[36] = "  sink(type = " + char(34) + "message" + char(34) + ")"
	s[37] = "  sink()"
	s[38] = "  close(logcon)"
	s[39] = "  quit(save = " + char(34) + "no" + char(34) + ", status = 0)"
	s[40] = "}"
	s[41] = ""
	s[42] = "if (!requireNamespace(" + char(34) + "fixest" + char(34) + ", quietly = TRUE)) fail(" + char(34) + "RFIXEST_NO_FIXEST" + char(34) + ")"
	s[43] = ""
	s[44] = "res <- tryCatch({"
	s[45] = "  meta <- strsplit(readLines(file.path(dir, " + char(34) + "meta.txt" + char(34) + ")), " + char(34) + " " + char(34) + ", fixed = TRUE)"
	s[46] = "  names(meta) <- vapply(meta, " + char(96) + "[" + char(96) + ", " + char(34) + char(34) + ", 1L)"
	s[47] = "  meta <- lapply(meta, " + char(96) + "[" + char(96) + ", -1L)"
	s[48] = "  one <- function(k) meta[[k]][1L]"
	s[49] = ""
	s[50] = "  suppressPackageStartupMessages(library(fixest))"
	s[51] = "  setFixest_notes(FALSE)"
	s[52] = "  t_lib <- proc.time()[[" + char(34) + "elapsed" + char(34) + "]]"
	s[53] = ""
	s[54] = "  d <- foreign::read.dta(file.path(dir, " + char(34) + "d.dta" + char(34) + "), convert.dates = FALSE, convert.factors = FALSE, missing.type = FALSE)"
	s[55] = "  t_read <- proc.time()[[" + char(34) + "elapsed" + char(34) + "]]"
	s[56] = ""
	s[57] = "  Y   <- meta" + char(36) + "y"
	s[58] = "  X   <- meta" + char(36) + "x"
	s[59] = "  FE  <- meta" + char(36) + "fe"
	s[60] = "  CL  <- meta" + char(36) + "cl"
	s[61] = "  W   <- meta" + char(36) + "w"
	s[62] = "  vce <- one(" + char(34) + "vce" + char(34) + ")"
	s[63] = "  fam <- one(" + char(34) + "family" + char(34) + ")"
	s[64] = "  nth <- as.integer(one(" + char(34) + "threads" + char(34) + "))"
	s[65] = "  tol <- as.numeric(one(" + char(34) + "tol" + char(34) + "))"
	s[66] = "  keep_single <- one(" + char(34) + "keepsingletons" + char(34) + ") == " + char(34) + "1" + char(34)
	s[67] = ""
	s[68] = "  rhs <- if (length(X)) paste(X, collapse = " + char(34) + "+" + char(34) + ") else " + char(34) + "1" + char(34)
	s[69] = "  make_fml <- function(ys) {"
	s[70] = "    lhs <- if (length(ys) > 1L) paste0(" + char(34) + "c(" + char(34) + ", paste(ys, collapse = " + char(34) + "," + char(34) + "), " + char(34) + ")" + char(34) + ") else ys"
	s[71] = "    f <- paste(lhs, " + char(34) + "~" + char(34) + ", rhs)"
	s[72] = "    if (length(FE)) f <- paste(f, " + char(34) + "|" + char(34) + ", paste(FE, collapse = " + char(34) + "+" + char(34) + "))"
	s[73] = "    as.formula(f, env = globalenv())"
	s[74] = "  }"
	s[75] = "  same_na <- length(Y) == 1L || all(vapply(Y[-1L], function(v) identical(is.na(d[[v]]), is.na(d[[Y[1L]]])), TRUE))"
	s[76] = "  fml <- make_fml(Y)"
	s[77] = ""
	s[78] = "  vc <- switch(vce,"
	s[79] = "    cluster = as.formula(paste(" + char(34) + "~" + char(34) + ", paste(CL, collapse = " + char(34) + "+" + char(34) + ")), env = globalenv()),"
	s[80] = "    robust  = " + char(34) + "hetero" + char(34) + ","
	s[81] = "    " + char(34) + "iid" + char(34) + ")"
	s[82] = ""
	s[83] = "  args <- list(fml = fml, data = d, vcov = vc, fixef.tol = tol)"
	s[84] = "  if (nth > 0) args" + char(36) + "nthreads <- nth"
	s[85] = "  if (length(W)) args" + char(36) + "weights <- as.formula(paste(" + char(34) + "~" + char(34) + ", W), env = globalenv())"
	s[86] = "  if (one(" + char(34) + "ppmlse" + char(34) + ") == " + char(34) + "1" + char(34) + ") args" + char(36) + "ssc <- ssc(adj = FALSE)"
	s[87] = ""
	s[88] = "  fit <- function(a) switch(fam,"
	s[89] = "    poisson = {"
	s[90] = "      a" + char(36) + "fixef.rm <- if (keep_single) " + char(34) + "perfect" + char(34) + " else " + char(34) + "both" + char(34)
	s[91] = "      do.call(fepois, a)"
	s[92] = "    },"
	s[93] = "    negbin = {"
	s[94] = "      a" + char(36) + "fixef.rm <- if (keep_single) " + char(34) + "perfect" + char(34) + " else " + char(34) + "both" + char(34)
	s[95] = "      do.call(fenegbin, a)"
	s[96] = "    },"
	s[97] = "    {"
	s[98] = "      a" + char(36) + "fixef.rm <- if (keep_single) " + char(34) + "none" + char(34) + " else " + char(34) + "singleton" + char(34)
	s[99] = "      do.call(feols, a)"
	s[100] = "    })"
	s[101] = "  if (same_na) {"
	s[102] = "    M <- fit(args)"
	s[103] = "    models <- if (inherits(M, " + char(34) + "fixest_multi" + char(34) + ")) as.list(M) else list(M)"
	s[104] = "  } else {"
	s[105] = "    models <- lapply(Y, function(v) {"
	s[106] = "      a <- args"
	s[107] = "      a" + char(36) + "fml <- make_fml(v)"
	s[108] = "      fit(a)"
	s[109] = "    })"
	s[110] = "  }"
	s[111] = "  t_fit <- proc.time()[[" + char(34) + "elapsed" + char(34) + "]]"
	s[112] = "  n_d <- nrow(d)"
	s[113] = "  num <- function(x) format(x, digits = 17L, trim = TRUE, scientific = TRUE)"
	s[114] = "  w_norm <- function(used) if (length(W)) length(used) / sum(d[[W]][used]) else 1"
	s[115] = ""
	s[116] = "  for (k in seq_along(models)) {"
	s[117] = "    m <- models[[k]]"
	s[118] = "    V <- vcov(m)"
	s[119] = "    b <- coef(m)"
	s[120] = "    N <- nobs(m)"
	s[121] = "    K <- m" + char(36) + "nparams"
	s[122] = "    nb <- length(b)"
	s[123] = ""
	s[124] = "    pos <- c(X, if (!length(FE)) " + char(34) + "(Intercept)" + char(34) + ")"
	s[125] = "    bf <- setNames(numeric(length(pos)), pos)"
	s[126] = "    Vf <- matrix(0, length(pos), length(pos), dimnames = list(pos, pos))"
	s[127] = "    est <- intersect(pos, names(b))"
	s[128] = "    bf[est] <- b[est]"
	s[129] = "    Vf[est, est] <- V[est, est]"
	s[130] = ""
	s[131] = "    used <- obs(m)"
	s[132] = "    G <- vapply(CL, function(v) length(unique(d[[v]][used])), 0)"
	s[133] = ""
	s[134] = "    if (fam == " + char(34) + "ols" + char(34) + ") {"
	s[135] = "      r <- r2(m, c(" + char(34) + "r2" + char(34) + ", " + char(34) + "ar2" + char(34) + ", " + char(34) + "wr2" + char(34) + ", " + char(34) + "war2" + char(34) + "))"
	s[136] = "      stat <- c(N = N, K = K, df_a = K - nb, ssr = m" + char(36) + "ssr,"
	s[137] = "                r2 = r[[" + char(34) + "r2" + char(34) + "]], r2_a = r[[" + char(34) + "ar2" + char(34) + "]], r2_within = r[[" + char(34) + "wr2" + char(34) + "]], r2_a_within = r[[" + char(34) + "war2" + char(34) + "]],"
	s[138] = "                rmse = sqrt(m" + char(36) + "ssr * w_norm(used) / (N - K)), ll = NA_real_, r2_p = NA_real_)"
	s[139] = "    } else {"
	s[140] = "      r <- r2(m, " + char(34) + "pr2" + char(34) + ")"
	s[141] = "      stat <- c(N = N, K = K, df_a = K - nb, ssr = NA_real_,"
	s[142] = "                r2 = NA_real_, r2_a = NA_real_, r2_within = NA_real_, r2_a_within = NA_real_,"
	s[143] = "                rmse = NA_real_, ll = as.numeric(logLik(m)), r2_p = r[[" + char(34) + "pr2" + char(34) + "]])"
	s[144] = "    }"
	s[145] = "    if (length(G)) stat <- c(stat, setNames(G, paste0(" + char(34) + "G" + char(34) + ", seq_along(G))))"
	s[146] = "    stat <- c(stat, n_collin = length(pos) - length(est))"
	s[147] = ""
	s[148] = "    writeLines(c("
	s[149] = "      paste(" + char(34) + "s" + char(34) + ", names(stat), num(unname(stat))),"
	s[150] = "      paste(" + char(34) + "b" + char(34) + ", num(unname(bf))),"
	s[151] = "      paste(" + char(34) + "v" + char(34) + ", num(as.vector(Vf)))"
	s[152] = "    ), file.path(dir, paste0(" + char(34) + "r" + char(34) + ", k, " + char(34) + ".txt" + char(34) + ")))"
	s[153] = ""
	s[154] = "    fesave <- meta" + char(36) + "fesave == " + char(34) + "1" + char(34)
	s[155] = "    doresid <- one(" + char(34) + "resid" + char(34) + ") == " + char(34) + "1" + char(34)
	s[156] = "    if (k == 1L && (any(fesave) || doresid)) {"
	s[157] = "      cols <- list(as.double(d[[one(" + char(34) + "rid" + char(34) + ")]][used]))"
	s[158] = "      if (doresid) cols <- c(cols, list(as.double(resid(m))))"
	s[159] = "      if (any(fesave)) {"
	s[160] = "        fv <- fixef(m, notes = FALSE, sorted = FALSE)"
	s[161] = "        for (j in which(fesave)) cols <- c(cols, list(as.double(fv[[j]][m" + char(36) + "fixef_id[[j]]])))"
	s[162] = "      }"
	s[163] = "      con <- file(file.path(dir, " + char(34) + "save1.bin" + char(34) + "), " + char(34) + "wb" + char(34) + ")"
	s[164] = "      writeBin(unlist(cols), con, size = 8L, endian = " + char(34) + "little" + char(34) + ")"
	s[165] = "      close(con)"
	s[166] = "      writeLines(as.character(c(length(used), length(cols))), file.path(dir, " + char(34) + "save1.txt" + char(34) + "))"
	s[167] = "    }"
	s[168] = ""
	s[169] = "    ex <- d[[one(" + char(34) + "rid" + char(34) + ")]][setdiff(seq_len(n_d), used)]"
	s[170] = "    if (length(ex)) writeLines(as.character(ex), file.path(dir, paste0(" + char(34) + "x" + char(34) + ", k, " + char(34) + ".txt" + char(34) + ")))"
	s[171] = "  }"
	s[172] = ""
	s[173] = "  t_end <- proc.time()[[" + char(34) + "elapsed" + char(34) + "]]"
	s[174] = "  writeLines(paste(c(" + char(34) + "lib" + char(34) + ", " + char(34) + "read" + char(34) + ", " + char(34) + "fit" + char(34) + ", " + char(34) + "write" + char(34) + "),"
	s[175] = "                   c(t_lib - t_start, t_read - t_lib, t_fit - t_read, t_end - t_fit)),"
	s[176] = "             file.path(dir, " + char(34) + "time.txt" + char(34) + "))"
	s[177] = "  " + char(34) + "ok" + char(34)
	s[178] = "}, error = function(e) conditionMessage(e))"
	s[179] = ""
	s[180] = "if (!identical(res, " + char(34) + "ok" + char(34) + ")) fail(paste(" + char(34) + "Error in R:" + char(34) + ", res))"
	s[181] = "writeLines(" + char(34) + "ok" + char(34) + ", file.path(dir, " + char(34) + "ok.txt" + char(34) + "))"
	s[182] = "sink(type = " + char(34) + "message" + char(34) + ")"
	s[183] = "sink()"
	s[184] = "close(logcon)"
	return(s)
}

void _rfx_writeR(string scalar fn)
{
	real scalar fh, i
	string colvector s
	s = _rfx_rcode()
	if (fileexists(fn)) unlink(fn)
	fh = fopen(fn, "w")
	for (i = 1; i <= rows(s); i++) fput(fh, s[i])
	fclose(fh)
}
end
* end embedded R code

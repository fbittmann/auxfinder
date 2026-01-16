*! version 1.0  16jan2026  Felix Bittmann
program define auxfinder, rclass
	syntax varlist(fv min=1) [if] [in], ///
	TESTvars(varlist fv min=1) ///
	[CORRValue(real 0.30)] ///
	[CORRMiss(real 0.15)] ///
	[MAXmiss(real 0.10)] ///
	[gain(real 0.15)] ///
	[CATlimit(integer 7)] ///	Testing whether a variable is categorical or not
	[SAVing(string)] ///
	[folds(integer 20)] ///
	[seed(integer -1)] ///
	[SKIPlasso] ///
	[ELASTICnet] ///
	[lassoopts(string)] ///
	[NONLinear(real -1)]
	
	if "`seed'" != "-1" {
		set seed `seed'
	}
	
	if `nonlinear' >= 0 {
		if `nonlinear' >= 1 {
			local temp = `nonlinear' * 100
			di as result "caution, your specification of nonlinear() is potentially incorrect. A value of `nonlinear' indicates that a nonlinear term is only considered if R² improves by `temp'%, compared to the linear model specification."
		}
		local nonlinear = `nonlinear' + 1	//Convert to multiplierfactor
	}
	
	local maxmiss = 1 + `maxmiss'		//Convert to multiplierfactor
	local targets `varlist'
	tempfile origdata
	qui save `origdata'
	
	marksample touse, novarlist		//IF+IN selection
	qui keep if `touse'

	*** Remove vars without cases ***
	local temp
	foreach VAR of local testvars {
		sum `VAR', meanonly
		if r(N) > 0 {
			local temp `temp' `VAR'
		}
		else {
			di as result "Removed `VAR' as it only contains missing values"
		}
	}
	local testvars `temp'
	local ntotaltestvars : list sizeof testvars
	
	*** Remove targetvars without missing cases ***
	local temp
	foreach VAR of local targets {
		if substr("`VAR'", 1, 2) == "i." {
			local t = substr("`VAR'", 3, .)
			qui count if `t' == .
			local nmiss = r(N)
			if `nmiss' > 0 {
			qui levelsof `t'
				if r(r) == 2 {
					di as result "caution: `VAR' has exactly 2 levels yet is declared as categorial. This is not recommended."
				}
			}
		}
		else {
			qui count if `VAR' == .
			local nmiss = r(N)
		}
		if `nmiss' > 0 {
			local temp `temp' `VAR'
		}
		else {
			di as result "`VAR' removed from targets because it has no missing values"
		}
	}
	local targets `temp'
	local n_targetvars : list sizeof targets	
	
	*** Sorting Testvars into continuous and categorical ***
	local testvars_cont
	local testvars_cat
	foreach VAR of local testvars {
		if substr("`VAR'", 1, 2) == "c." {
			local res = substr("`VAR'", 3, .)
			local testvars_cont `testvars_cont' `res'
		}
		else {
			if substr("`VAR'", 1, 2) == "i." {
				local res = substr("`VAR'", 3, .)
				local testvars_cat `testvars_cat' `res'
			}
			else {
				qui levelsof `VAR'
				if r(r) > `catlimit' {
					local testvars_cont `testvars_cont' `VAR'
				}
				else {
					local testvars_cat `testvars_cat' `VAR'
				}
			}
		}
	}
	
	*** Sorting Targetvars into continuous and categorial ***
	local targetvars_cont
	local targetvars_cat
	foreach VAR of local targets {
		if substr("`VAR'", 1, 2) == "i." {
			local res = substr("`VAR'", 3, .)
			local targetvars_cat `targetvars_cat' `res'
		}
		else {
			if substr("`VAR'", 1, 2) == "c." {
				local res = substr("`VAR'", 3, .)
				local targetvars_cont `targetvars_cont' `res'
			}
			else {
				local targetvars_cont `targetvars_cont' `VAR'
			}
		}
	}	
	
	*** Cleaning ***
	local targetvars_cont : list uniq targetvars_cont
	local targetvars_cat : list uniq targetvars_cat
	local testvars_cont : list uniq testvars_cont
	local testvars_cat : list uniq testvars_cat
	local testvars_cont : list testvars_cont - targetvars_cont
	local testvars_cont : list testvars_cont - targetvars_cat
	local testvars_cat : list testvars_cat - targetvars_cont
	local testvars_cat : list testvars_cat - targetvars_cat
	
	local t1 : list sizeof testvars_cont
	local t2 : list sizeof testvars_cat
	local n_testvars = `t1' + `t2'
	
	*** Count valid cases ***
	tempfile missdata
	tempname name
	postfile `name' str16 testvar nvalid using `missdata'
	local temp `testvars_cont' `testvars_cat'
	foreach VAR of local temp {
		qui count if !missing(`VAR')
		post `name' ("`VAR'") (r(N))
	}
	postclose `name'	
	di "Part 1/5 done"	
	
	************************************************************************
	*** Prediction of missingness ***
	************************************************************************
	local tvars `targetvars_cont' `targetvars_cat'
	foreach VAR of local tvars {
		qui recode `VAR' (.=1) (nonmiss=0), gen(missing_x_`VAR')
	}
	local miss_catvars
	local miss_contvars
	local miss_sqterms
	tempfile corr_missingness
	tempname mname	
	postfile `mname' str16 targetvar str32 testvar corr_miss mode using `corr_missingness'
	foreach VAR of local tvars {
		local lasso1vars_`VAR'
		qui count if !missing(`VAR')
		local b = r(N)
		foreach t of local testvars_cont {
			qui count if !missing(`t')
			local a = r(N) * `maxmiss'
			if `a' >= `b' {
				qui reg missing_x_`VAR' `t'
				local c1 = sqrt(e(r2))
				if `c1' >= `corrmiss' & `c1' != . {
					post `mname' ("`VAR'") ("`t'") (`c1') (0)
					local miss_contvars `miss_contvars' `t'
					local lasso1vars_`VAR' `lasso1vars_`VAR'' `t'
					if `nonlinear' >= 1 {
						cap drop `t'²
						qui gen `t'² = `t' * `t'
						qui regress missing_x_`VAR' `t' `t'²
						local c2 = sqrt(e(r2))
						if `c2' > `c1' * `nonlinear' & `c2' != . {
							post `mname' ("`VAR'") ("`t'²") (`c2') (0)
							local miss_sqterms `miss_sqterms' `t'²
							local lasso1vars_`VAR' `lasso1vars_`VAR'' `t'²
						}
					
					}
				}
			}
		}
		foreach t of local testvars_cat {
			qui count if !missing(`t')
			local a = r(N) * `maxmiss'
			if `a' >= `b' {
				qui reg missing_x_`VAR' i.`t'
				local c1 = sqrt(e(r2))
				if `c1' >= `corrmiss' & `c1' != . {
					post `mname' ("`VAR'") ("`t'") (`c1') (0)
					local miss_catvars `local miss_catvars' `t'
					local lasso1vars_`VAR' `lasso1vars_`VAR'' i.`t'
				}
			}
		}
	}
	postclose `mname'
	di "Part 2/5 done"
	
	* Lasso *
	local miss_sqterms : list uniq miss_sqterms		//remove duplicates
	local miss_contvars : list uniq miss_contvars
	local miss_catvars : list uniq miss_catvars
	*local lasso1vars_`VAR' : list uniq lasso1vars_`VAR'

	if "`skiplasso'" == "" {
		if "`elasticnet'" == "" {
			local lmodel lasso
		}
		else {
			local lmodel elasticnet
		}
	macro drop miss_catvarsx
	macro drop test0
	if "`miss_catvars'" != "" {
		qui vl create test0 = (`miss_catvars')
		qui vl substitute miss_catvarsx = i.test0
	}

	tempfile lasso_missing
	tempname name
	postfile `name' str16 targetvar str32 testvar selected_missingness using `lasso_missing'
	foreach VAR of varlist missing_x_* {
		local tempname = substr("`VAR'",11,.)
		if "`lasso1vars_`tempname''" == "" {
			post `name' ("`VAR'") ("") (-9)
			continue
		}
		*di "`VAR' `lasso1vars_`tempname''"
		cap `lmodel' linear `VAR' `lasso1vars_`tempname'', folds(`folds') `lassoopts' nolog
		if _rc == 0 {
			local temp = e(post_sel_vars)
			foreach element of local temp {
				if "`element'" != "`VAR'" {
					post `name' ("`tempname'") ("`element'") (1)
				}
			}
			if "`noisily'" != "" {
				di as result "Dependent variable: `VAR'"
				`lmodel'
				lassocoef
				di ""
			}
		}
		else {
			post `name' ("`tempname'") ("") (-9)
		}
	}
	postclose `name'
	}
	di "Part 3/5 done"
	************************************************************************
	*** Testing for correlations and missingness ***
	************************************************************************
	tempfile file
	tempname name
	postfile `name' str16 targetvar str16 testvar corr_value infogain mode using `file'
	local keepvars_cont
	local keepvars_cat
	local sqterms
	
	*Init empty values*
	foreach VAR of local tvars {
		local lasso2vars_`VAR'
	}
	
	* CAT (Target) - CAT *
	foreach TAVAR of local targetvars_cat {
		foreach TEVAR of local testvars_cat {
			qui count if `TAVAR' == . & !missing(`TEVAR')
			local nextra = r(N)
			qui count if `TAVAR' == .
			local nmissing = r(N)
			qui count if !missing(`TEVAR')
			local a = r(N) * `maxmiss'
			qui count if !missing(`TAVAR')
			local b = r(N)			
			local propextra = `nextra' / `nmissing'
			
			if `propextra' > `gain' & `a' >= `b' {
				cap mlogit `TAVAR' i.`TEVAR'			//MLOGIT
				if _rc == 0 {
					local c1 = sqrt(abs(e(r2_p)))
					if "`noisily'" != ""{
						mlogit
					}
				}
				else {
					local c1 = 0
				}	
				if `c1' >= `corrvalue' & `c1' != . {
					local keepvars_cat `keepvars_cat' `TEVAR'
					local lasso2vars_`TAVAR' `lasso2vars_`TAVAR'' i.`TEVAR'
					post `name' ("`TAVAR'") ("`TEVAR'") (`c1') (`propextra') (1)
				}
			}
		}
	}
	* CONT (Target) - CONT *
	foreach TAVAR of local targetvars_cont {
		foreach TEVAR of local testvars_cont {
			qui count if `TAVAR' == . & !missing(`TEVAR')
			local nextra = r(N)
			qui count if `TAVAR' == .
			local nmissing = r(N)
			qui count if !missing(`TEVAR')
			local a = r(N) * `maxmiss'
			qui count if !missing(`TAVAR')
			local b = r(N)			
			local propextra = `nextra' / `nmissing'
			
			if `propextra' > `gain' & `a' >= `b' {
				cap regress `TAVAR' `TEVAR'			//REGRESS
				if _rc == 0 {
					local c1 = sqrt(e(r2))
					if "`noisily'" != ""{
						regress
					}
				}
				else {
					local c1 = 0
				}	
				if `c1' >= `corrvalue' & `c1' != . {
					local keepvars_cont `keepvars_cont' `TEVAR'
					local lasso2vars_`TAVAR' `lasso2vars_`TAVAR'' `TEVAR'
					post `name' ("`TAVAR'") ("`TEVAR'") (`c1') (`propextra') (1)
					if `nonlinear' >= 1 {
						cap drop `TEVAR'²
						qui gen `TEVAR'² = `TEVAR' * `TEVAR'
						qui regress `TAVAR' `TEVAR' `TEVAR'²
						local c2 = sqrt(e(r2))
						if `c2' > `c1' * `nonlinear' & `c2' != . {
							local lasso2vars_`TAVAR' `lasso2vars_`TAVAR'' `TEVAR'²
							post `name' ("`TAVAR'") ("`TEVAR'²") (`c2') (`propextra') (1)
							local sqterms `sqterms' `TEVAR'²
						}
					
					}
				}

			}
		}
	}	
	* CAT (Target) - CONT *
	foreach TAVAR of local targetvars_cat {
		foreach TEVAR of local testvars_cont {
			qui count if `TAVAR' == . & !missing(`TEVAR')
			local nextra = r(N)
			qui count if `TAVAR' == .
			local nmissing = r(N)
			qui count if !missing(`TEVAR')
			local a = r(N) * `maxmiss'
			qui count if !missing(`TAVAR')
			local b = r(N)			
			local propextra = `nextra' / `nmissing'
			
			if `propextra' > `gain' & `a' >= `b' {
				cap mlogit `TAVAR' `TEVAR'		//MLOGIT
				if _rc == 0 {
					local c1 = sqrt(abs(e(r2_p)))
					if "`noisily'" != ""{
						mlogit
					}
				}
				else {
					local c1 = 0
				}	
				if `c1' >= `corrvalue' & `c1' != . {
					local keepvars_cont `keepvars_cont' `TEVAR'
					local lasso2vars_`TAVAR' `lasso2vars_`TAVAR'' `TEVAR'
					post `name' ("`TAVAR'") ("`TEVAR'") (`c1') (`propextra') (1)
					if `nonlinear' >= 1 {
						cap drop `TEVAR'²
						qui gen `TEVAR'² = `TEVAR' * `TEVAR'
						qui mlogit `TAVAR' `TEVAR' `TEVAR'²
						local c2 = sqrt(e(r2_p))
						if `c2' > `c1' * `nonlinear' & `c2' != . {
							local lasso2vars_`TAVAR' `lasso2vars_`TAVAR'' `TEVAR'²
							post `name' ("`TAVAR'") ("`TEVAR'²") (`c2') (`propextra') (1)
							local sqterms `sqterms' `TEVAR'²
						}
					
					}
				}

			}
		}
	}	
	* CONT (Target) - CAT *
	foreach TAVAR of local targetvars_cont {
		foreach TEVAR of local testvars_cat {
			qui count if `TAVAR' == . & !missing(`TEVAR')
			local nextra = r(N)
			qui count if `TAVAR' == .
			local nmissing = r(N)
			qui count if !missing(`TEVAR')
			local a = r(N) * `maxmiss'
			qui count if !missing(`TAVAR')
			local b = r(N)			
			local propextra = `nextra' / `nmissing'
			
			if `propextra' > `gain' & `a' >= `b' {
				cap regress `TAVAR' i.`TEVAR'		//REGRESS
				if _rc == 0 {
					local c1 = sqrt(e(r2))
					if "`noisily'" != ""{
						regress
					}
				}
				else {
					local c1 = 0
				}	
				if `c1' >= `corrvalue' & `c1' != . {
					local keepvars_cat `keepvars_cat' `TEVAR'
					local lasso2vars_`TAVAR' `lasso2vars_`TAVAR'' i.`TEVAR'
					post `name' ("`TAVAR'") ("`TEVAR'") (`c1') (`propextra') (1)
				}
			}
		}
	}
	di "Part 4/5 done"
	postclose `name'
	local sqterms : list uniq sqterms		//remove duplicates
	local keepvars_cont : list uniq keepvars_cont
	local keepvars_cat : list uniq keepvars_cat
	
	local t1 : list sizeof keepvars_cont
	local t2 : list sizeof keepvars_cat
	local n_corrvars = `t1' + `t2'	
	************************************************************************
	*** LASSO ***
	************************************************************************
	
	if "`skiplasso'" == "" {
		if "`elasticnet'" == "" {
			local lmodel lasso
		}
		else {
			local lmodel elasticnet
		}
	macro drop catvars1
	macro drop test1
	if "`keepvars_cat'" != "" {
		qui vl create test1 = (`keepvars_cat')
		qui vl substitute catvars1 = i.test1
	}

	tempfile lasso
	tempname name
	postfile `name' str16 targetvar str32 testvar selected_value mode using `lasso'
	foreach VAR of local targetvars_cont {
		*local testing : list keepvars_cont - VAR		//Indepvar cannot be depvar
		if "`lasso2vars_`VAR''" == "" {
			post `name' ("`VAR'") ("") (-9) (1)
			continue
		}
		*di "`lasso2vars_`VAR''"
		cap `lmodel' linear `VAR' `lasso2vars_`VAR'', folds(`folds') `lassoopts' nolog
		if inlist(_rc, 0) {
			local temp = e(post_sel_vars)
			foreach element of local temp {
				post `name' ("`VAR'") ("`element'") (1) (1)
			}
			if "`noisily'" != "" {
				di as result "Dependent variable: `VAR'"
				`lmodel'
				lassocoef
				di ""
			}
		}
		else {
			post `name' ("`VAR'") ("") (-9) (1)
		}
	}
	postclose `name'
	}
	di "Part 5/5 done"
	************************************************************************
	*** Present results ***
	************************************************************************
	use `file', clear
	qui merge 1:1 targetvar testvar using `corr_missingness', nogen
	qui merge m:1 testvar using `missdata', nogen keep(1 3)
	qui count
	if r(N) == 0 {
		di as result "No relevant auxiliaries found for any of the target variables"
		di as result "Try lowering the values of corrmiss() and corrvalue() for more results"
		exit
	}
	if "`skiplasso'" == "" {		//LASSO specified
		qui merge 1:1 targetvar testvar using `lasso_missing', nogen
		qui merge 1:1 targetvar testvar using `lasso', nogen
	}
	else {
		gen selected_value = 0
		gen selected_missing = 0
	}
		qui bysort targetvar: egen fail_value = min(selected_value)
		qui bysort targetvar: egen fail_missing = min(selected_missing)
		qui replace selected_missing = fail_missing if missing(selected_missing) & fail_missing < 0
		qui replace selected_value = fail_value if missing(selected_value) & fail_value < 0 
		sort targetvar testvar
		qui keep if !missing(corr_value) | !missing(corr_miss)
		format %6.3f corr_*
		format %6.2f infogain
		qui tostring corr_miss corr_value, gen(a b) force usedisplay
		rename (corr_miss corr_value) (corr_miss_numerical corr_value_numerical)
		qui gen catvar = 0
		foreach VAR of local testvars_cat {
			qui replace catvar = 1 if testvar == "`VAR'"
		}
		label define catvar 0 "No" 1 "Yes"
		label value catvar catvar
		if "`skiplasso'" == "" {
			qui gen symbol_miss = ""
			qui replace symbol_miss = "*" if selected_missing == 1
			qui replace symbol_miss = "†" if selected_missing == -9
			qui egen a_1 = concat(a symbol_miss) if a != "."
			qui replace a_1 = "." if a == "."
			
			qui gen symbol_value = ""
			qui replace symbol_value = "*" if selected_value == 1
			qui replace symbol_value = "†" if selected_value == -9
			qui egen b_1 = concat(b symbol_value) if b != "."
			qui replace b_1 = "." if b == "."
			order targetvar testvar catvar a_1 b_1 infogain
			rename a_1 corr_miss
			rename b_1 corr_value
			list targetvar testvar nvalid catvar corr_miss corr_value infogain ///
				, noobs sepby(targetvar) abbrev(16)
			di as result "    * Selected by `lmodel' 	† `lmodel' failed	"
		}
		else {	
			order targetvar testvar catvar corr_miss corr_value infogain
			list targetvar testvar nvalid catvar corr* infogain ///
				, noobs sepby(targetvar) abbrev(16)
			di as result "    Lasso selection not computed"
		}


	*** Macros ***
	gen squared = substr(testvar,-2,.) == "²"
	
	local strong_cont
	cap drop take
	gen take = !missing(corr_miss_numerical) & !missing(corr_value_numerical) ///
		& selected_value == 1 & selected_missing == 1 & cat == 0 & squared == 0
	qui count if take == 1
	local n = r(N)
	if `n' > 0 {
		preserve
		qui keep if take == 1
		forvalues i = 1/`n' {
			local res = testvar[`i']
			local strong_cont `strong_cont' `res'
		}
	restore	
	}
	
	
	local strong_cat
	cap drop take
	gen take = !missing(corr_miss_numerical) & !missing(corr_value_numerical) ///
		& selected_value == 1 & selected_missing == 1 & cat == 1 & squared == 0
	qui count if take == 1
	local n = r(N)
	if `n' > 0 {
		preserve
		qui keep if take == 1
		forvalues i = 1/`n' {
			local res = testvar[`i']
			local strong_cat `strong_cat' `res'
		}
	restore
	}
	
	
	local strong_squared
	cap drop take
	gen take = !missing(corr_miss_numerical) & !missing(corr_value_numerical) ///
		& selected_value == 1 & selected_missing == 1 & cat == 0 & squared == 1
	qui count if take == 1
	local n = r(N)
	if `n' > 0 {
		preserve
		qui keep if take == 1
		forvalues i = 1/`n' {
			local res = testvar[`i']
			local res = subinstr("`res'", "²", "",1)
			local strong_squared `strong_squared' substr(`res'
		}
	restore
	}
	
	local strong_cont : list uniq strong_cont
	local strong_cat : list uniq strong_cat
	local strong_squared : list uniq strong_squared	
	local n_strong_cont : word count `strong_cont'
	local n_strong_cat : word count `strong_cat'
	local n_strong_squared : word count `strong_squared'
	if "`skiplasso'" == "" {
		di as text "Strong continuous auxiliary variables (`n_strong_cont'):"
		di as result "`strong_cont'"
		di as text "Strong categorical auxiliary variables (`n_strong_cat'):"
		di as result "`strong_cat'"
		if `nonlinear' >= 1 {
			di as text "Strong squared auxiliary variables (`n_strong_squared'):"
			di as result "`strong_squared'"
		}
	}
	
	*** Weak ***
	local weak_cont
	cap drop take
	gen take = ((!missing(corr_miss_numerical) & selected_missing == 1) | ///
		(!missing(corr_value_numerical) & selected_value == 1)) & ///
		cat == 0 & squared == 0
	qui count if take == 1
	local n = r(N)
	if `n' > 0 {
		preserve
		qui keep if take == 1
		forvalues i = 1/`n' {
			local res = testvar[`i']
			local weak_cont `weak_cont' `res'
		}
	restore
	}
	
	
	local weak_cat
	cap drop take
	gen take = ((!missing(corr_miss_numerical) & selected_missing == 1) | ///
		(!missing(corr_value_numerical) & selected_value == 1)) & ///
		cat == 1 & squared == 0
	qui count if take == 1
	local n = r(N)
	if `n' > 0 {
		preserve
		qui keep if take == 1
		forvalues i = 1/`n' {
			local res = testvar[`i']
			local weak_cat `weak_cat' `res'
		}
	restore
	}
	
	
	local weak_squared
	cap drop take
	gen take = ((!missing(corr_miss_numerical) & selected_missing == 1) | ///
		(!missing(corr_value_numerical) & selected_value == 1)) & ///
		cat == 0 & squared == 1
	qui count if take == 1
	local n = r(N)
	if `n' > 0 {
		preserve
		qui keep if take == 1
		forvalues i = 1/`n' {
			local res = testvar[`i']
			local res = subinstr("`res'", "²", "",1)
			local weak_squared `weak_squared' `res'
		}
	restore
	}
	
	local weak_cont : list uniq weak_cont
	local weak_cat : list uniq weak_cat
	local weak_squared : list uniq weak_squared	
	local n_weak_cat : word count `weak_cat'
	local n_weak_cont : word count `weak_cont'
	local n_weak_squared : word count `weak_squared'
	if "`skiplasso'" == "" {
		di as text "Weak continuous auxiliary variables (`n_weak_cont'):"
		di as result "`weak_cont'"
		di as text "Weak categorical auxiliary variables (`n_weak_cat'):"
		di as result "`weak_cat'"
		if `nonlinear' >= 1 {
			di as text "Weak squared auxiliary variables (`n_weak_squared'):"
			di as result "`weak_squared'"
		}
	}
	
	* Returns *
	return scalar n_weak_cat = `n_weak_cat'
	return scalar n_weak_cont = `n_weak_cont'
	return scalar n_weak_squared = `n_weak_squared'
	return scalar n_strong_cat = `n_strong_cat'
	return scalar n_strong_cont = `n_strong_cont'
	return scalar n_strong_squared = `n_strong_squared'
	return local weak_cat `weak_cat'
	return local weak_cont `weak_cont'
	return local weak_squared `weak_squared'
	return local strong_cat `strong_cat'
	return local strong_cont `strong_cont'
	return local strong_squared `strong_squared'
	return local testvars `testvars'
	return scalar n_total_tested = `ntotaltestvars'

	* Restore original data *
	if "`saving'" != "" {
		save `saving'
	}
	qui use `origdata', clear	
	end

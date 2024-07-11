global path "D:\Projects\#2034_ICON_chemo_radio"
global res "N:\ICON_all\PROJECTS\#2034_ICON_chemo_radio\Results\paper2\resubmit"
cd "$path\Data\analysis2"
use stage3colon, clear
codebook age
drop if age >90
cap drop age_grp
recode age (18/50=1) (50/60=2) (60/70=3) (70/80=4) (80/90=5), gen(age_grp)
label define age_grp 1 "18.0-49.9" 2 "50.0-59.9" 3 "60.0-69.9" 4 "70.0-79.9" 5 "80.0-90.0", modify
label values age_grp age_grp
egen float missing = rowmiss(age income2015_5 edu2015_5 employ2015_5 ses2015_5 ///
							 sex cci1 nonwhite t_stage n_stage final_route)
tab missing, m
misstable pattern age ses2015_5 sex cci1 nonwhite t_stage n_stage final_route, freq
replace missing = 1 if missing>0

tab final_route, m
gen ep = 1 if final_route == 1
replace ep = 0 if final_route !=1 & final_route !=.
gen tstage = 0 if t_stage == 1 | t_stage == 2
replace tstage = 1 if t_stage == 3 | t_stage == 4
label define tstage 0 "1-2" 1 "3-4"
label values tstage tstage
replace n_stage = n_stage - 1
label define n_stage 0 "1" 1 "2"
label values n_stage n_stage

table1_mc, total(before) by(missing) ///
vars(ses2015_5 cat %4.1f\ income2015_5 cat %4.1f\ employ2015_5 cat %4.1f\ 		///
edu2015_5 cat %4.1f\ yod cat %4.1f\ age conts %4.1f\ age_grp cat %4.1f\ 		///
sex cat %4.1f\ nonwhite cat %4.1f\ ep bin %4.1f\ cci1 bin %4.1f\ 		///
tstage cat %4.1f\ n_stage cat %4.1f\ optimal cat %4.1f\ surgery cat %4.1f\) 	///
nospace onecol miss saving("$res/table1.xlsx", sheet("missing", replace))

gen death1y = 1 if (deathdate - diagdate) !=. & death == 1 & (deathdate - diagdate) <365
replace death1y = 0 if death1y == .
gen death6m_adjuvant = 1 if (deathdate - surgdate) !=. & death == 1 & surgery == 1 & (deathdate - surgdate) <180
replace death6m_adjuvant = 0 if death6m_adjuvant == .

keep if missing == 0
foreach e in income employ edu ses {
	table1_mc, total(before) by(`e'2015_5) ///
	vars(yod cat %4.1f\ age conts %4.1f\ age_grp cat %4.1f\ sex cat %4.1f\ ///
	nonwhite cat %4.1f\ ep bin %4.1f\ cci1 bin %4.1f\ tstage cat %4.1f\ ///
	n_stage cat %4.1f\ optimal cat %4.1f\ surgery cat %4.1f\) ///
	nospace miss onecol saving("$res/table1.xlsx", sheet("`e'_complete", replace))
}


foreach e in income employ edu ses {
	table1_mc if surgery == 1, total(before) by(`e'2015_5) ///
	vars(postop_chemo cat %4.1f\ ) ///
	nospace miss onecol saving("$res/table1.xlsx", sheet("`e'_adj", replace))
	table1_mc if surgery == 0, total(before) by(`e'2015_5) ///
	vars(death1y bin %4.1f\ ) ///
	nospace miss onecol saving("$res/table1.xlsx", sheet("`e'_death6m", replace))
	table1_mc if surgery == 1 & postop_chemo == 0, total(before) by(`e'2015_5) ///
	vars(death6m_adjuvant bin %4.1f\ ) ///
	nospace miss onecol saving("$res/table1.xlsx", sheet("`e'_death6m_adj", replace))
}

gen dif = surgdate -  diagdate if surgery == 1
gen dif1 = postop_chemodate - surgdate if postop_chemo == 1
distplot dif, xlabel(0(30)360) xtitle(Days between diagnosis and surgery) name(surg, replace)
distplot dif1, xlabel(0(30)180) xtitle(Days between surgery and adjuvant chemotherapy) name(chemo, replace)
graph combine surg chemo, rows(1)
graph export "C:\Users\lshsl7\OneDrive - London School of Hygiene and Tropical Medicine\BK_CRC_paper\results\resubmit\displot.svg", as(svg)  replace
	
/*multinominal logistic for optimal, logistic surgery*/
preserve
clear
tempfile res
save `res', replace emptyok
restore

mkspline zage = age, cubic nknots(4)
matrix knots = r(knots)
local knots

forvalues j = 1 / 4 {
	local k = knots[1, `j']
	local knots `knots' `k'
}

sort age
gen ages = _n + 24 if _n <= 66
mkspline wage = ages, cubic knots(`knots')
foreach e in income edu employ ses {
	mlogit optimal 									///
	i.`e'2015_5#c.(zage1 zage2 zage3) 				///
	zage1 zage2 zage3 i.sex i.`e'2015_5 i.nonwhite  ///
	i.cci1 i.tstage i.n_stage i.ep ///
	, vce(cluster creg_code)
	
	sum `e'2015_5, d
	marginscontplot age (zage1 zage2 zage3) `e'2015_5, var1(ages(wage1 wage2 wage3)) ///
	at2(`r(min)' `r(max)') ci margopts(predict(outcome(0))) nograph saving("$res/opt0_prob_`e'1", replace)

	sum `e'2015_5, d
	marginscontplot age (zage1 zage2 zage3) `e'2015_5, var1(ages(wage1 wage2 wage3)) ///
	at2(`r(min)' `r(max)') ci margopts(predict(outcome(1))) nograph saving("$res/opt1_prob_`e'1", replace)
	
	sum `e'2015_5, d
	marginscontplot age (zage1 zage2 zage3) `e'2015_5, var1(ages(wage1 wage2 wage3)) ///
	at2(`r(min)' `r(max)') ci margopts(predict(outcome(2))) nograph saving("$res/opt2_prob_`e'1", replace)
	

	mlogit optimal 									///
	i.`e'2015_5#c.(zage1 zage2 zage3) 				///
	zage1 zage2 zage3 i.sex i.`e'2015_5 i.nonwhite ///
	, vce(cluster creg_code)
	
	sum `e'2015_5, d
	marginscontplot age (zage1 zage2 zage3) `e'2015_5, var1(ages(wage1 wage2 wage3)) ///
	at2(`r(min)' `r(max)') ci margopts(predict(outcome(0))) nograph saving("$res/opt0_prob_`e'2", replace)

	sum `e'2015_5, d
	marginscontplot age (zage1 zage2 zage3) `e'2015_5, var1(ages(wage1 wage2 wage3)) ///
	at2(`r(min)' `r(max)') ci margopts(predict(outcome(1))) nograph saving("$res/opt1_prob_`e'2", replace)
	
	sum `e'2015_5, d
	marginscontplot age (zage1 zage2 zage3) `e'2015_5, var1(ages(wage1 wage2 wage3)) ///
	at2(`r(min)' `r(max)') ci margopts(predict(outcome(2))) nograph saving("$res/opt2_prob_`e'2", replace)
	
	
	logit surgery									///
	i.`e'2015_5#c.(zage1 zage2 zage3) 				///
	zage1 zage2 zage3 i.sex i.`e'2015_5 i.nonwhite ///
	, vce(cluster creg_code)
		
	sum `e'2015_5, d
	marginscontplot age (zage1 zage2 zage3) `e'2015_5, var1(ages(wage1 wage2 wage3)) ///
	at2(`r(min)' `r(max)') ci nograph saving("$res/surg_prob_`e'", replace)

	sum `e'2015_5, d
	local k = `r(max)'
	preserve
	partpred or if `e'2015_5 == `k', ///
	for(`k'.`e'2015_5 `k'.`e'2015_5#c.zage1 `k'.`e'2015_5#c.zage2 `k'.`e'2015_5#c.zage3) ///
	ci(or_lci or_uci) eform
	keep if or!=.
	keep age or*
	gen strata = "`e'`k'"
	gen outcome = "Surgery"
	append using `res'
	save `res', replace
	restore	
	
	foreach mf in tstage cci1 n_stage ep {
		mediate (surgery i.`e'2015_5#c.(zage1 zage2 zage3) zage1 zage2 zage3 i.sex i.nonwhite, logit) ///
		(`mf' zage1 zage2 zage3 i.sex i.nonwhite, logit) (`e'2015_5), control(1)
		return list
		matrix m = r(table)
		matrix list m
		matselrc m m, r(1 5 6) c(4 8 12)
		putexcel set "$res/mediate_effect.xlsx", sheet("`e'_surg_`mf'", replace) modify
		putexcel A1=matrix(m)

		estat proportion, force
		return list
		matrix m = r(table)
		matrix list m
		matselrc m m, r(1 5 6) c(4)
		putexcel set "$res/mediate_prop.xlsx", sheet("`e'_surg_`mf'", replace) modify
		putexcel A1=matrix(m)
	}
}

**#////explore mediation of treatment and survival
preserve
drop if diagdate > date("06/02/2016","DMY") | (enddate - diagdate <365.24*3 & death !=1)
gen fup3y = diagdate + 365.24*3
egen enddate3 = rowmin(fup3y enddate deathdate)
gen death3y =1 if death == 1
replace death3y = 0 if death3y == .
tab death3y, m
tab death3y ses2015_5, m col
cap drop zage* wage* ages 
mkspline zage = age, cubic nknots(4)
matrix knots = r(knots)
local knots

forvalues j = 1 / 4 {
	local k = knots[1, `j']
	local knots `knots' `k'
}

foreach e in income edu employ ses {
	foreach mf in tstage cci1 n_stage ep surgery postop_chemo {
		mediate (death3y zage1 zage2 zage3 i.sex i.nonwhite, logit) ///
		(`mf' i.`e'2015_5#c.(zage1 zage2 zage3) zage1 zage2 zage3 i.sex i.nonwhite, logit) (`e'2015_5), control(1)	
		return list
		matrix m = r(table)
		matrix list m
		matselrc m m, r(1 5 6) c(4 8 12)
		putexcel set "$res/mediate_effect.xlsx", sheet("`e'_death3y_`mf'", replace) modify
		putexcel A1=matrix(m)

		estat proportion, force
		return list
		matrix m = r(table)
		matrix list m
		matselrc m m, r(1 5 6) c(4)
		putexcel set "$res/mediate_prop.xlsx", sheet("`e'_death3y_`mf'", replace) modify
		putexcel A1=matrix(m)
	}
}
restore

/*logistic adjuvent chemo*/
keep if surgery == 1
cap drop zage* wage* ages 
mkspline zage = age, cubic nknots(4)
matrix knots = r(knots)
local knots

forvalues j = 1 / 4 {
	local k = knots[1, `j']
	local knots `knots' `k'
}

sort age
gen ages = _n + 24 if _n <= 66
mkspline wage = ages, cubic knots(`knots')
foreach e in income edu employ ses {
	logit postop_chemo								///
	i.`e'2015_5#c.(zage1 zage2 zage3) 				///
	zage1 zage2 zage3 i.sex i.`e'2015_5 i.nonwhite ///
	, vce(cluster creg_code)
	
	sum `e'2015_5, d
	marginscontplot age (zage1 zage2 zage3) `e'2015_5, var1(ages(wage1 wage2 wage3)) ///
	at2(`r(min)' `r(max)') ci nograph saving("$res/adjuvant_prob_`e'", replace)

	sum `e'2015_5, d
	local k = `r(max)'
	preserve
	partpred or if `e'2015_5 == `k', ///
	for(`k'.`e'2015_5 `k'.`e'2015_5#c.zage1 `k'.`e'2015_5#c.zage2 `k'.`e'2015_5#c.zage3) ///
	ci(or_lci or_uci) eform
	keep if or!=.
	keep age or*
	gen strata = "`e'`k'"
	gen outcome = "Adjuvant Chemo"
	append using `res'
	save `res', replace
	restore		
	
	foreach mf in tstage cci1 n_stage ep {
		mediate (postop_chemo i.`e'2015_5#c.(zage1 zage2 zage3) zage1 zage2 zage3 i.sex i.nonwhite, logit) ///
		(`mf' zage1 zage2 zage3 i.sex i.nonwhite, logit) (`e'2015_5), control(1)	
		return list
		matrix m = r(table)
		matrix list m
		matselrc m m, r(1 5 6) c(4 8 12)
		putexcel set "$res/mediate_effect.xlsx", sheet("`e'_chemo_`mf'", replace) modify
		putexcel A1=matrix(m)

		estat proportion, force
		return list
		matrix m = r(table)
		matrix list m
		matselrc m m, r(1 5 6) c(4)
		putexcel set "$res/mediate_prop.xlsx", sheet("`e'_chemo_`mf'", replace) modify
		putexcel A1=matrix(m)
	}
}

use `res', clear
save "$res\or.dta", replace
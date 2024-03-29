global data "D:\Projects\#2034_ICON_chemo_radio"
global res "N:\ICON_all\PROJECTS\#2034_ICON_chemo_radio\Results\paper1\BJC_R1"
cd "$data\Data\analysis1"

cap log close
log using "$data\Results\logfiles\crc_treatment_analysis1_multistate_2015_16_R1.log", text replace
display "$S_TIME  $S_DATE"
set seed 1521462
/*S Ling Nov 2021, data analysis for CRC treatment access, multistate*/
///*complete-case analysis*/
///*update on 3rd May 2022: adjust for the same set of covariates for every stage*/
preserve
clear
tempfile hr
save `hr', emptyok replace
restore

preserve
clear
tempfile df
save `df', emptyok replace
restore

preserve
clear
tempfile dif
save `dif', emptyok replace
restore

foreach c in colon rectal {
	use `c', clear
	gen ydiag = year(diagdate)
	tab ydiag, m
	keep if ydiag>=2015
	tab ydiag tnm_stage, m
	egen float missing = rowmiss(tnm_stage age female nonwhite final_route income2015)
	tab missing, m
	keep if missing == 0
	gen dead = 1 if death == 1 & (deathdate <= diagdate + 365.24 & deathdate!=.) ///
	| (enddate <= diagdate + 365.24 & enddate!=.)
	replace dead = 0 if dead ==.
	gen censoring = diagdate + 365.24
	egen censor_date = rowmin(censoring enddate deathdate)
	egen treat_date = rowmin(surgdate chemodate radiodate) if treat == 1
	replace treat_date = censor_date if treat == 0 
	gen r1 = runiform(0.1,0.9)
	gen r2 = runiform(0.1,0.5)
	gen treat_time = treat_date - diagdate
	replace treat_time = treat_time + r1 if ( treat_time == 0 & treat == 1)
	gen dead_date = censor_date
	gen dead_time = dead_date - diagdate
	replace dead_time = dead_time + r1 if (dead_time == 0 & dead == 1)
	tab tnm_stage if dead_date<= treat_date & treat == 1 & dead == 1
	replace dead_time = treat_time + r2 if (dead_time == treat_time & dead == 1 & treat == 1)
	tab tnm_stage if dead_date<= treat_date & treat == 1 & dead == 1
	replace treat_time = dead_time if (dead_date == treat_date & dead == 1 & treat == 0)
	format *date %td

	
	tempfile `c'
	save ``c'', replace
	
	////*updated on 3rd May 2022: interaction between stage and income2015_5 is difficult to implement using stmerlin*/
	////*so I stick to stratified analyses*/	
	////*updated on 18th May 2022: we decide to stratified by early stages (I or II) and advanced stages (III or IV)*/
	
forvalues i = 1(1)5 {
	use ``c'', clear
	keep if tnm_stage == `i'
	msset, id(pseudo_patientid) states(treat dead) times(treat_time dead_time)
	matrix tmat = r(transmatrix)
	mat list tmat
		stset _stop, enter(_start) failure(_status==1)
		forvalues k = 1(1)3 {
			forvalues j = 1(1)5 {
				cap noisily stpm2 if _trans`k' == 1, df(`j')  scale(h)
				if _rc == 0 {
					preserve
					estat ic
					return list
					matrix m = r(S)
					clear
					svmat long m, names(matcol) 
					gen cancer = "`c'"
					gen stage = `i'
					gen trans = `k'
					gen df = `j'
					append using `df'
					save `df', replace
					restore
				}
			}
		}
	}

}


use `df', clear
sort cancer stage trans df
bysort cancer stage trans: gen n_test = _N
bysort cancer stage trans: egen max_df = max(df)
drop if max_df != n_test & max_df == df
bysort cancer stage trans: egen minAIC = min(mAIC)
bysort cancer stage trans: egen minBIC = min(mBIC)
save "$res/df_1516.dta", replace
keep if minAIC == mAIC | minBIC == mBIC
bysort cancer stage trans: egen mindf = min(df)
keep if df == mindf
save `df', replace
drop mll0 mll mdf n_test max_df minAIC minBIC mindf
rename m* *
reshape wide N AIC BIC df, i(cancer stage) j(trans)
replace cancer = proper(cancer)
label define stage 1 "I" 2 "II" 3 "III" 4 "IV" 5 "Missing"
label values stage stage
export excel using "$res\tables.xlsx", sheet("TableS5_df_1516") sheetmodify firstrow(variables)

foreach c in colon rectal {
	forvalues i = 1(1)5 {
		use ``c'', clear
		keep if tnm_stage == `i'
		foreach var of varlist income2015_5 final_route  {
			tab `var', gen(`var')
		}
		rcsgen age, gen(ages) orthog knots(5 35 65 95)
		global Kage `r(knots)'
		matrix Mage = r(R)

		msset, id(pseudo_patientid) states(treat dead) times(treat_time dead_time)
		matrix tmat = r(transmatrix)
		mat list tmat
		
		/*
		msboxes, transmat(tmat) id(pseudo_patientid) ///
		xvalues(0.2 0.7 0.45) yvalues(0.7 0.7 0.2) ///
		statenames("Diagnosis" "Treatment" "Dead") ///
		boxwidth(0.3)
		graph save Graph "`c'`i'.gph", replace
		*/
		stset _stop, enter(_start) failure(_status==1)
		forvalues k = 1(1)3 {
			
			preserve
			use `df', clear
			sum df if cancer == "`c'" & stage == `i' & trans == `k'
			local df1 = `r(mean)'
			restore
			
			stmerlin income2015_52 income2015_53 income2015_54 income2015_55 ///
			ages1 ages2 ages3 female hf_p mi_p dmc_p pulmonary_p nonwhite  ///
			final_route1 final_route3 final_route4 final_route5 final_route6 ///
			if _trans`k' == 1, distribution(rp)	df(`df1')
		
			estimates store m`k'
								
			preserve
			parmest, eform fast
			gen stage = `i'
			gen cancer = "`c'"
			gen trans = `k'
			append using `hr'
			save `hr', replace
			restore
			}
		
		forvalues k = 75(5)75 { 
		display "`c' cancer" " stage `i'"
		
		cap drop tt
		range tt 0 365 366
		preserve
		rcsgen, scalar(`k') knots(${Kage}) rmatrix(Mage) gen(v)
			
		predictms, transmatrix(tmat) models(m1 m2 m3) ///
		probability timevar(tt) los diff ci ///
		at1(ages1 `=v1' ages2 `=v2' ages3 `=v3')  ///
		at2(income2015_55 1 ages1 `=v1' ages2 `=v2' ages3 `=v3') 
		
		keep _*prob* _*los* tt
		drop if tt ==.
		gen stage = `i'
		gen cancer = "`c'"
		gen age = `k'
		append using `dif'
		save `dif', replace
		restore
/*		
		cap drop tt
		range tt 0 360 13
		foreach f in female nonwhite cci final_route1 final_route6 {
			/*each factor yes: most deprived vs. most affluent*/
			preserve
			rcsgen, scalar(`k') knots(${Kage}) rmatrix(Mage) gen(v)
				
			predictms, transmatrix(tmat) models(m1 m2 m3) ///
			probability los timevar(tt) ci diff ///
			at1(`f' 1 ages1 `=v1' ages2 `=v2' ages3 `=v3')  ///
			at2(`f' 1 income2015_55 1 ages1 `=v1' ages2 `=v2' ages3 `=v3') 
			
			keep _*prob* *_los* tt
			drop if tt ==.
			gen stage = `i'
			gen cancer = "`c'"
			gen age = `k'
			gen analysis = "`f'"
			append using `dif'
			save `dif', replace
			restore
			*/
			}
		}

	
/*	
graph combine "`c'1" "`c'2" "`c'3" "`c'4" "`c'5", row(5) xsize(4.125) ysize(11.75) imargin(zero) title("`c'")
graph save Graph "`c'.gph", replace
*/
}
/*
graph combine "colon" "rectal", col(2) xsize(8.25) ysize(11.75) imargin(zero)
graph export "$res\msbox_1516.svg", as(svg) replace
*/
display "$S_TIME  $S_DATE"

use `hr', clear
save "$res/hr_1516.dta", replace

use `dif', clear
save "$res/dif_1516.dta", replace
display "$S_TIME  $S_DATE"

log close

cap log close
log using "$data\Results\logfiles\crc_treatment_analysis1_multistate_by_screening_R1.log", text replace
display "$S_TIME  $S_DATE"
set seed 1521462
/*S Ling Nov 2021, data analysis for CRC treatment access, multistate*/
///*complete-case analysis*/
///*update on 3rd May 2022: adjust for the same set of covariates for every stage*/
preserve
clear
tempfile hr
save `hr', emptyok replace
restore

preserve
clear
tempfile df
save `df', emptyok replace
restore

preserve
clear
tempfile dif
save `dif', emptyok replace
restore

foreach c in colon rectal {
	use `c', clear
	tab final_route, m
	egen float missing = rowmiss(tnm_stage age female nonwhite final_route income2015)
	tab missing, m
	keep if missing == 0
	gen screening = 1 if final_route == 5
	replace screening = 0 if final_route !=5 & final_route !=.
	gen dead = 1 if death == 1 & (deathdate <= diagdate + 365.24 & deathdate!=.) ///
	| (enddate <= diagdate + 365.24 & enddate!=.)
	replace dead = 0 if dead ==.
	gen censoring = diagdate + 365.24
	egen censor_date = rowmin(censoring enddate deathdate)
	egen treat_date = rowmin(surgdate chemodate radiodate) if treat == 1
	replace treat_date = censor_date if treat == 0 
	gen r1 = runiform(0.1,0.9)
	gen r2 = runiform(0.1,0.5)
	gen treat_time = treat_date - diagdate
	replace treat_time = treat_time + r1 if ( treat_time == 0 & treat == 1)
	gen dead_date = censor_date
	gen dead_time = dead_date - diagdate
	replace dead_time = dead_time + r1 if (dead_time == 0 & dead == 1)
	tab tnm_stage if dead_date<= treat_date & treat == 1 & dead == 1
	replace dead_time = treat_time + r2 if (dead_time == treat_time & dead == 1 & treat == 1)
	tab tnm_stage if dead_date<= treat_date & treat == 1 & dead == 1
	replace treat_time = dead_time if (dead_date == treat_date & dead == 1 & treat == 0)
	format *date %td
	
	gen dead_before_treat = 1 if dead == 1 & treat == 0
	replace dead_before_treat = 0 if dead == 1 & treat == 1
	replace dead_before_treat = 2 if dead_before_treat ==.
	label define ddd 0 "Died after some treatment" 1 "Died before any treatment" 2 "Survived at 1 year after diagnosis"
	label values dead_before_treat ddd
	bysort screening: tab dead_before_treat tnm_stage, m
	replace tnm_stage = 1 if tnm_stage <=3
	replace tnm_stage =2 if tnm_stage == 4
	replace tnm_stage = 3 if tnm_stage == 5
	label define stage 1 "I-III" 2 "IV" 3 "Missing", modify
	bysort screening: tab dead_before_treat tnm_stage, m

	forvalues s = 0(1)1 {
	preserve
	drop if tnm_stage == 3   /*only in stage I-IV*/
	keep if screening == `s'
	tempfile `c'`s'
	save ``c'`s'', replace
	restore
	}

forvalues s = 0(1)1 {
	forvalues i = 1(1)2 {
		use ``c'`s'', clear
		keep if tnm_stage == `i'
		msset, id(pseudo_patientid) states(treat dead) times(treat_time dead_time)
		matrix tmat = r(transmatrix)
		mat list tmat
			stset _stop, enter(_start) failure(_status==1)
			forvalues k = 1(1)3 {
				forvalues j = 1(1)5 {
					cap noisily stpm2 if _trans`k' == 1, df(`j')  scale(h)
					if _rc == 0 {
						preserve
						estat ic
						return list
						matrix m = r(S)
						clear
						svmat long m, names(matcol) 
						gen cancer = "`c'"
						gen stage = `i'
						gen trans = `k'
						gen screening = `s'
						gen df = `j'
						append using `df'
						save `df', replace
						restore
					}
				}
			}
		}
	}
}

use `df', clear
sort screening cancer stage trans df
bysort screening cancer stage trans: gen n_test = _N
bysort screening cancer stage trans: egen max_df = max(df)
drop if max_df != n_test & max_df == df
bysort screening cancer stage trans: egen minAIC = min(mAIC)
bysort screening cancer stage trans: egen minBIC = min(mBIC)
save "$res/df_screening.dta", replace
keep if minAIC == mAIC | minBIC == mBIC
bysort screening cancer stage trans: egen mindf = min(df)
keep if df == mindf
save `df', replace
drop mll0 mll mdf n_test max_df minAIC minBIC mindf
rename m* *
reshape wide N AIC BIC df, i(screening cancer stage) j(trans)
replace cancer = proper(cancer)
label define stage 1 "I-III" 2 "IV" 3 "Missing"
label values stage stage
export excel using "$res\tables.xlsx", sheet("TableS5_df_screen") sheetmodify firstrow(variables)

foreach c in colon rectal {
	forvalues s = 0(1)1 {
		forvalues i = 1(1)2 {
				display "`c' cancer" " stage `i' screening `s'"
				use ``c'`s'', clear
				keep if tnm_stage == `i'
			
				if `s' == 0 {
					foreach var of varlist income2015_5 final_route  {
						tab `var', gen(`var')
					}
				}
			
				if `s' == 1 {
					foreach var of varlist income2015_5  {
						tab `var', gen(`var')
					}	
				}
				rcsgen age, gen(ages) orthog knots(5 35 65 95)
				global Kage `r(knots)'
				matrix Mage = r(R)

				msset, id(pseudo_patientid) states(treat dead) times(treat_time dead_time)
				matrix tmat = r(transmatrix)
				mat list tmat
			/*
			msboxes, transmat(tmat) id(pseudo_patientid) ///
			xvalues(0.2 0.7 0.45) yvalues(0.7 0.7 0.2) ///
			statenames("Diagnosis" "Treatment" "Dead") ///
			boxwidth(0.3)
			graph save Graph "`c'`i'.gph", replace
			*/
			stset _stop, enter(_start) failure(_status==1)
			forvalues k = 1(1)3 {
				
				preserve
				use `df', clear
				sum df if screening == `s' & cancer == "`c'" & stage == `i' & trans == `k'
				local df1 = `r(mean)'
				restore
				
				if `s' == 0 {
				stmerlin income2015_52 income2015_53 income2015_54 income2015_55 ///
				ages1 ages2 ages3 female hf_p mi_p dmc_p pulmonary_p nonwhite  ///
				final_route1 final_route3 final_route4 final_route5  ///
				if _trans`k' == 1, distribution(rp)	df(`df1')
			
				estimates store m`k'
				}
				
				if `s' == 1 {
				stmerlin income2015_52 income2015_53 income2015_54 income2015_55 ///
				ages1 ages2 ages3 female hf_p mi_p dmc_p pulmonary_p nonwhite  ///
				if _trans`k' == 1, distribution(rp)	df(`df1')
			
				estimates store m`k'
				}
				
				preserve
				parmest, eform fast
				gen screening = `s'
				gen stage = `i'
				gen cancer = "`c'"
				gen trans = `k'
				append using `hr'
				save `hr', replace
				restore
				}
			
			forvalues k = 75(5)75 { 
			display "`c' cancer" " stage `i' screening `s'"
			
			cap drop tt
			range tt 0 365 366
			preserve
			rcsgen, scalar(`k') knots(${Kage}) rmatrix(Mage) gen(v)
				
			predictms, transmatrix(tmat) models(m1 m2 m3) ///
			probability timevar(tt) los diff ci ///
			at1(ages1 `=v1' ages2 `=v2' ages3 `=v3')  ///
			at2(income2015_55 1 ages1 `=v1' ages2 `=v2' ages3 `=v3') 
			
			keep _*prob* _*los* tt
			drop if tt ==.
			gen screening = `s'
			gen stage = `i'
			gen cancer = "`c'"
			gen age = `k'
			append using `dif'
			save `dif', replace
			restore
	/*		
			cap drop tt
			range tt 0 360 13
			foreach f in female nonwhite cci final_route1 final_route6 {
				/*each factor yes: most deprived vs. most affluent*/
				preserve
				rcsgen, scalar(`k') knots(${Kage}) rmatrix(Mage) gen(v)
					
				predictms, transmatrix(tmat) models(m1 m2 m3) ///
				probability los timevar(tt) ci diff ///
				at1(`f' 1 ages1 `=v1' ages2 `=v2' ages3 `=v3')  ///
				at2(`f' 1 income2015_55 1 ages1 `=v1' ages2 `=v2' ages3 `=v3') 
				
				keep _*prob* *_los* tt
				drop if tt ==.
				gen stage = `i'
				gen cancer = "`c'"
				gen age = `k'
				gen analysis = "`f'"
				append using `dif'
				save `dif', replace
				restore
				*/
				}
			}

		
	/*	
	graph combine "`c'1" "`c'2" "`c'3" "`c'4" "`c'5", row(5) xsize(4.125) ysize(11.75) imargin(zero) title("`c'")
	graph save Graph "`c'.gph", replace
	*/
	}
}

/*
graph combine "colon" "rectal", col(2) xsize(8.25) ysize(11.75) imargin(zero)
graph export "$res\msbox_1516.svg", as(svg) replace
*/
display "$S_TIME  $S_DATE"

use `hr', clear
save "$res/hr_screening.dta", replace

use `dif', clear
save "$res/dif_screening.dta", replace
display "$S_TIME  $S_DATE"

log close

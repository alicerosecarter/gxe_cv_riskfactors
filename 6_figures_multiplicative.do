/*

Alice R Carter

18/06/2021

Sub do file for figures on the additive scale
*/

* Set local variables and display args

local threshold = "`1'"
local date = "`2'"

di "`threshold'"
di "`date'"

* change to data directory
cd "$resDir/data/results/genetic/202108/"

* Load in results data

foreach var in sd_grs_alcohol_`threshold' sd_grs_af_`threshold' sd_grs_bmi_`threshold'  sd_grs_cad_`threshold' sd_grs_diabetes_t2_`threshold' sd_grs_ldl_`threshold' sd_grs_stroke_`threshold' sd_grs_smoking_`threshold' sd_grs_sbp_`threshold' {

import excel "`date'_p`threshold'_trait_multiplicative", sheet("`var'") firstrow clear
tempfile `var'
save "`var'_multiplicative", replace

}

use sd_grs_alcohol_`threshold'_multiplicative
foreach var in  sd_grs_af_`threshold'_multiplicative sd_grs_bmi_`threshold'_multiplicative  sd_grs_cad_`threshold'_multiplicative sd_grs_diabetes_t2_`threshold'_multiplicative sd_grs_ldl_`threshold'_multiplicative sd_grs_stroke_`threshold'_multiplicative  sd_grs_smoking_`threshold'_multiplicative sd_grs_sbp_`threshold'_multiplicative {

append using `var'


}

* Set figures directory
cd "$resDir/data/results/genetic/figures/202108/"

* Format and label data
rename PValue p_value
lab var p_value "P Value"

rename SampleSize n
lab var n "N"

rename Exposure exposure_label
gen exposure = 1 if exposure_label=="Drinks per week P=5x10-8 (SD)"
replace exposure = 2 if exposure_label=="BMI P=5x10-8 (SD)"
replace exposure = 3 if exposure_label=="LDL-C P=5x10-8 (SD)"
replace exposure = 4 if exposure_label=="Smoking combined P=5x10-8 (SD)"
replace exposure = 5 if exposure_label=="SBP combined P=5x10-8 (SD)"

replace exposure = 6 if exposure_label=="Atrial fibrillation P=5x10-8 (SD)"
replace exposure = 7 if exposure_label=="CAD P=5x10-8 (SD)"
replace exposure = 8 if exposure_label=="Diabetes (T2) P=5x10-8 (SD)"
replace exposure = 9 if exposure_label=="Stroke P=5x10-8 (SD)"

gen scale = 1 if exposure <=5
replace scale = 2 if exposure >=6
lab def scale 1 "Mean difference" 2 "Risk difference", modify
lab val scale scale

rename Interaction interaction
gen eduyears = 1 if interaction=="NONE"
replace eduyears = 7 if interaction=="Education 7"
replace eduyears = 10 if interaction=="Education 10" 
replace eduyears = 13 if interaction=="Education 13"
replace eduyears = 15 if interaction=="Education 15"
replace eduyears = 19 if interaction=="Education 19"
replace eduyears = 20 if interaction=="Education 20"
replace eduyears = 0 if interaction=="Interaction"

lab def eduyears 1 "No interaction" 7 "7 years" 10 "10 years" 13 "13 years" 15 "15 years" 19 "19 years" 20 "20 years" 0 "Interaction coef"
lab val eduyears eduyears
lab var eduyears "Years of education"

rename Outcome outcome
rename BetaOR beta
rename LCI lci
rename UCI uci
rename Numberofcases n_case

gen p_value_int = p_value if eduyears==0

format n %-9.0gc
format beta %8.4f
format lci %8.4f
format uci %8.4f
format exposure  %-8.0f
format eduyears  %-8.0f
format p_value_int %-4.0g 
format p_value %-4.0g 

preserve
drop if scale==2
gen row = _n	


twoway (rspike lci uci row if eduyears==1, lcolor(erose) vertical) ///
		(scatter beta row if eduyears==1, mcolor(erose) msize(.8) msymbol(S))  ///
		(rspike lci uci row if eduyears==7, lcolor(navy*0.5) vertical) ///
		(scatter  beta row  if eduyears==7, mcolor(navy*0.5) msize(.8)) ///
		(rspike lci uci row if eduyears==10, lcolor(navy*0.6) vertical) ///
		(scatter  beta row if eduyears==10, mcolor(navy*0.6) msize(.8)) ///
		(rspike lci uci row if eduyears==13, lcolor(navy*0.7) vertical) ///
		(scatter  beta row if eduyears==13, mcolor(navy*0.7) msize(.8)) ///
		(rspike lci uci row if eduyears==15, lcolor(navy*0.8) vertical) ///
		(scatter  beta row if eduyears==15, mcolor(navy*0.8) msize(.8)) ///
		(rspike lci uci row if eduyears==19, lcolor(navy*0.9) vertical) ///
		(scatter  beta row if eduyears==19, mcolor(navy*0.9) msize(.8)) ///
		(rspike lci uci row if eduyears==20, lcolor(navy) vertical) ///
		(scatter  beta row if eduyears==20, mcolor(navy) msize(.8)), ///		
		legend(row(1) order(2 "No interaction" 4 "7 years" 6 "10 years" 8 "13 years" 10 "15 years" 12 "19 years" 14 "20 years")) ///
		legend(size(vsmall)) ///
		ylabel( -0.05 "-0.05" 0 "0.00" 0.05 "0.05" 0.1 "0.10" 0.15 "0.15" 0.2 "0.20" 0.25 "0.25" 0.3 "0.30", angle(0) labsize(vsmall)) ///
		xlabel(4 "Alcohol"  12 "BMI" 20  "LDL-C" 28 "Smoking" 36  "SBP", angle(0) noticks labsize(vsmall)) ///
		xtitle("Continuous Trait",color(black) size(small)) ytitle("Mean difference in natural log of phenotypic trait",color(black) size(small)) ///
		aspect(.5) graphregion(color(white))

graph export PGS_`threshold'_continuous_multiplicative_`date'.eps, replace

restore, preserve
keep if scale==2

gen row = _n


twoway (rspike lci uci row if eduyears==1, lcolor(erose) vertical) ///
		(scatter beta row if eduyears==1, mcolor(erose) msize(.8) msymbol(S))  ///
		(rspike lci uci row if eduyears==7, lcolor(navy*0.5) vertical) ///
		(scatter  beta row  if eduyears==7, mcolor(navy*0.5) msize(.8)) ///
		(rspike lci uci row if eduyears==10, lcolor(navy*0.6) vertical) ///
		(scatter  beta row if eduyears==10, mcolor(navy*0.6) msize(.8)) ///
		(rspike lci uci row if eduyears==13, lcolor(navy*0.7) vertical) ///
		(scatter  beta row if eduyears==13, mcolor(navy*0.7) msize(.8)) ///
		(rspike lci uci row if eduyears==15, lcolor(navy*0.8) vertical) ///
		(scatter  beta row if eduyears==15, mcolor(navy*0.8) msize(.8)) ///
		(rspike lci uci row if eduyears==19, lcolor(navy*0.9) vertical) ///
		(scatter  beta row if eduyears==19, mcolor(navy*0.9) msize(.8)) ///
		(rspike lci uci row if eduyears==20, lcolor(navy) vertical) ///
		(scatter  beta row if eduyears==20, mcolor(navy) msize(.8)) , ///
		legend(row(1) order(2 "No interaction" 4 "7 years" 6 "10 years" 8 "13 years" 10 "15 years" 12 "19 years" 14 "20 years")) ///
		legend(size(vsmall)) ///
		ylabel( 0.8 "0.8" 1.0 "1.0" 1.2 "1.2" 1.4 "1.4" 1.6 "1.6"1.8 "1.8"  , angle(0) labsize(vsmall)) ///
		xlabel(4.5 "AF" 12 "CHD" 20 "T2D" 27.5 "Stroke", angle(0) noticks labsize(vsmall)) ///
		xtitle("Binary Trait",color(black) size(small)) ytitle("Odds ratio of outcome",color(black) size(small)) ///
		aspect(.5) graphregion(color(white))

graph export PGS_`threshold'_binary_multiplicative_`date'.eps, replace

restore

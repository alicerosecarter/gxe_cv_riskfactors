/*

Alice R Carter

18/06/2021

Sub do file for figures for interaction coefficients
*/
* Set local variables and display args

local threshold = "`1'"
local date = "`2'"

di "`threshold'"
di "`date'"

* change to data directory
cd "$resDir/data/results/genetic/202108/"


use sd_grs_alcohol_`threshold'_additive, clear
foreach var in sd_grs_alcohol_`threshold'_multiplicative sd_grs_af_`threshold'_additive sd_grs_bmi_`threshold'_additive  sd_grs_cad_`threshold'_additive sd_grs_diabetes_t2_`threshold'_additive sd_grs_ldl_`threshold'_additive sd_grs_stroke_`threshold'_additive  sd_grs_smoking_`threshold'_additive sd_grs_sbp_`threshold'_additive sd_grs_af_`threshold'_multiplicative sd_grs_bmi_`threshold'_multiplicative  sd_grs_cad_`threshold'_multiplicative sd_grs_diabetes_t2_`threshold'_multiplicative sd_grs_ldl_`threshold'_multiplicative sd_grs_stroke_`threshold'_multiplicative  sd_grs_smoking_`threshold'_multiplicative sd_grs_sbp_`threshold'_multiplicative {

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

gen scale = 1 if Scale =="additive"
replace scale = 2 if scale==.
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

keep if eduyears==0

sort exposure scale
gen row = _n
gen row_2 = _n-10

twoway (rspike lci uci row if scale==1 & exposure <6 & eduyears==0, lcolor(erose) vertical) ///
		(scatter beta row if scale==1 & exposure <6 & eduyears==0, mcolor(erose) msize(.8) msymbol(S))  ///
		(rspike lci uci row if scale==2 & exposure <6 & eduyears==0, lcolor(navy) vertical) ///
		(scatter beta row if scale==2 & exposure <6 & eduyears==0, mcolor(navy) msize(.8)),  ///
				legend(row(1) order(2 "Additive Interaction" 4 "Multiplicative interaction")) ///
		legend(size(vsmall)) ///
		ylabel(-0.002 "-0.002" -0.001 "-0.001" 0 "0.000"  0.001 "0.001"  0.002 "0.002"  0.003 "0.003"  0.004 "0.004", angle(0) labsize(vsmall)) ///
		xlabel(1.5 "Alcohol"  3.5 "BMI"  5.5 "LDL-C" 7.5 "Smoking" 9.5  "SBP", angle(0) noticks labsize(vsmall)) ///
		xtitle("Continuous Trait",color(black) size(small)) ytitle("Interaction coefficient",color(black) size(small))  ///
		aspect(.5) graphregion(color(white))
		
graph export PGS_`threshold'_interaction_coef_continuous_`date'.eps, replace


twoway (rspike lci uci row_2 if scale==1 & exposure >=6 , lcolor(erose) vertical) ///
		(scatter beta row_2 if scale==1 & exposure >=6, mcolor(erose) msize(0.8) msymbol(S))  ///
		(rspike lci uci row_2 if scale==2 & exposure >=6 , lcolor(navy) vertical) ///
		(scatter beta row_2 if scale==2 & exposure >=6, mcolor(navy) msize(0.8)),  ///
				legend(row(1) order(2 "Additive Interaction" 4 "Multiplicative interaction")) ///
		legend(size(vsmall)) ///
		ylabel(-0.008 "-0.008"  -0.004 "-0.004" 0 "0.000"  0.004 "0.004" 0.008 "0.008" 0.012 "0.012", angle(0) labsize(vsmall)) ///
		xlabel(1.5 "AF"  3.5 "CHD"  5.5 "T2D" 7.5 "Stroke" , angle(0) noticks labsize(vsmall)) ///
		xtitle("Binary Trait",color(black) size(small)) ytitle("Interaction coefficient",color(black) size(small))  ///
		aspect(.5) graphregion(color(white))

graph export PGS_`threshold'_interaction_coef_binary_`date'.eps, replace
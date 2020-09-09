
mark touse
markout touse eduyears statin grs_*
keep if touse==1
drop if CVD_to_exclude==1

global all_grs grs_*
global PCs PC_1-PC_40

foreach var of varlist $all_grs {
	
	egen sd_`var' = std(`var')
}

global pheno_risk drinks_per_week af bmi ckd t2d IHD ldl_c mdd migraine ra Stroke lupus total_chol
global primary_grs sd_grs_alcohol_5e8 sd_grs_af_5e8 sd_grs_bmi_5e8 sd_grs_ckd_5e8 sd_grs_cad_5e8 sd_grs_diabetes_t2_5e8 sd_grs_ldl_5e8 sd_grs_mdd_5e8 sd_grs_migraine_5e8 sd_grs_ra_5e8 sd_grs_sle_5e8 sd_grs_stroke_5e8 sd_grs_total_chol_5e8 

lab var sd_grs_alcohol_5e8 "Drinks per week P=5x10-8 (SD)"
lab var sd_grs_af_5e8 "Atrial fibrillation P=5x10-8 (SD)"
lab var sd_grs_bmi_5e8 "BMI P=5x10-8 (SD)"
lab var sd_grs_cad_5e8 "CAD P=5x10-8 (SD)"
lab var sd_grs_ckd_5e8 "CKD P=5x10-8 (SD)"
lab var sd_grs_diabetes_t2_5e8 "Diabetes (T2) P=5x10-8 (SD)"
lab var sd_grs_ldl_5e8 "LDL-C P=5x10-8 (SD)"
lab var sd_grs_mdd_5e8 "Major Depressive Disorder P=5x10-8 (SD)"
lab var sd_grs_ra_5e8 "Rheumatoid Arthritis P=5x10-8 (SD)"
lab var sd_grs_sle_5e8 "Systemic Lupus Erythematosus P=5x10-8 (SD)"
lab var sd_grs_stroke_5e8 "Stroke P=5x10-8 (SD)"
lab var sd_grs_total_chol_5e8 "Total Cholesterol P=5x10-8 (SD)"
lab var sd_grs_migraine_5e8 "Migraine P=5x10-8 (SD)"

foreach var of varlist $primary_grs {
	
	gen `var'_int = `var'*eduyears
}

rename drinks_per_week alcohol
lab var alcohol "drinks per week"

rename t2_inc diabetes_t2
rename IHD_inc cad
rename Stroke_inc stroke
rename ldl_c ldl
rename csi smoking
rename af_inc af
rename lupus_inc sle
rename ckd_inc ckd
rename mdd_inc mdd
rename arthritis_inc ra
rename migraine_inc migraine

global pheno_risk alcohol af bmi cad ckd diabetes_t2 ldl mdd migraine ra sle stroke total_chol smoking sbp

foreach var of varlist $pheno_risk {
	
	mark touse_sd_grs_`var'_5e8
	markout touse_sd_grs_`var'_5e8 eduyears statin age sex `var'
}

global cont_trait alcohol bmi ldl total_chol smoking sbp

foreach var of varlist $cont_trait {

egen sd_`var' = std(`var')

}


********************************************************************************
							* Setting up results file *
********************************************************************************

foreach exp in $primary_grs {
	putexcel set 202008_p5e8_trait_multiplicative, sheet(`exp') modify
	putexcel A1="Exposure" B1="Outcome" C1="Interaction" D1="Scale" ///
		E1="Beta/OR" F1="LCI" G1="UCI" H1="P Value" I1 = "P value for interaction (with EA)" ///
		J1 = "Sample size"
		
}
********************************************************************************
						* Assessing the role of genetic risk scores *
********************************************************************************


********************************************************************************
*Per unit increase in GRS for continuous traits

** Not including split sample variables, SBP and CSI - will run separately in script later

local pheno_out sd_alcohol sd_bmi sd_ldl sd_total_chol 
local primary_grs sd_grs_alcohol_5e8 sd_grs_bmi_5e8 sd_grs_ldl_5e8 sd_grs_total_chol_5e8 

foreach var of varlist `pheno_out' {

	gen ln_`var' = ln(`var')
	
}	

local pheno_out ln_sd_alcohol ln_sd_bmi ln_sd_ldl ln_sd_total_chol 

local n : word count `pheno_out'

forvalues i = 1/`n' {

	local out : word `i' of `pheno_out'
	local exp : word `i' of `primary_grs'

	
	local x=1
	local x=`x'+1
	
	putexcel set 202008_p5e8_trait_multiplicative, sheet(`exp') modify	
		
		regress `out' `exp' eduyears age sex $PCs if touse_`exp'==1
	
	matrix results = r(table)

	local beta = results[1,1]
	local lci = results[5,1]
	local uci = results[6,1]
	local p_value = results[4,1]
	local out_label : var label `out'
	local exp_label : var label `exp'
	local n = e(N)
	
	putexcel A`x'="`exp_label'" B`x'="`out_label'" C`x'="NONE" D`x'="Differnce" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_value' J`x'=`n'
	
	
	}

********************************************************************************
* Per year of educational attainment
local pheno_out ln_sd_alcohol ln_sd_bmi ln_sd_ldl ln_sd_total_chol 
local primary_grs sd_grs_alcohol_5e8 sd_grs_bmi_5e8 sd_grs_ldl_5e8 sd_grs_total_chol_5e8 
local n : word count `pheno_out'
levelsof eduyears, local(levels)

forvalues i = 1/`n' {
	local out : word `i' of `pheno_out'
	local exp : word `i' of `primary_grs'
local x=2
foreach l of local levels {

	local x=`x'+1
		

	putexcel set 202008_p5e8_trait_multiplicative, sheet(`exp') modify	

regress `out' `exp' age sex  $PCs if eduyears==`l' & touse_`exp'==1

	matrix results = r(table)

	local beta = results[1,1]
	local lci = results[5,1]
	local uci = results[6,1]
	local p_value = results[4,1]
	local out_label : var label `out'
	local exp_label : var label `exp'
	local n = e(N)
	
	putexcel A`x'="`exp_label'" B`x'="`out_label'" C`x'="Education `l'" D`x'="Difference" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_value' J`x'=`n'
	


}		
}		

********************************************************************************
** Testing for interaction with education

local pheno_out ln_sd_alcohol ln_sd_bmi ln_sd_ldl ln_sd_total_chol 
local primary_grs sd_grs_alcohol_5e8 sd_grs_bmi_5e8 sd_grs_ldl_5e8 sd_grs_total_chol_5e8 
local n : word count `pheno_out'

forvalues i = 1/`n' {
	local out : word `i' of `pheno_out'
	local exp : word `i' of `primary_grs'
local x=9


	putexcel set 202008_p5e8_trait_multiplicative, sheet(`exp') modify	

regress `out' `exp' eduyears `exp'_int age sex  $PCs if touse_`exp'==1

	matrix results = r(table)
	local beta_int = results[1,3]
	local lci_int = results[5,3]
	local uci_int = results[6,3]
	local p_value = results[4,3]
	local out_label : var label `out'
	local exp_label : var label `exp'
	
putexcel A`x'="`exp_label'" B`x'="`out_label'" C`x'="Interaction" D`x'="Multiplicative" E`x'=`beta_int' F`x'=`lci_int' G`x'=`uci_int' H`x'=`p_value' 


}	

** Smoking - split sample analysis_sample

gen ln_sd_smoking = ln(sd_smoking)

foreach out in ln_sd_smoking {
local x = 1
local x = `x'+1

	
	putexcel set csi_results_trait_5e8_multiplicative, sheet(all) modify
	putexcel A1="sample" B1="out" C1="log_or" D1="lci" E1="uci" F1="n" 

	
	regress `out' sd_grs_smoking_5e8_gwas2 eduyears age sex $PCs if sample==1 & touse_sd_grs_smoking_5e8==1
	
	local out_label : var label `out'
	local lci = _b[sd_grs_smoking_5e8_gwas2]-1.96*_se[sd_grs_smoking_5e8_gwas2]
	local uci = _b[sd_grs_smoking_5e8_gwas2]+1.96*_se[sd_grs_smoking_5e8_gwas2]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[sd_grs_smoking_5e8_gwas2] D`x'=`lci' E`x'=`uci' 
	
	regress `out' sd_grs_smoking_5e8_gwas1 eduyears age sex $PCs if sample==2 & touse_sd_grs_smoking_5e8==1

	
	local x = `x'+1

	local out_label : var label `out'
	local lci = _b[sd_grs_smoking_5e8_gwas1]-1.96*_se[sd_grs_smoking_5e8_gwas1]
	local uci = _b[sd_grs_smoking_5e8_gwas1]+1.96*_se[sd_grs_smoking_5e8_gwas1]
	
	putexcel A`x'="2" B`x'="`out_label'" C`x'=_b[sd_grs_smoking_5e8_gwas1] D`x'=`lci' E`x'=`uci'  
		
	
}	

levelsof eduyears, local(levels)
local x =3	
foreach out in ln_sd_smoking  {


foreach l of local levels {

	
	putexcel set csi_results_trait_5e8_multiplicative, sheet(all) modify
	putexcel A1="sample" B1="out" C1="or" D1="lci" E1="uci" F1="n" G1="eduyears_csi"

	
	regress `out' sd_grs_smoking_5e8_gwas2 age sex $PCs if eduyears==`l' & sample==1 & touse_sd_grs_smoking_5e8==1
	
	local x=`x'+1
	local out_label : var label `out'
	local lci = _b[sd_grs_smoking_5e8_gwas2]-1.96*_se[sd_grs_smoking_5e8_gwas2]
	local uci = _b[sd_grs_smoking_5e8_gwas2]+1.96*_se[sd_grs_smoking_5e8_gwas2]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[sd_grs_smoking_5e8_gwas2] D`x'=`lci' E`x'=`uci' G`x'="Education `l'"
	
	regress `out' sd_grs_smoking_5e8_gwas1 age sex $PCs if eduyears==`l' & sample==2 & touse_sd_grs_smoking_5e8==1

	local x=`x'+1
	
	local out_label : var label `out'
	local lci = _b[sd_grs_smoking_5e8_gwas1]-1.96*_se[sd_grs_smoking_5e8_gwas1]
	local uci = _b[sd_grs_smoking_5e8_gwas1]+1.96*_se[sd_grs_smoking_5e8_gwas1]
	
	putexcel A`x'="2" B`x'="`out_label'" C`x'=_b[sd_grs_smoking_5e8_gwas1] D`x'=`lci' E`x'=`uci'  G`x'="Education `l'"
		
	
}	
}

** Testing for interaction with education

gen sd_grs_smoking_5e8_gwas2_int = sd_grs_smoking_5e8_gwas2*eduyears if sample==1
gen sd_grs_smoking_5e8_gwas1_int = sd_grs_smoking_5e8_gwas1*eduyears if sample==2


foreach out in ln_sd_smoking  {

	local x = 15
	
	putexcel set csi_results_trait_5e8_multiplicative, sheet(all) modify
	
	local x=`x'+1
	regress `out' sd_grs_smoking_5e8_gwas2 eduyears sd_grs_smoking_5e8_gwas2_int age sex $PCs if touse_sd_grs_smoking_5e8==1 & sample==1
	
	local lci = _b[sd_grs_smoking_5e8_gwas2_int]-1.96*_se[sd_grs_smoking_5e8_gwas2_int]
	local uci = _b[sd_grs_smoking_5e8_gwas2_int]+1.96*_se[sd_grs_smoking_5e8_gwas2_int]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[sd_grs_smoking_5e8_gwas2_int] D`x'=`lci' E`x'=`uci' G`x'="EA interact"

	regress `out' sd_grs_smoking_5e8_gwas1 eduyears sd_grs_smoking_5e8_gwas1_int age sex $PCs if touse_sd_grs_smoking_5e8==1 & sample==2
	
	
	local x=`x'+1
	local lci = _b[sd_grs_smoking_5e8_gwas1_int]-1.96*_se[sd_grs_smoking_5e8_gwas1_int]
	local uci = _b[sd_grs_smoking_5e8_gwas1_int]+1.96*_se[sd_grs_smoking_5e8_gwas1_int]
	
	putexcel A`x'="2" B`x'="`out_label'" C`x'=_b[sd_grs_smoking_5e8_gwas1_int] D`x'=`lci' E`x'=`uci' G`x'="EA interact"

	
}
	
** Systolic blood pressure - split sample analysis_sample
gen ln_sd_sbp = ln(sd_sbp)

foreach out in ln_sd_sbp {
local x = 1
local x = `x'+1

	
	putexcel set sbp_results_trait_5e8_multiplicative, sheet(all) modify
	putexcel A1="sample" B1="out" C1="or" D1="lci" E1="uci" F1="n" 

	
	regress `out' sd_grs_sbp_5e8_gwas2 eduyears age sex $PCs if sample==1 & touse_sd_grs_sbp_5e8==1
	
	local out_label : var label `out'
	local lci = _b[sd_grs_sbp_5e8_gwas2]-1.96*_se[sd_grs_sbp_5e8_gwas2]
	local uci = _b[sd_grs_sbp_5e8_gwas2]+1.96*_se[sd_grs_sbp_5e8_gwas2]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[sd_grs_sbp_5e8_gwas2] D`x'=`lci' E`x'=`uci' 
	
	regress `out' sd_grs_sbp_5e8_gwas1 eduyears age sex $PCs if sample==2 & touse_sd_grs_sbp_5e8==1

	
	local x = `x'+1

	local out_label : var label `out'
	local lci = _b[sd_grs_sbp_5e8_gwas1]-1.96*_se[sd_grs_sbp_5e8_gwas1]
	local uci = _b[sd_grs_sbp_5e8_gwas1]+1.96*_se[sd_grs_sbp_5e8_gwas1]
	
	putexcel A`x'="2" B`x'="`out_label'" C`x'=_b[sd_grs_sbp_5e8_gwas1] D`x'=`lci' E`x'=`uci'  
		
	
}	

levelsof eduyears, local(levels)
local x =3	
foreach out in ln_sd_sbp  {


foreach l of local levels {

	
	putexcel set sbp_results_trait_5e8_multiplicative, sheet(all) modify
	putexcel A1="sample" B1="out" C1="or" D1="lci" E1="uci" F1="n" G1="eduyears_sbp"

	
	regress `out' sd_grs_sbp_5e8_gwas2 age sex $PCs if eduyears==`l' & sample==1 & touse_sd_grs_sbp_5e8==1
	
	local x=`x'+1
	local out_label : var label `out'
	local lci = _b[sd_grs_sbp_5e8_gwas2]-1.96*_se[sd_grs_sbp_5e8_gwas2]
	local uci = _b[sd_grs_sbp_5e8_gwas2]+1.96*_se[sd_grs_sbp_5e8_gwas2]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[sd_grs_sbp_5e8_gwas2] D`x'=`lci' E`x'=`uci' G`x'="Education `l'"
	
	regress `out' sd_grs_sbp_5e8_gwas1 age sex $PCs if eduyears==`l' & sample==2 & touse_sd_grs_sbp_5e8==1

	local x=`x'+1
	
	local out_label : var label `out'
	local lci = _b[sd_grs_sbp_5e8_gwas1]-1.96*_se[sd_grs_sbp_5e8_gwas1]
	local uci = _b[sd_grs_sbp_5e8_gwas1]+1.96*_se[sd_grs_sbp_5e8_gwas1]
	
	putexcel A`x'="2" B`x'="`out_label'" C`x'=_b[sd_grs_sbp_5e8_gwas1] D`x'=`lci' E`x'=`uci'  G`x'="Education `l'"
		
	
}	
}

** Testing for interaction with education

gen sd_grs_sbp_5e8_gwas2_int = sd_grs_sbp_5e8_gwas2*eduyears if sample==1
gen sd_grs_sbp_5e8_gwas1_int = sd_grs_sbp_5e8_gwas1*eduyears if sample==2


foreach out in ln_sd_sbp  {

	local x = 15
	
	putexcel set sbp_results_trait_5e8_multiplicative, sheet(all) modify
	
	local x=`x'+1
	
	regress `out' sd_grs_sbp_5e8_gwas2 eduyears sd_grs_sbp_5e8_gwas2_int age sex $PCs if touse_sd_grs_sbp_5e8==1 & sample==1
	
	local lci = _b[sd_grs_sbp_5e8_gwas2_int]-1.96*_se[sd_grs_sbp_5e8_gwas2_int]
	local uci = _b[sd_grs_sbp_5e8_gwas2_int]+1.96*_se[sd_grs_sbp_5e8_gwas2_int]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[sd_grs_sbp_5e8_gwas2_int] D`x'=`lci' E`x'=`uci' G`x'="EA interact"

	regress `out' sd_grs_sbp_5e8_gwas1 eduyears sd_grs_sbp_5e8_gwas1_int age sex $PCs if touse_sd_grs_sbp_5e8==1 & sample==2
	
	
	local x=`x'+1
	local lci = _b[sd_grs_sbp_5e8_gwas1_int]-1.96*_se[sd_grs_sbp_5e8_gwas1_int]
	local uci = _b[sd_grs_sbp_5e8_gwas1_int]+1.96*_se[sd_grs_sbp_5e8_gwas1_int]
	
	putexcel A`x'="2" B`x'="`out_label'" C`x'=_b[sd_grs_sbp_5e8_gwas1_int] D`x'=`lci' E`x'=`uci' G`x'="EA interact"

	
}

********************************************************************************
* Binary traits - logistic scale for outcome

local pheno_out  af  cad ckd diabetes_t2  mdd migraine ra sle stroke  
local primary_grs  sd_grs_af_5e8   sd_grs_cad_5e8 sd_grs_ckd_5e8 sd_grs_diabetes_t2_5e8  sd_grs_mdd_5e8 sd_grs_migraine_5e8 sd_grs_ra_5e8 sd_grs_sle_5e8 sd_grs_stroke_5e8  
local n : word count `pheno_out'

forvalues i = 1/`n' {

	local out : word `i' of `pheno_out'
	local exp : word `i' of `primary_grs'

	
	local x=1
	local x=`x'+1
	
	putexcel set 202008_p5e8_trait_multiplicative, sheet(`exp') modify	
		
		logistic `out' `exp' eduyears age sex $PCs if touse_`exp'==1
	
	matrix results = r(table)

	local beta = results[1,1]
	local lci = results[5,1]
	local uci = results[6,1]
	local p_value = results[4,1]
	local out_label : var label `out'
	local exp_label : var label `exp'
	local n = e(N)
	
	putexcel A`x'="`exp_label'" B`x'="`out_label'" C`x'="NONE" D`x'="Logistic" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_value' J`x'=`n'
	
	
	}

********************************************************************************
* Per year of educational attainment
local pheno_out  af  cad ckd diabetes_t2  mdd migraine ra sle stroke  
local primary_grs  sd_grs_af_5e8   sd_grs_cad_5e8 sd_grs_ckd_5e8 sd_grs_diabetes_t2_5e8  sd_grs_mdd_5e8 sd_grs_migraine_5e8 sd_grs_ra_5e8 sd_grs_sle_5e8 sd_grs_stroke_5e8  
local n : word count `pheno_out'
levelsof eduyears, local(levels)

forvalues i = 1/`n' {
	local out : word `i' of `pheno_out'
	local exp : word `i' of `primary_grs'
local x=2
foreach l of local levels {

	local x=`x'+1
		

	putexcel set 202008_p5e8_trait_multiplicative, sheet(`exp') modify	

logistic `out' `exp' age sex  $PCs if eduyears==`l' & touse_`exp'==1

	matrix results = r(table)

	local beta = results[1,1]
	local lci = results[5,1]
	local uci = results[6,1]
	local p_value = results[4,1]
	local out_label : var label `out'
	local exp_label : var label `exp'
	local n = e(N)
	
	putexcel A`x'="`exp_label'" B`x'="`out_label'" C`x'="Education `l'" D`x'="Logistic" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_value' J`x'=`n'
	


}		
}		

********************************************************************************
** Testing for interaction with education

local pheno_out  af  cad ckd diabetes_t2  mdd migraine ra sle stroke  
local primary_grs  sd_grs_af_5e8   sd_grs_cad_5e8 sd_grs_ckd_5e8 sd_grs_diabetes_t2_5e8  sd_grs_mdd_5e8 sd_grs_migraine_5e8 sd_grs_ra_5e8 sd_grs_sle_5e8 sd_grs_stroke_5e8  
local n : word count `pheno_out'

forvalues i = 1/`n' {
	local out : word `i' of `pheno_out'
	local exp : word `i' of `primary_grs'
local x=9


	putexcel set 202008_p5e8_trait_multiplicative, sheet(`exp') modify	

logit `out' `exp' eduyears `exp'_int age sex  $PCs if touse_`exp'==1


	matrix results = r(table)
	local beta_int = results[1,3]
	local lci_int = results[5,3]
	local uci_int = results[6,3]
	local p_value = results[4,3]
	local out_label : var label `out'
	local exp_label : var label `exp'
	
putexcel A`x'="`exp_label'" B`x'="`out_label'" C`x'="Interaction" D`x'="Multiplicative" E`x'=`beta_int' F`x'=`lci_int' G`x'=`uci_int' H`x'=`p_value' 

}	


local pheno_out  af  cad ckd diabetes_t2  mdd migraine ra sle stroke  
local primary_grs  sd_grs_af_5e8   sd_grs_cad_5e8 sd_grs_ckd_5e8 sd_grs_diabetes_t2_5e8  sd_grs_mdd_5e8 sd_grs_migraine_5e8 sd_grs_ra_5e8 sd_grs_sle_5e8 sd_grs_stroke_5e8  
local n : word count `pheno_out'

forvalues i = 1/`n' {
	local out : word `i' of `pheno_out'
	local exp : word `i' of `primary_grs'
local x=3


	putexcel set 202008_p5e8_trait_multiplicative, sheet(`exp') modify	

logistic `out' `exp' eduyears `exp'_int age sex  $PCs if touse_`exp'==1


	matrix results = r(table)
	local p_value = results[4,3]
	
	putexcel I`x'=`p_value'

}	
*** Tab sample size 

foreach touse in  touse_sd_grs_smoking_5e8 touse_sd_grs_sbp_5e8  {

		local sheet = substr("`touse'",7,.)
	
		putexcel set 202008_p5e8_trait_multiplicative, sheet(`sheet') modify
		
		tab statin eduyears if `touse'==1, matcell(numbers)
		 
		local n_7 = sum(numbers[2,1]+numbers[1,1])
		local n_10 = sum(numbers[2,2]+numbers[1,2])
		local n_13 = sum(numbers[2,3]+numbers[1,3])
		local n_15 = sum(numbers[2,4]+numbers[1,4])
		local n_19 = sum(numbers[2,5]+numbers[1,5])
		local n_20 = sum(numbers[2,6]+numbers[1,6])
		local n_all = sum(`n_7'+`n_10'+`n_13'+`n_15'+`n_19'+`n_20')
		
	
		putexcel J2=`n_all' J3=`n_7' J4=`n_10' J5=`n_13' J6=`n_15' J7=`n_19' J8=`n_20'

		 
		 }

	
** Meta-analysing two samples - smoking

import excel "csi_results_trait_5e8_multiplicative.xlsx", sheet("all") firstrow clear
destring sample, replace
lab var sample "Sample"
lab def sample 1 "Sample 1" 2 "Sample 2"
lab val sample sample

lab var n "N"

generate eduyears_level=.
replace eduyears_level = 7 if eduyears_csi=="Education 7"
replace eduyears_level = 10 if eduyears_csi=="Education 10"
replace eduyears_level = 13 if eduyears_csi=="Education 13"
replace eduyears_level = 15 if eduyears_csi=="Education 15"
replace eduyears_level = 19 if eduyears_csi=="Education 19"
replace eduyears_level = 20 if eduyears_csi=="Education 20"

gen eduyears_interact = 1 if eduyears_csi=="EA interact"

	metan or lci uci if  eduyears_level==. & eduyears_interact==., lcols(sample ) effect(OR) null(1) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	putexcel set 202008_p5e8_trait_multiplicative, sheet(sd_grs_smoking_5e8) modify
	putexcel A1="Exposure" B1="Outcome" C1="Interaction" D1="Scale" ///
		E1="Beta/OR" F1="LCI" G1="UCI" H1="P Value" I1 = "P value for interaction (with EA)" ///
		J1 = "Sample size" 
	local x = 2
	local beta = (r(ES))
	local lci = (r(ci_low))
	local uci = (r(ci_upp))
	local p_val = r(p_z)	
	putexcel A`x'="Smoking combined P=5x10-8 (SD)" B`x'="Smoking" C`x'="NONE" D`x'="Difference" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_val'

levelsof eduyears_level, local(levels)
local x = 2
foreach l of local levels {

	putexcel set 202008_p5e8_trait_multiplicative, sheet(sd_grs_smoking_5e8) modify


	metan or lci uci if  eduyears_level==`l', lcols(sample ) effect(OR) null(1) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	
	local x = `x'+1
	local beta = (r(ES))
	local lci = (r(ci_low))
	local uci = (r(ci_upp))
	local p_val = r(p_z)	
	putexcel A`x'="Smoking combined P=5x10-8 (SD)" B`x'="Smoking" C`x'="Education `l'" D`x'="Difference" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_val'
	
}

metan or lci uci if  eduyears_interact==1, lcols(sample ) effect(OR) null(1) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	putexcel set 202008_p5e8_trait_multiplicative, sheet(sd_grs_smoking_5e8) modify

	local x=9
	local beta_int = r(ES)
	local lci_int = r(ci_low)
	local uci_int = r(ci_upp)
	local p_value = r(p_z)

	
putexcel A`x'="SMoking combined P=5x10-8 (SD)" B`x'="Smoking" C`x'="Interaction" D`x'="Multiplicative" E`x'=`beta_int' F`x'=`lci_int' G`x'=`uci_int' H`x'=`p_value' 



	** Meta-analysing two samples - SBP

import excel "sbp_results_trait_5e8_multiplicative.xlsx", sheet("all") firstrow clear
destring sample, replace
lab var sample "Sample"
lab def sample 1 "Sample 1" 2 "Sample 2"
lab val sample sample

lab var n "N"


generate eduyears_level=.
replace eduyears_level = 7 if eduyears_sbp=="Education 7"
replace eduyears_level = 10 if eduyears_sbp=="Education 10"
replace eduyears_level = 13 if eduyears_sbp=="Education 13"
replace eduyears_level = 15 if eduyears_sbp=="Education 15"
replace eduyears_level = 19 if eduyears_sbp=="Education 19"
replace eduyears_level = 20 if eduyears_sbp=="Education 20"

gen eduyears_interact = 1 if eduyears_sbp=="EA interact"


	metan or lci uci if  eduyears_level==. & eduyears_interact==., lcols(sample ) effect(OR) null(1) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	putexcel set 202008_p5e8_trait_multiplicative, sheet(sd_grs_sbp_5e8) modify
	putexcel A1="Exposure" B1="Outcome" C1="Interaction" D1="Scale" ///
		E1="Beta/OR" F1="LCI" G1="UCI" H1="P Value" I1 = "P value for interaction (with EA)" ///
		J1 = "Sample size" K1="N for statin use"
	local x = 2
	local beta = (r(ES))
	local lci = (r(ci_low))
	local uci = (r(ci_upp))
	local p_val = r(p_z)	
	putexcel A`x'="SBP combined P=5x10-8 (SD)" B`x'="SBP" C`x'="NONE" D`x'="Difference" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_val'

levelsof eduyears_level, local(levels)
local x = 2
foreach l of local levels {

	putexcel set 202008_p5e8_trait_multiplicative, sheet(sd_grs_sbp_5e8) modify


	metan or lci uci if  eduyears_level==`l', lcols(sample ) effect(OR) null(1) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	
	local x = `x'+1
	local beta = (r(ES))
	local lci = (r(ci_low))
	local uci = (r(ci_upp))
	local p_val = r(p_z)	
	putexcel A`x'="SBP combined P=5x10-8 (SD)" B`x'="SBP" C`x'="Education `l'" D`x'="Difference" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_val'
	
}

metan or lci uci if  eduyears_interact==1, lcols(sample ) effect(OR) null(1) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	putexcel set 202008_p5e8_trait_multiplicative, sheet(sd_grs_sbp_5e8) modify

	local x=9
	local beta_int = (r(ES))
	local lci_int = (r(ci_low))
	local uci_int = (r(ci_upp))
	local p_value = r(p_z)

	
putexcel A`x'="SBP combined P=5x10-8 (SD)" B`x'="SBP (SD)" C`x'="Interaction" D`x'="Multiplicative" E`x'=`beta_int' F`x'=`lci_int' G`x'=`uci_int' H`x'=`p_value' 



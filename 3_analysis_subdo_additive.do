/*

Alice R Carter

18/06/2021

Sub-do file for carrying out analyses on the additive scale
*/

* Set local variables and display args

local threshold = "`1'"
local date = "`2'"

di "`threshold'"
di "`date'"

* Set environment
use "$resDir/data/stata/genetic/GxE_analysis_202108.dta", clear
cd "$resDir/data/results/genetic/202108"

global pheno_risk alcohol af bmi cad diabetes_t2 ldl stroke smoking sbp

foreach var of varlist $pheno_risk {
	
	mark touse_sd_grs_`var'_`threshold'
	markout touse_sd_grs_`var'_`threshold' eduyears age sex `var'
}


********************************************************************************
							* Setting up results file *
********************************************************************************

global grs sd_grs_alcohol_`threshold' sd_grs_af_`threshold' sd_grs_bmi_`threshold'  sd_grs_cad_`threshold' sd_grs_diabetes_t2_`threshold' sd_grs_ldl_`threshold' sd_grs_stroke_`threshold' 

foreach exp in $grs {
	putexcel set `date'_p`threshold'_trait_additive, sheet(`exp') modify
	putexcel A1="Exposure" B1="Outcome" C1="Interaction" D1="Scale" ///
		E1="Beta/OR" F1="LCI" G1="UCI" H1="P Value" I1 = "" ///
		J1 = "Sample Size" K1 = "Number of cases"
		
}
********************************************************************************
						* Assessing the role of genetic risk scores *
********************************************************************************


********************************************************************************
*Per SD increase in GRS

** Not including split sample variables, SBP and CSI - will run separately in script later

local pheno_out sd_alcohol af sd_bmi cad diabetes_t2 sd_ldl stroke  
local primary_grs sd_grs_alcohol_`threshold' sd_grs_af_`threshold' sd_grs_bmi_`threshold'  sd_grs_cad_`threshold' sd_grs_diabetes_t2_`threshold' sd_grs_ldl_`threshold' sd_grs_stroke_`threshold' 

local n : word count `pheno_out'

forvalues i = 1/`n' {

	local out : word `i' of `pheno_out'
	local exp : word `i' of `primary_grs'

	
	local x=1
	local x=`x'+1
	
	putexcel set `date'_p`threshold'_trait_additive, sheet(`exp') modify	
		
		regress `out' `exp' eduyears age sex $PCs `covariate' if touse_`exp'==1
	
	matrix results = r(table)

	local beta = results[1,1]
	local lci = results[5,1]
	local uci = results[6,1]
	local p_value = results[4,1]
	local out_label : var label `out'
	local exp_label : var label `exp'
	local n = e(N)
	
	putexcel A`x'="`exp_label'" B`x'="`out_label'" C`x'="NONE" D`x'="additive" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_value' J`x'=`n'
	
	
	}

********************************************************************************
* Per year of educational attainment
local pheno_out sd_alcohol af sd_bmi cad diabetes_t2 sd_ldl stroke  
local primary_grs sd_grs_alcohol_`threshold' sd_grs_af_`threshold' sd_grs_bmi_`threshold'  sd_grs_cad_`threshold' sd_grs_diabetes_t2_`threshold' sd_grs_ldl_`threshold' sd_grs_stroke_`threshold' 

local n : word count `pheno_out'
levelsof eduyears, local(levels)

forvalues i = 1/`n' {
	local out : word `i' of `pheno_out'
	local exp : word `i' of `primary_grs'
local x=2
foreach l of local levels {

	local x=`x'+1
		

	putexcel set `date'_p`threshold'_trait_additive, sheet(`exp') modify	

regress `out' `exp' age sex  $PCs `covariate' if eduyears==`l' & touse_`exp'==1

	matrix results = r(table)

	local beta = results[1,1]
	local lci = results[5,1]
	local uci = results[6,1]
	local p_value = results[4,1]
	local out_label : var label `out'
	local exp_label : var label `exp'
	local n = e(N)
	
	putexcel A`x'="`exp_label'" B`x'="`out_label'" C`x'="Education `l'" D`x'="additive" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_value' J`x'=`n'
	


}		
}		

********************************************************************************
** Testing for interaction with education

local pheno_out sd_alcohol af sd_bmi cad diabetes_t2 sd_ldl stroke  
local primary_grs sd_grs_alcohol_`threshold' sd_grs_af_`threshold' sd_grs_bmi_`threshold'  sd_grs_cad_`threshold' sd_grs_diabetes_t2_`threshold' sd_grs_ldl_`threshold' sd_grs_stroke_`threshold' 
local n : word count `pheno_out'

forvalues i = 1/`n' {
	local out : word `i' of `pheno_out'
	local exp : word `i' of `primary_grs'
local x=9


	putexcel set `date'_p`threshold'_trait_additive, sheet(`exp') modify	

regress `out' `exp' eduyears `exp'_int age sex  $PCs `covariate' if touse_`exp'==1


	matrix results = r(table)
	local beta_int = results[1,3]
	local lci_int = results[5,3]
	local uci_int = results[6,3]
	local p_value = results[4,3]
	local out_label : var label `out'
	local exp_label : var label `exp'
	
putexcel A`x'="`exp_label'" B`x'="`out_label'" C`x'="Interaction" D`x'="additive" E`x'=`beta_int' F`x'=`lci_int' G`x'=`uci_int' H`x'=`p_value' 

	

}	

** Smoking - split sample analysis_sample

foreach out in sd_smoking {
local x = 1
local x = `x'+1

	
	putexcel set csi_results_trait_`threshold'_additive, sheet(all) modify
	putexcel A1="sample" B1="out" C1="log_or" D1="lci" E1="uci" F1="n" 

	
	regress `out' sd_grs_smoking_`threshold'_gwas2 eduyears age sex $PCs `covariate' if sample==1 & touse_sd_grs_smoking_`threshold'==1
	
	local out_label : var label `out'
	local lci = _b[sd_grs_smoking_`threshold'_gwas2]-1.96*_se[sd_grs_smoking_`threshold'_gwas2]
	local uci = _b[sd_grs_smoking_`threshold'_gwas2]+1.96*_se[sd_grs_smoking_`threshold'_gwas2]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[sd_grs_smoking_`threshold'_gwas2] D`x'=`lci' E`x'=`uci' 
	
	regress `out' sd_grs_smoking_`threshold'_gwas1 eduyears age sex $PCs `covariate' if sample==2 & touse_sd_grs_smoking_`threshold'==1

	
	local x = `x'+1

	local out_label : var label `out'
	local lci = _b[sd_grs_smoking_`threshold'_gwas1]-1.96*_se[sd_grs_smoking_`threshold'_gwas1]
	local uci = _b[sd_grs_smoking_`threshold'_gwas1]+1.96*_se[sd_grs_smoking_`threshold'_gwas1]
	
	putexcel A`x'="2" B`x'="`out_label'" C`x'=_b[sd_grs_smoking_`threshold'_gwas1] D`x'=`lci' E`x'=`uci'  
		
	
}	

levelsof eduyears, local(levels)
local x =3	
foreach out in sd_smoking  {


foreach l of local levels {

	
	putexcel set csi_results_trait_`threshold'_additive, sheet(all) modify
	putexcel A1="sample" B1="out" C1="or" D1="lci" E1="uci" F1="n" G1="eduyears_csi"

	
	regress `out' sd_grs_smoking_`threshold'_gwas2 age sex $PCs  `covariate' if eduyears==`l' & sample==1 & touse_sd_grs_smoking_`threshold'==1
	
	local x=`x'+1
	local out_label : var label `out'
	local lci = _b[sd_grs_smoking_`threshold'_gwas2]-1.96*_se[sd_grs_smoking_`threshold'_gwas2]
	local uci = _b[sd_grs_smoking_`threshold'_gwas2]+1.96*_se[sd_grs_smoking_`threshold'_gwas2]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[sd_grs_smoking_`threshold'_gwas2] D`x'=`lci' E`x'=`uci' G`x'="Education `l'"
	
	regress `out' sd_grs_smoking_`threshold'_gwas1 age sex $PCs `covariate' if eduyears==`l' & sample==2 & touse_sd_grs_smoking_`threshold'==1

	local x=`x'+1
	
	local out_label : var label `out'
	local lci = _b[sd_grs_smoking_`threshold'_gwas1]-1.96*_se[sd_grs_smoking_`threshold'_gwas1]
	local uci = _b[sd_grs_smoking_`threshold'_gwas1]+1.96*_se[sd_grs_smoking_`threshold'_gwas1]
	
	putexcel A`x'="2" B`x'="`out_label'" C`x'=_b[sd_grs_smoking_`threshold'_gwas1] D`x'=`lci' E`x'=`uci'  G`x'="Education `l'"
		
	
}	
}

** Testing for interaction with education

foreach out in sd_smoking  {

	local x = 15
	
	putexcel set csi_results_trait_`threshold'_additive, sheet(all) modify
	
	local x=`x'+1
	regress `out' sd_grs_smoking_`threshold'_gwas2 eduyears sd_grs_smoking_`threshold'_gwas2_int age sex $PCs `covariate' if touse_sd_grs_smoking_`threshold'==1 & sample==1
	
	local lci = _b[sd_grs_smoking_`threshold'_gwas2_int]-1.96*_se[sd_grs_smoking_`threshold'_gwas2_int]
	local uci = _b[sd_grs_smoking_`threshold'_gwas2_int]+1.96*_se[sd_grs_smoking_`threshold'_gwas2_int]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[sd_grs_smoking_`threshold'_gwas2_int] D`x'=`lci' E`x'=`uci' G`x'="EA interact"

	regress `out' sd_grs_smoking_`threshold'_gwas1 eduyears sd_grs_smoking_`threshold'_gwas1_int age sex $PCs `covariate' if touse_sd_grs_smoking_`threshold'==1 & sample==2
	
	
	local x=`x'+1
	local lci = _b[sd_grs_smoking_`threshold'_gwas1_int]-1.96*_se[sd_grs_smoking_`threshold'_gwas1_int]
	local uci = _b[sd_grs_smoking_`threshold'_gwas1_int]+1.96*_se[sd_grs_smoking_`threshold'_gwas1_int]
	
	putexcel A`x'="2" B`x'="`out_label'" C`x'=_b[sd_grs_smoking_`threshold'_gwas1_int] D`x'=`lci' E`x'=`uci' G`x'="EA interact"

	
}


	
** Systolic blood pressure - split sample analysis_sample

foreach out in sd_sbp {
local x = 1
local x = `x'+1

	
	putexcel set sbp_results_trait_`threshold'_additive, sheet(all) modify
	putexcel A1="sample" B1="out" C1="or" D1="lci" E1="uci" F1="n" 

	
	regress `out' sd_grs_sbp_`threshold'_gwas2 eduyears age sex $PCs `covariate' if sample==1 & touse_sd_grs_sbp_`threshold'==1
	
	local out_label : var label `out'
	local lci = _b[sd_grs_sbp_`threshold'_gwas2]-1.96*_se[sd_grs_sbp_`threshold'_gwas2]
	local uci = _b[sd_grs_sbp_`threshold'_gwas2]+1.96*_se[sd_grs_sbp_`threshold'_gwas2]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[sd_grs_sbp_`threshold'_gwas2] D`x'=`lci' E`x'=`uci' 
	
	regress `out' sd_grs_sbp_`threshold'_gwas1 eduyears age sex $PCs `covariate' if sample==2 & touse_sd_grs_sbp_`threshold'==1

	
	local x = `x'+1

	local out_label : var label `out'
	local lci = _b[sd_grs_sbp_`threshold'_gwas1]-1.96*_se[sd_grs_sbp_`threshold'_gwas1]
	local uci = _b[sd_grs_sbp_`threshold'_gwas1]+1.96*_se[sd_grs_sbp_`threshold'_gwas1]
	
	putexcel A`x'="2" B`x'="`out_label'" C`x'=_b[sd_grs_sbp_`threshold'_gwas1] D`x'=`lci' E`x'=`uci'  
		
	
}	

levelsof eduyears, local(levels)
local x =3	
foreach out in sd_sbp  {


foreach l of local levels {

	
	putexcel set sbp_results_trait_`threshold'_additive, sheet(all) modify
	putexcel A1="sample" B1="out" C1="or" D1="lci" E1="uci" F1="n" G1="eduyears_sbp"

	
	regress `out' sd_grs_sbp_`threshold'_gwas2 age sex $PCs `covariate' if eduyears==`l' & sample==1 & touse_sd_grs_sbp_`threshold'==1
	
	local x=`x'+1
	local out_label : var label `out'
	local lci = _b[sd_grs_sbp_`threshold'_gwas2]-1.96*_se[sd_grs_sbp_`threshold'_gwas2]
	local uci = _b[sd_grs_sbp_`threshold'_gwas2]+1.96*_se[sd_grs_sbp_`threshold'_gwas2]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[sd_grs_sbp_`threshold'_gwas2] D`x'=`lci' E`x'=`uci' G`x'="Education `l'"
	
	regress `out' sd_grs_sbp_`threshold'_gwas1 age sex $PCs `covariate' if eduyears==`l' & sample==2 & touse_sd_grs_sbp_`threshold'==1

	local x=`x'+1
	
	local out_label : var label `out'
	local lci = _b[sd_grs_sbp_`threshold'_gwas1]-1.96*_se[sd_grs_sbp_`threshold'_gwas1]
	local uci = _b[sd_grs_sbp_`threshold'_gwas1]+1.96*_se[sd_grs_sbp_`threshold'_gwas1]
	
	putexcel A`x'="2" B`x'="`out_label'" C`x'=_b[sd_grs_sbp_`threshold'_gwas1] D`x'=`lci' E`x'=`uci'  G`x'="Education `l'"
		
	
}	
}


** Testing for interaction with education

foreach out in sd_sbp  {

	local x = 15
	
	putexcel set sbp_results_trait_`threshold'_additive, sheet(all) modify
	
	local x=`x'+1
	
	regress `out' sd_grs_sbp_`threshold'_gwas2 eduyears sd_grs_sbp_`threshold'_gwas2_int age sex $PCs `covariate' if touse_sd_grs_sbp_`threshold'==1 & sample==1
	
	local lci = _b[sd_grs_sbp_`threshold'_gwas2_int]-1.96*_se[sd_grs_sbp_`threshold'_gwas2_int]
	local uci = _b[sd_grs_sbp_`threshold'_gwas2_int]+1.96*_se[sd_grs_sbp_`threshold'_gwas2_int]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[sd_grs_sbp_`threshold'_gwas2_int] D`x'=`lci' E`x'=`uci' G`x'="EA interact"

	regress `out' sd_grs_sbp_`threshold'_gwas1 eduyears sd_grs_sbp_`threshold'_gwas1_int age sex $PCs `covariate' if touse_sd_grs_sbp_`threshold'==1 & sample==2
	
	
	local x=`x'+1
	local lci = _b[sd_grs_sbp_`threshold'_gwas1_int]-1.96*_se[sd_grs_sbp_`threshold'_gwas1_int]
	local uci = _b[sd_grs_sbp_`threshold'_gwas1_int]+1.96*_se[sd_grs_sbp_`threshold'_gwas1_int]
	
	putexcel A`x'="2" B`x'="`out_label'" C`x'=_b[sd_grs_sbp_`threshold'_gwas1_int] D`x'=`lci' E`x'=`uci' G`x'="EA interact"

	
}


	

*** Tab sample size 

foreach touse in  touse_sd_grs_smoking_`threshold' touse_sd_grs_sbp_`threshold'  {

		local sheet = substr("`touse'",7,.)
	
		putexcel set `date'_p`threshold'_trait_additive, sheet(`sheet') modify
		
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

*** Tab number of cases 

local pheno_out af cad diabetes_t2 stroke  
local primary_grs sd_grs_af_`threshold' sd_grs_cad_`threshold' sd_grs_diabetes_t2_`threshold' sd_grs_stroke_`threshold'  

local n : word count `pheno_out'

forvalues i = 1/`n' {

	local out : word `i' of `pheno_out'
	local exp : word `i' of `primary_grs'

	putexcel set `date'_p`threshold'_trait_additive, sheet(`exp') modify	
		
		
		tab `out' eduyears if touse_`exp'==1, matcell(numbers)
		 
		local n_7 = numbers[2,1]
		local n_10 = numbers[2,2]
		local n_13 = numbers[2,3]
		local n_15 = numbers[2,4]
		local n_19 = numbers[2,5]
		local n_20 = numbers[2,6]
		local n_all = `n_7'+`n_10'+`n_13'+`n_15'+`n_19'+`n_20'
		
	
		putexcel K2=`n_all' K3=`n_7' K4=`n_10' K5=`n_13' K6=`n_15' K7=`n_19' K8=`n_20'
	
	
	}
** Meta-analysing two samples - smoking

import excel "$resDir/data/results/genetic/202108/csi_results_trait_`threshold'_additive.xlsx", sheet("all") firstrow clear
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
	putexcel set `date'_p`threshold'_trait_additive, sheet(sd_grs_smoking_`threshold') modify
	putexcel A1="Exposure" B1="Outcome" C1="Interaction" D1="Scale" ///
		E1="Beta/OR" F1="LCI" G1="UCI" H1="P Value" I1 = "" ///
		J1 = "Sample Size" K1 = "Number of cases"
	local x = 2
	local beta = (r(ES))
	local lci = (r(ci_low))
	local uci = (r(ci_upp))
	local p_val = r(p_z)	
	putexcel A`x'="Smoking combined P=5x10-8 (SD)" B`x'="Smoking" C`x'="NONE" D`x'="additive" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_val'

levelsof eduyears_level, local(levels)
local x = 2
foreach l of local levels {

	putexcel set `date'_p`threshold'_trait_additive, sheet(sd_grs_smoking_`threshold') modify


	metan or lci uci if  eduyears_level==`l', lcols(sample ) effect(OR) null(1) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	
	local x = `x'+1
	local beta = (r(ES))
	local lci = (r(ci_low))
	local uci = (r(ci_upp))
	local p_val = r(p_z)	
	putexcel A`x'="Smoking combined P=5x10-8 (SD)" B`x'="Smoking" C`x'="Education `l'" D`x'="additive" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_val'
	
}

metan or lci uci if  eduyears_interact==1, lcols(sample ) effect(OR) null(1) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	putexcel set `date'_p`threshold'_trait_additive, sheet(sd_grs_smoking_`threshold') modify
	
	local x=9
	local beta_int = r(ES)
	local lci_int = r(ci_low)
	local uci_int = r(ci_upp)
	local p_value = r(p_z)

	
putexcel A`x'="Smoking combined P=5x10-8 (SD)" B`x'="SBP (SD)" C`x'="Interaction" D`x'="additive" E`x'=`beta_int' F`x'=`lci_int' G`x'=`uci_int' H`x'=`p_value' 



	** Meta-analysing two samples - SBP

import excel "$resDir/data/results/genetic/202108/sbp_results_trait_`threshold'_additive.xlsx", sheet("all") firstrow clear
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
	putexcel set `date'_p`threshold'_trait_additive, sheet(sd_grs_sbp_`threshold') modify
	putexcel A1="Exposure" B1="Outcome" C1="Interaction" D1="Scale" ///
		E1="Beta/OR" F1="LCI" G1="UCI" H1="P Value" I1 = "" ///
		J1 = "Sample Size" K1 = "Number of cases"
	local x = 2
	local beta = (r(ES))
	local lci = (r(ci_low))
	local uci = (r(ci_upp))
	local p_val = r(p_z)	
	putexcel A`x'="SBP combined P=5x10-8 (SD)" B`x'="SBP (SD)" C`x'="NONE" D`x'="additive" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_val'

levelsof eduyears_level, local(levels)
local x = 2
foreach l of local levels {

	putexcel set `date'_p`threshold'_trait_additive, sheet(sd_grs_sbp_`threshold') modify


	metan or lci uci if  eduyears_level==`l', lcols(sample ) effect(OR) null(1) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	
	local x = `x'+1
	local beta = (r(ES))
	local lci = (r(ci_low))
	local uci = (r(ci_upp))
	local p_val = r(p_z)	
	putexcel A`x'="SBP combined P=5x10-8 (SD)" B`x'="SBP (SD)" C`x'="Education `l'" D`x'="additive" E`x'=`beta' F`x'=`lci' G`x'=`uci' H`x'=`p_val'
	
}

metan or lci uci if  eduyears_interact==1, lcols(sample ) effect(OR) null(1) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	putexcel set `date'_p`threshold'_trait_additive, sheet(sd_grs_sbp_`threshold') modify

	local x=9
	local beta_int = r(ES)
	local lci_int = r(ci_low)
	local uci_int = r(ci_upp)
	local p_value = r(p_z)

	
putexcel A`x'="SBP combined P=5x10-8 (SD)" B`x'="SBP (SD)" C`x'="Interaction" D`x'="additive" E`x'=`beta_int' F`x'=`lci_int' G`x'=`uci_int' H`x'=`p_value' 


/*

Alice R Carter

18/06/2020

Code to create analysis data set for use in Educational attainment as a modifier of the effect of polygenic scores for cardiovascular risk factors: cross-sectional and prospective analysis of UK Biobank
*/

*do "$scriptDir/genetic/Code to select variables" 

cd "$resDir/data/stata/genetic"	

use "analysis_variables", clear

** Merge old linker file **
rename n_eid id_phe
merge 1:1 id_phe using "./polygenic_scores/old_linker.dta"
capture drop _merge 
rename id_phe n_eid

rename id_ieu_old appid8786

** Merge current linker file **
capture drop _merge
merge 1:1 n_eid using "linking_file"
keep if _m==3
capture drop _merge
rename v2 id_ieu
lab var id_ieu "Genetic Linker"

rename n_eid id_phe

** Withdrawals (and exclusions) **
merge 1:1 id_phe using "./exclusions/withdrawals.dta"
keep if _merge == 1
drop _merge
rename id_phe n_eid

** Genetic QC exlcuions **

*Remove recommended drops
foreach var in _recommended _highly_related _relateds _non_white_british {
    merge 1:1 id_ieu using "./exclusions/exclusions`var'.dta"
    keep if _merge == 1
    drop _merge
}

drop if n_3140_0_0 ==1
drop if n_3140_0_0 ==2

** Merge in sample 1 and sample 2 IDs for split sample GWAS analyses **

rename appid8786 IID
merge 1:1 IID using "sample_1_BP_GWAS.dta", keepusing(IID)
gen sample = 1 if _merge==3
drop if _merge==2
capture drop _merge

merge 1:1 IID using "sample_2_BP_GWAS.dta", keepusing(IID)
replace sample = 2 if _merge==3
drop if _merge==2
capture drop _merge

rename IID appid8786

** Phenotypic risk factor variables and covariates **

* Sex
rename n_31_0_0 sex	

* Age	
rename n_21003_0_0 age
lab var age "Age at baseline"

* BMI	
rename n_21001_0_0 bmi
lab var bmi "BMI"

* SBP
gen sbp = (n_4080_0_1+n_4080_0_0)/2
lab var sbp "Systolic blood pressure"

* DBP
gen dbp = (n_4079_0_1+n_4079_0_0)/2
lab var dbp "Diastolic blood pressure"
capture drop n_4079_*

* Total cholesterol
rename n_30690_0_0 total_chol
lab var total_chol "Total cholesterol"

* HDL-C (not currently in analysis plan)
rename n_30760_0_0 hdl_c
lab var hdl_c "HDL cholesterol"

* LDL-C (not currently in analysis sample data, add variable code)
rename n_30780_0_0 ldl_c
lab var ldl_c "LDL cholesterol"

* Smoking - may want to use a different phenotypic measure
gen smoke_cat = 0 if n_20116_0_0==0 
replace smoke_cat = 1 if n_20116_0_0==1 & smoke_cat==.
replace smoke_cat = 2 if n_3456_0_0 >=-10 & n_3456_0_0 <=9 & smoke_cat==.
replace smoke_cat = 3 if n_3456_0_0 >=10 & n_3456_0_0 <=19 &smoke_cat==.
replace smoke_cat = 4 if n_3456_0_0 >=20 & smoke_cat==.
lab var smoke_cat "Smoking status"
lab def smoke_cat 0 "Never" 1 "Former" 2 "Light 1-9/day" 3"Moderate 10-19/day" 4 "Heavy >20/day"
lab val smoke_cat smoke_cat

merge 1:1 appid8786 using "full_sample_csi", keepusing (csi)
drop if _merge==2
capture drop _merge

* drinks per week - using code from Laurence Howe from doi: 10.1038/s41467-019-12424-x
foreach var in n_1568_0_0 n_1578_0_0 n_1588_0_0 n_1598_0_0 n_1608_0_0 {
	
	replace `var' = 0 if `var'==-1
	replace `var' = . if `var'==-3
}

*Summing units as estimates of unit per each drink type
gen units_per_week = (2*n_1568_0_0)+(2*n_1578_0_0)+(2.5*n_1588_0_0)+(n_1598_0_0)+(2*n_1608_0_0)

summ units_per_week
gen units_fivesd = (abs((units_per_week-r(mean))/r(sd))>5) if units_per_week<.
replace units_per_week=. if units_fivesd==1

*replacing those who drink 1-3 times per week, special occasions or never to 0 units
replace units_per_week = 0 if n_1558_0_0==4 | n_1558_0_0==5 |n_1558_0_0==6

* Summing drinks per week not accounting for units consumed
gen drinks_per_week = (n_1568_0_0)+(n_1578_0_0)+(n_1588_0_0)+(n_1598_0_0)+(n_1608_0_0)
replace drinks_per_week = 0 if n_1558_0_0==4 | n_1558_0_0==5 |n_1558_0_0==6

summ drinks_per_week
gen drinks_fivesd = (abs((drinks_per_week-r(mean))/r(sd))>5) if drinks_per_week<.
replace drinks_per_week=. if drinks_fivesd==1

* Family history of CVD - not included in analysis plan but may be useful
gen fh_cvd=.
foreach var of varlist n_20111_0_1-n_20111_0_11{
replace fh_cvd =1 if `var'==1
}

foreach var of varlist n_20107_0_1-n_20107_0_9{
replace fh_cvd =1 if `var'==1
}

foreach var of varlist n_20110_0_1-n_20110_0_10{
replace fh_cvd =1 if `var'==1
}

replace fh_cvd=0 if fh_cvd==.
lab var fh_cvd "CVD in first degree relative"
lab def fh_cvd 0 "No CVD" 1 "CVD"
lab val fh_cvd fh_cvd


** Socioeconomic variables ** 
gen eduyears=.
ds n_6138_*
foreach i in `r(varlist)'{
	replace eduyears=20 if `i'==1
	replace eduyears=19 if `i'==5 & (eduyears<19|eduyears==.) 
	replace eduyears=15 if `i'==6 & (eduyears<15|eduyears==.)
	replace eduyears=13 if `i'==2 & (eduyears<13|eduyears==.)
	replace eduyears=10 if (`i'==3|`i'==4) & (eduyears<10|eduyears==.)
	replace eduyears=7 if `i'==-7 & (eduyears<7|eduyears==.)
	}
lab var eduyears "Years of education"

gen edu_low = 0 if eduyears ==7 | eduyears ==10
replace edu_low = 1 if eduyears >=13 & edu_low==.

rename n_738_0_0 income
lab var income "Household income (pre tax)"

rename n_680_0_0 home_owner
lab var home_owner "Own or rent home"

rename n_4674_0_0 private_health
lab var private_health "Use of private healthcare"

rename n_189_0_0 tdi
xtile tdi_cat=tdi, n(5)
lab var tdi_cat "Quintiles of TDI"	
rename tdi town	


********************************************************************************
*							CVD Cases Data									   *
********************************************************************************
*Adding HES dates
drop if n_eid==.
merge 1:1 n_eid using "hes_exclusions_and_cases_202007.dta", keepusing(CVD_to_exclude *_inc *_prev)
drop if _m==2
drop _merge

local var


foreach var of varlist Stroke_inc AMI_inc t1_inc ckd_inc angina_inc TIA_inc PAD_inc FH_inc CVD_inc IHD_inc mdd_inc af_inc arthritis_inc t2_inc lupus_inc impotence_inc migraine_inc hiv_inc {
	local var2 =substr("`var'",1, length("`var'")-4)
	replace `var' = 0 if `var'==. &  `var2'_prev!=1
	}
	

foreach var of varlist Stroke_prev AMI_prev t1_prev ckd_prev angina_prev TIA_prev PAD_prev FH_prev CVD_prev IHD_prev mdd_prev af_prev arthritis_prev t2_prev lupus_prev impotence_prev migraine_prev hiv_prev {
	replace `var'=0 if `var'==.
}	

label define casecontrol 0 "control" 1 "case", modify


lab var AMI_inc "Incident AMI Case"
lab var AMI_prev "Prevalent AMI Case"
*lab var AMI_all "All AMI cases"

lab var Stroke_inc "Incident Stroke Case"
lab var Stroke_prev "Prevalent Stroke Case"
*lab var Stroke_all "All Stroke cases"

lab var IHD_inc "Incident IHD Case"
lab var IHD_prev "Prevalent IHD Case"
*lab var IHD_all "All IHD cases"

lab var CVD_inc "Incident CVD Case"
lab var CVD_prev "Prevalent CVD Case"
*lab var CVD_all "All CVD cases"

lab val *_inc casecontrol
lab val *_prev casecontrol
*lab val *_all casecontrol

replace CVD_to_exclude=0 if CVD_to_exclude==.
lab val CVD_to_exclude casecontrol
drop if n_eid==.
capture drop _merge

*Adding diabetes cases
replace t2_prev = 1 if n_2443_0_0==1
replace t2_inc=1 if n_2443_1_0==1 & t2_prev!=1
replace t2_inc=1 if n_2443_2_0==1 & t2_prev!=1
						
* Merge in medication data from phenotypic "1B - Medications script.do" for use in QRISK scores

merge 1:1 n_eid using "V:/data/stata/phenotypic/medication_data", ///
	keepusing(statin medication_female medication_male statin_selfreport statin_cat medications b_treatedhyp b_atypicalantipsy erectile_dysfunction b_corticosteroids)
drop if _m==2
capture drop _merge

* Adjust BP measurements for antihypertensive use - not used in analyses
gen sbp_plus10 = sbp
replace sbp_plus10 = sbp_plus10+10 if b_treatedhyp==1
lab var sbp_plus10 "SBP Adjusted for antihypertensive use"

gen dbp_plus5 = dbp
replace dbp_plus5 =dbp_plus5+5 if b_treatedhyp==1
lab var dbp_plus5 "DBP adjusted for antihypertensive use - 5mmHg"


replace impotence_prev=1 if erectile_dysfunction==1

* Merging in polygenic scores

merge 1:1 id_ieu using "./polygenic_scores/all_scores/stata_files/pgs_all.dta"
keep if _merge==3
capture drop _merge

merge 1:1 id_ieu using "./polygenic_scores/education_score.dta"
keep if _merge==3
capture drop _merge

rename n_22009_0_* PC*

drop if eduyears==.
drop if CVD_to_exclude==1

* Keep required variables
keep n_eid sex n_54_0_0 town home_owner income private_health bmi age PC* total_chol hdl_c ldl_c id_ieu sample sbp dbp smoke_cat csi units_per_week units_fivesd drinks_per_week drinks_fivesd fh_cvd eduyears edu_low tdi_cat Stroke_prev AMI_prev t1_prev ckd_prev angina_prev TIA_prev PAD_prev FH_prev CVD_prev IHD_prev mdd_prev af_prev t2_prev arthritis_prev lupus_prev impotence_prev migraine_prev hiv_prev Stroke_inc AMI_inc t1_inc ckd_inc angina_inc TIA_inc PAD_inc FH_inc CVD_inc IHD_inc mdd_inc af_inc t2_inc arthritis_inc lupus_inc impotence_inc migraine_inc hiv_inc CVD_to_exclude statin medication_female medication_male statin_selfreport statin_cat b_treatedhyp b_atypicalantipsy erectile_dysfunction b_corticosteroids medications sbp_plus10 dbp_plus5 grs_af_05 grs_af_5 grs_af_5e8 grs_alcohol_05 grs_alcohol_5 grs_alcohol_5e8 grs_bmi_05 grs_bmi_5 grs_bmi_5e8 grs_cad_05 grs_cad_5 grs_cad_5e8 grs_diabetes_t2_05 grs_diabetes_t2_5 grs_diabetes_t2_5e8 grs_ldl_05 grs_ldl_5 grs_ldl_5e8  grs_sbp_05_gwas1 grs_sbp_05_gwas2  grs_sbp_5e8_gwas1 grs_sbp_5e8_gwas2 grs_sbp_5_gwas1 grs_sbp_5_gwas2 grs_sle_5e8  grs_smoking_05_gwas1 grs_smoking_05_gwas2 grs_smoking_5e8_gwas1 grs_smoking_5e8_gwas2 grs_smoking_5_gwas1 grs_smoking_5_gwas2 grs_stroke_05 grs_stroke_5 grs_stroke_5e8 grs_total_chol_05 grs_total_chol_5 grs_total_chol_5e8 ea_weighted

* Save master data
save  "GxE_master_202106.dta", replace

********************************************************************************
				* Setting up analysis variables *
********************************************************************************
			
mark touse
markout touse eduyears grs_*
keep if touse==1

global all_grs grs_* ea_weighted
global PCs PC1-PC40

* Create SDs for polygenic scores and label
foreach var of varlist $all_grs {
	
	egen sd_`var' = std(`var')
}

lab var sd_ea_weighted "Educational attainment PGS (SD)"

lab var sd_grs_alcohol_5e8 "Drinks per week P=5x10-8 (SD)"
lab var sd_grs_af_5e8 "Atrial fibrillation P=5x10-8 (SD)"
lab var sd_grs_bmi_5e8 "BMI P=5x10-8 (SD)"
lab var sd_grs_cad_5e8 "CAD P=5x10-8 (SD)"
lab var sd_grs_diabetes_t2_5e8 "Diabetes (T2) P=5x10-8 (SD)"
lab var sd_grs_ldl_5e8 "LDL-C P=5x10-8 (SD)"
lab var sd_grs_stroke_5e8 "Stroke P=5x10-8 (SD)"
lab var sd_grs_total_chol_5e8 "Total Cholesterol P=5x10-8 (SD)"

lab var sd_grs_alcohol_05 "Drinks per week P=0.05 (SD)"
lab var sd_grs_af_05 "Atrial fibrillation P=0.05 (SD)"
lab var sd_grs_bmi_05 "BMI P=0.05 (SD)"
lab var sd_grs_cad_05 "CAD P=0.05 (SD)"
lab var sd_grs_diabetes_t2_05 "Diabetes (T2) P=0.05 (SD)"
lab var sd_grs_ldl_05 "LDL-C P=0.05 (SD)"
lab var sd_grs_stroke_05 "Stroke P=0.05 (SD)"
lab var sd_grs_total_chol_05 "Total Cholesterol P=0.05 (SD)"

lab var sd_grs_alcohol_5 "Drinks per week P=0.5 (SD)"
lab var sd_grs_af_5 "Atrial fibrillation P=0.5 (SD)"
lab var sd_grs_bmi_5 "BMI P=0.5 (SD)"
lab var sd_grs_cad_5 "CAD P=0.5 (SD)"
lab var sd_grs_diabetes_t2_5 "Diabetes (T2) P=0.5 (SD)"
lab var sd_grs_ldl_5 "LDL-C P=0.5 (SD)"
lab var sd_grs_stroke_5 "Stroke P=0.5 (SD)"
lab var sd_grs_total_chol_5 "Total Cholesterol P=0.5 (SD)"


* Create interaction parameter for polygenic scores and education (observational [obs] and genetic [gen])
global sd_cv_grs  sd_grs_af_05  sd_grs_af_5  sd_grs_af_5e8  sd_grs_alcohol_05  sd_grs_alcohol_5  sd_grs_alcohol_5e8  sd_grs_bmi_05  sd_grs_bmi_5  sd_grs_bmi_5e8  sd_grs_cad_05  sd_grs_cad_5  sd_grs_cad_5e8  sd_grs_diabetes_t2_05  sd_grs_diabetes_t2_5  sd_grs_diabetes_t2_5e8  sd_grs_ldl_05  sd_grs_ldl_5  sd_grs_ldl_5e8 sd_grs_sle_5e8   sd_grs_stroke_05  sd_grs_stroke_5  sd_grs_stroke_5e8  sd_grs_total_chol_05  sd_grs_total_chol_5  sd_grs_total_chol_5e8 

foreach var of varlist $sd_cv_grs {
	
	gen `var'_int = `var'*eduyears
	
}

global sample_1 sd_grs_smoking_05_gwas1 sd_grs_smoking_5e8_gwas1 sd_grs_smoking_5_gwas1 sd_grs_sbp_05_gwas1 sd_grs_sbp_5e8_gwas1 sd_grs_sbp_5_gwas1

global sample_2 sd_grs_smoking_05_gwas2 sd_grs_smoking_5e8_gwas2    sd_grs_smoking_5_gwas2 sd_grs_sbp_05_gwas2 sd_grs_sbp_5e8_gwas2 sd_grs_sbp_5_gwas2

foreach var of varlist $sample_1 {
	
	gen `var'_int = `var'*eduyears if sample==2
}

foreach var of varlist $sample_2 {
	
	gen `var'_int = `var'*eduyears if sample==1
}

* Rename/label phenotypic variables to match exactly the polygenic score label
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

rename ea_weighted grs_eduyears
rename sd_ea_weighted sd_grs_eduyears

* Create SDs for observed continuous traits
global cont_trait alcohol bmi ldl total_chol smoking sbp

foreach var of varlist $cont_trait {

egen sd_`var' = std(`var')

}

* Create interaction terms for confounders/covariates
* Observational education
foreach confounder of varlist age sex $PCs {
	
	gen eduyears_obs_`confounder'_X = eduyears*`confounder'
	
	}

* Genetic education	
foreach confounder of varlist age sex $PCs {
	
	gen eduyears_gen_`confounder'_X = sd_grs_eduyears*`confounder'
	
	}

* CV polygenic scores	
foreach var of varlist $sd_cv_grs {
	
	foreach confounder of varlist age sex $PCs {
	
	gen `var'_`confounder'_X = `var'*`confounder'
	
	}
}

* Save analysis data
save  "GxE_analysis_202108.dta", replace
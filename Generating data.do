
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

merge 1:1 n_eid using "medication_data", ///
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

** Merging Statin primary care data **
merge 1:1 n_eid using  "baseline_receiver", keepusing (baseline_receiver_* *_event_statin statin_incident statin_end_prescribing drug_name)
gen statin_primary = 1 if _merge==3
drop if _merge==2
capture drop _merge

merge 1:1 n_eid using "analysis_variables_inc_primary_care"
drop if _merge==2
capture drop _merge

gen primary_care_presc = 1 if  n_42039_0_0!=.
replace primary_care_presc=0 if  n_42039_0_0==. & primary_care_presc==.
lab var primary_care_presc "Data in primary care prescription records"
lab def primary_care 0 "Not in primary care" 1 "Data in primary care", modify
lab val primary_care_presc primary_care

replace statin_primary=0 if statin_primary==.
lab var statin_primary "Statin prescribed in primary care data"
lab def statin_primary 0 "No prescription" 1 "Statin prescription"
lab val statin_primary statin_primary

replace statin_incident=0 if statin_incident==. & first_event_statin==. & statin==0
lab var statin_incident "Incident statin prescription from primanry care"

*replace baseline_receiver_3month=0 if baseline_receiver==. 
lab var baseline_receiver_3month "Prescription +/- 3 months of baseline"

*replace baseline_receiver_6month=0 if baseline_receiver_6month==.
lab var baseline_receiver_6month "Prescription +/- 6 months of baseline"

*replace baseline_receiver_9month=0 if baseline_receiver_9month==.
lab var baseline_receiver_9month "Prescription +/- 9 months of baseline"

*replace baseline_receiver_12month=0 if baseline_receiver_12month==.
lab var baseline_receiver_12month "Prescription +/- 12 months of baseline"

gen non_complier_3month = 1 if baseline_receiver_3month==1 & statin==0
replace non_complier_3month = 0 if baseline_receiver_3month==1 & statin==1
replace non_complier_3month = -9 if primary_care_presc==0
lab var non_complier_3month "Prescription +/- 3 months of baseline, but not reported at baselines "
lab def non_complier 0 "Complier" 1 "Non complier" -9 "Not in primary care data", modify
lab val non_complier_3month non_complier

lab var first_event_statin "Date of first statin prescription"
lab var last_event_statin "Date of last statin prescription"
lab var drug_name "Statin name"
lab var statin_incident "Statin prescription from primary care - only after baseline"
lab var statin_end_prescribing "No statin prescription within 12 months prior to baseline"

*capture drop duplicate id  issue_date read_2 bnf_code quantity baseline_receiver dmd_code data_provider

* Generate variable for over the counter statin use *
gen statin_otc = 1 if statin==1 & statin_primary==0 & primary_care==1
replace statin_otc = 1 if statin_primary==1 & statin==1 & first_event_statin>ts_53_0_0 & primary_care==1 & statin_otc==.
replace statin_otc=0 if statin==0 & statin_primary==0 & primary_care==1 /* Coded as 0 if no statin at all */


* Merging in polygenic scores

merge 1:1 id_ieu using "/pgs_all.dta"
keep if _merge==3
capture drop _merge

rename n_22009_0_* PC_*

drop if eduyears==.
drop if CVD_to_exclude==1




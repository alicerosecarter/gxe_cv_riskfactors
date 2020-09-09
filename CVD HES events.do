
rename eid n_eid
*** ICD-10 codes

generate Stroke_diag = 0
lab var Stroke_diag "Diagnosis - Stroke"
generate AMI_diag = 0
lab var AMI_diag "Diagnosis - Acute MI"
generate IHD_diag = 0
lab var IHD_diag "Diagnosis - IHD"
generate CVD_diag = 0
lab var CVD_diag "Diagnosis - CVD"
generate t1_diag = 0
lab var t1_diag "Diagnosis - Diabetes Type I"
generate ckd_diag = 0
lab var ckd_diag "Diagnosis - Chronic Kidney Disease"
generate angina_diag = 0
lab var angina_diag "Diagnosis - Angina"
generate TIA_diag = 0
lab var TIA_diag "Diagnosis - Transient Ischaemic Disease"
generate PAD_diag = 0
lab var PAD_diag "Diagnosis - Peripheral Arterial Disease"
generate FH_diag = 0
lab var FH_diag "Diagnosis - Familial Hypercholesterolaemia"
gen mdd_diag = 0
lab var mdd_diag "Diagnosis - Severe mental illness"
generate af_diag = 0
lab var af_diag "Diagnosis - Atrial Fibrillation"
generate arthritis_diag = 0
lab var arthritis_diag "Diagnosis - Arthritis"
generate t2_diag = 0
lab var t2_diag "Diagnosis - Diabetes Type II"
generate lupus_diag = 0
lab var lupus_diag "Diagnosis - Lupus"
generate migraine_diag = 0
lab var migraine_diag "Diagnosis - Migraine"
generate hiv_diag = 0
lab var hiv_diag "Diagnosis - HIV/AIDs"
generate ed_diag = 0
lab var ed_diag "Diagnosis - Erectile Dysfunction"


label define casecontrol 0 "Control" 1 "Case", modify
label values *_diag casecontrol

global diagnoses Stroke_diag AMI_diag t1_diag ckd_diag angina_diag TIA_diag PAD_diag FH_diag CVD_diag IHD_diag mdd_diag af_diag arthritis_diag t2_diag lupus_diag ed_diag migraine_diag hiv_diag 

foreach cause of varlist $diagnoses {
    generate IHD_code_`cause' =.
	generate Stroke_code_`cause' =.
	generate CVD_code_`cause' =.
	generate AMI_code_`cause' =.
	generate t1_code_`cause' =.
	generate ckd_code_`cause' =.
	generate angina_code_`cause' =.
	generate tia_code_`cause' =.
	generate pad_code_`cause' =.
	generate fh_code_`cause' =.
	generate mdd_code_`cause' =.
	generate af_code_`cause' =.
	generate arthritis_code_`cause' =.
	generate t2_code_`cause' =.
	generate lupus_code_`cause' =.
	generate migraine_code_`cause' =.
	generate hiv_code_`cause' =.
	generate ed_code_`cause'=.
	}
	
* ICD 10 codes
foreach diag of varlist diag_icd10_* {
foreach cause of varlist $diagnoses {
	* generate variable to capture ICD codes relating to IHD (ICDI20-I25)
    replace IHD_code_`cause' = strpos(`diag', "I20") > 0 
	replace IHD_code_`cause' = strpos(`diag', "I21") > 0 if IHD_code_`cause' == 0
	replace IHD_code_`cause' = strpos(`diag', "I22") > 0 if IHD_code_`cause' == 0
	replace IHD_code_`cause' = strpos(`diag', "I23") > 0 if IHD_code_`cause' == 0
	replace IHD_code_`cause' = strpos(`diag', "I24") > 0 if IHD_code_`cause' == 0
	replace IHD_code_`cause' = strpos(`diag', "I25") > 0 if IHD_code_`cause' == 0

	* generate variable to capture ICD codes relating to stroke (ICDI60-69, ICDG45*)
	replace Stroke_code_`cause' = strpos(`diag', "I6") > 0 
	replace Stroke_code_`cause' = strpos(`diag', "G45") > 0 if Stroke_code_`cause' == 0
	
	* generate variable to capture ICD codes relating to CVD (ICDI*)
	replace CVD_code_`cause' = strpos(`diag', "I") > 0 
	replace CVD_code_`cause' = strpos(`diag', "G45") > 0 if CVD_code_`cause' == 0

	* generate variable to capture ICD codes relating to AMI (ICDI21*)	
	replace AMI_code_`cause' = strpos(`diag', "I21") > 0 
	replace AMI_code_`cause' = strpos(`diag', "I22") > 0 if AMI_code_`cause' == 0

	* generate variable to capture ICD codes relating to T1 diabetes (E10)
	replace t1_code_`cause' = strpos(`diag', "E10") > 0

	* generate variable to capture ICD codes relating to Chronic kidney disease stages 3-5 (N183-N185)
	replace ckd_code_`cause' = strpos(`diag', "N183") > 0
	replace ckd_code_`cause' = strpos(`diag', "N184") > 0 if ckd_code_`cause' == 0
	replace ckd_code_`cause' = strpos(`diag', "N185") > 0 if ckd_code_`cause' == 0
	
	* generate variable to capture ICD codes relating to Angina (I20.0-I20.9)
	replace angina_code_`cause' = strpos(`diag', "I20") >0
	
	* generate variable to capture ICD codes relating to Transient Arterial Attack (G45)
	replace tia_code_`cause' = strpos(`diag', "G45") >0
	
	* generate variable to capture ICD codes relating to Peripheral arterial disease (I73.9)
	replace pad_code_`cause' = strpos(`diag', "I73.9") > 0 
	
	* generate variable to capture ICD codes relating to familial hypercholesterolaemia (I78)
	replace fh_code_`cause' = strpos(`diag', "I78.0") > 0
	
	* generate variable to capture ICD codes relating to severe mental illness/major depressive disorder
	replace mdd_code_`cause' = strpos(`diag', "F20") > 0
	replace mdd_code_`cause' = strpos(`diag', "F23") > 0 if mdd_code_`cause' ==0
	replace mdd_code_`cause' = strpos(`diag', "F31") > 0 if mdd_code_`cause' ==0
	replace mdd_code_`cause' = strpos(`diag', "F32") > 0 if mdd_code_`cause' ==0
	replace mdd_code_`cause' = strpos(`diag', "F33") > 0 if mdd_code_`cause' ==0

	* generate variable to capture ICD codes relating to atrial fibrillation
	replace af_code_`cause' = strpos(`diag', "I48") > 0
	
	* generate variable to capture ICD codes relating to arthritis
	replace arthritis_code_`cause' = strpos(`diag', "M05") > 0
	
	* generate variable to capture ICD codes relating to type2 diabetes
	replace t2_code_`cause' = strpos(`diag', "E11") > 0
	
	* generate variable to capture ICD codes relating to lupus
	replace lupus_code_`cause' = strpos(`diag', "M329") > 0
	
	* generate variable to capture ICD codes relating to migraine	
	replace migraine_code_`cause' = strpos(`diag', "G43") > 0

	* generate variable to capture ICD codes relating to HIV
	replace hiv_code_`cause' = strpos(`diag', "B20") > 0

	* generate variable to capture ICD codes relating to erectile dysfunction/impotence
	replace ed_code_`cause' = strpos(`diag', "N52") > 0
	
	* update case status
	replace Stroke_diag = 1 if Stroke_code_`cause' == 1
	replace AMI_diag = 1 if AMI_code_`cause' == 1
	replace IHD_diag = 1 if IHD_code_`cause' == 1
	replace CVD_diag = 1 if CVD_code_`cause' == 1
	replace t1_diag 	= 1 if t1_code_`cause' ==1
	replace ckd_diag 	 = 1 if ckd_code_`cause' ==1
	replace angina_diag = 1 if angina_code_`cause' ==1
	replace TIA_diag = 1 if tia_code_`cause' ==1
	replace PAD_diag = 1 if pad_code_`cause' ==1
	replace FH_diag = 1 if fh_code_`cause' ==1
	replace mdd_diag = 1 if mdd_code_`cause'==1
	replace af_diag = 1 if af_code_`cause'==1
	replace arthritis_diag 	= 1 if arthritis_code_`cause' ==1
	replace t2_diag 	= 1 if t2_code_`cause' ==1
	replace lupus_diag 	= 1 if lupus_code_`cause' ==1
	replace migraine_diag = 1 if migraine_code_`cause' ==1
	replace hiv_diag 	 = 1 if hiv_code_`cause' ==1
	replace ed_diag 	 = 1 if ed_code_`cause' ==1
	
	}
}
drop *_code_*


*ICD 9 codes

* generate variable to capture ICD codes relating to IHD (ICD 410-414 & 4100-4149)
gen byte IHD_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 410/414 {
        replace  IHD_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 
								
								
* generate variable to capture ICD codes relating to stroke (ICDI60-69, ICDG45*)
gen byte Stroke_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 430/438 {
        replace Stroke_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 
		               
* generate variable to capture ICD codes relating to CVD (ICDI*)
gen byte CVD_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 390/459 {
        replace CVD_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 
		
* generate variable to capture ICD codes relating to AMI (ICDI21*)             
gen byte AMI_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 410/410 {
        replace AMI_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
	}
foreach var of varlist diag_icd9* {	
	forvalues i = 412/412 {
        replace AMI_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
} 
}
										
* generate variable to capture ICD codes relating to Type 1 diabetes (ICD 250*)
gen byte t1d_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 2500/2509 {
        replace t1d_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 	
foreach var of varlist diag_icd9* {	
    forvalues i = v5867/v5867 {
replace t1d_code_ICD9 = 1 if substr(`var',1,4) == "`i'"
}
}					
* generate variable to capture ICD codes relating to CKD stage 3-5 (ICD 585.3-585.5)             
gen byte ckd_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 5853/5855 {
        replace ckd_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 		
													
* generate variable to capture ICD codes relating to Angina
gen byte angina_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 4139/4139 {
        replace angina_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 

* generate variable to capture ICD codes relating to TIA
gen byte TIA_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 4359/4359 {
        replace TIA_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 
								
* generate variable to capture ICD codes relating to PAD
gen byte PAD_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 4439/4439 {
        replace PAD_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 

* generate variable to capture ICD codes relating to FH							
gen byte FH_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 2720/2720 {
        replace FH_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
}

* generate variable to capture ICD codes relating to severe mental diagnosis (ICD 295-296)             
gen byte mdd_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 2953/2953   {
        replace mdd_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 

foreach var of varlist diag_icd9* {
    forvalues i = 29581/29581   {
        replace mdd_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 

foreach var of varlist diag_icd9* {
    forvalues i = 2964/2965   {
        replace mdd_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 

foreach var of varlist diag_icd9* {
    forvalues i = 2967/2967   {
        replace mdd_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 

foreach var of varlist diag_icd9* {
    forvalues i = 2958/2958    {
        replace mdd_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 

foreach var of varlist diag_icd9* {
    forvalues i = 2962/2963   {
        replace mdd_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 
* generate variable to capture ICD codes relating to Atrial fibrilation (ICD 427.31)             
gen byte AF_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 42731/42731 {
        replace AF_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
}
 
* generate variable to capture ICD codes relating to Rheumatois arthritis (ICD 714)
gen byte arthritis_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 714/714 {
        replace arthritis_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 

* generate variable to capture ICD codes relating to Type 2 diabetes (ICD 250* )
gen byte t2d_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 2500/2509 {
        replace t2d_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 
* generate variable to capture ICD codes relating to lupus (ICD 710)             
gen byte lupus_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 710/710 {
        replace lupus_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 
								
* generate variable to capture ICD codes relating to migraine (ICD 346)             
gen byte migraine_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 346/346 {
        replace migraine_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 
								
* generate variable to capture ICD codes relating to HIV (ICD 042, v08, 795.71 **TO ADD  V08)             
gen byte hiv_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 042/042 {
        replace hiv_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
	}	
    }
foreach var of varlist diag_icd9* {

	forvalues i = 79571/79571 {
        replace hiv_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
	}
	}
							
* generate variable to capture ICD codes relating to erectile dysfunction (ICD 607.84)             
gen byte ed_code_ICD9 = 0
foreach var of varlist diag_icd9* {
    forvalues i = 60784/60784 {
        replace ed_code_ICD9 = 1 if substr(`var',1,3) == "`i'"
    }
} 
													
replace CVD_diag = 1 if CVD_code_ICD9 ==1
replace Stroke_diag = 1 if Stroke_code_ICD9 ==1
replace IHD_diag = 1 if IHD_code_ICD9 ==1
replace AMI_diag = 1 if AMI_code_ICD9 ==1
replace t1_diag 	= 1 if t1d_code_ICD9 ==1
replace ckd_diag 	 = 1 if ckd_code_ICD9 ==1
replace angina_diag = 1 if angina_code_ICD9==1
replace TIA_diag = 1 if TIA_code_ICD9 ==1
replace PAD_diag = 1 if PAD_code_ICD9 ==1
replace FH_diag = 1 if FH_code_ICD9==1
replace mdd_diag = 1 if mdd_code_ICD9==1
replace af_diag = 1 if AF_code_ICD9==1
replace arthritis_diag 	= 1 if arthritis_code_ICD9 ==1
replace t2_diag 	= 1 if t2d_code_ICD9 ==1
replace lupus_diag 	= 1 if lupus_code_ICD9 ==1
replace migraine_diag = 1 if migraine_code_ICD9 ==1
replace hiv_diag 	 = 1 if hiv_code_ICD9 ==1
replace ed_diag 	 = 1 if ed_code_ICD9 ==1

merge m:1  n_eid using "analysis_variables.dta", keepusing(ts_53_0_0) nogen keep(3)

gen date = epistart
replace date = epiend if date==.
replace date = date(admidate, "YMD") if date==.
format date %td

global diagnoses Stroke_diag AMI_diag t1_diag ckd_diag angina_diag TIA_diag PAD_diag FH_diag CVD_diag IHD_diag mdd_diag af_diag t2_diag arthritis_diag lupus_diag ed_diag migraine_diag hiv_diag 
local x=0
foreach var of varlist $diagnoses { 

	local x = `x'+1
	gen v`x'_prev = `var' if date <=ts_53_0_0
	lab var v`x'_prev `var'
	replace v`x'_prev=0 if v`x'_prev==.

}

local x=0
foreach var of varlist $diagnoses { 

	local x = `x'+1
	gen v`x'_inc = `var' if date >ts_53_0_0
	lab var v`x'_inc `var'

}

local icd_total=18

forvalues i = 1/`icd_total' {
	bysort n_eid: egen x = max(v`i'_prev)
	replace v`i'_prev = x 

	drop x
}

local icd_total=18

forvalues i = 1/`icd_total' {
	bysort n_eid: egen x = max(v`i'_inc)
	replace v`i'_inc = x 
	drop x
}


forvalues i = 1/`icd_total' {

	replace v`i'_inc = . if  v`i'_prev==1

}

duplicates tag n_eid, gen(duplicate)

duplicates drop n_eid, force
egen count = rowtotal(v*)
drop if count == 0 
drop count


rename v1_* Stroke_*
rename v2_* AMI_*
rename v3_* t1_*
rename v4_* ckd_*
rename v5_* angina_*
rename v6_* TIA_*
rename v7_* PAD_*
rename v8_* FH_*
rename v9_* CVD_*
rename v10_* IHD_*
rename v11_* mdd_*
rename v12_* af_*
rename v13_* t2_*
rename v14_* arthritis_*
rename v15_* lupus_*
rename v16_* impotence_*
rename v17_* migraine_*
rename v18_* hiv_*


gen CVD_to_exclude = 1 if Stroke_prev==1 | AMI_prev==1 | t1_prev==1 | ckd_prev==1 | angina_prev==1 | TIA_prev==1 | PAD_prev==1 | FH_prev==1
lab var CVD_to_exclude "Individuals to exclude based on previous diagnoses"

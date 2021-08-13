/*

Alice R Carter

18/06/2021

Main do file for specifying analyses
*/

********************************************************
* Additive interactions
********************************************************

* Genome-wide significant (P<5x10-8)
do "$scriptDir/3_analysis_subdo_additive.do" "5e8" "202108" 

* P<0.05
do "$scriptDir/3_analysis_subdo_additive.do" "05" "202108" 

* P<0.5
do "$scriptDir/3_analysis_subdo_additive.do" "5" "202108" 

********************************************************
* Multiplicative interactions
********************************************************

* Genome-wide significant (P<5x10-8)
do "$scriptDir/4_analysis_subdo_multiplicative.do" "5e8" "202108" 

* P<0.05
do "$scriptDir/4_analysis_subdo_multiplicative.do" "05" "202108" 

* P<0.5
do "$scriptDir/4_analysis_subdo_multiplicative.do" "5" "202108" 

********************************************************
* Main paper figures
********************************************************

* Additive scale
do "$scriptDir/5_figures_additive.do" "5e8" "202108" 

* Multiplicative scale
do "$scriptDir/6_figures_multiplicative.do" "5e8" "202108" 

* Interaction coefficients 
do "$scriptDir/7_figures_interaction_coef.do" "5e8" "202108" 

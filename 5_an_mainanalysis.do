/*=========================================================================
DO FILE NAME:		5_an_mainanalysis

AUTHOR:				Jennifer Davidson (edited do file provided by Kate 
					Mansfield)		
DATE:				13/05/2021
			
PURPOSE:			Regression results from SCCS analysis
			
			
DATASET USED: 		analysis_macelabconfirmsarscov2 & 
					analysis_maceclinreportcovid19
													
NEXT STEP:			6_an_cv+singlestratified.do	
					
*=========================================================================*/

clear all
capture log close
log using "$logdir\mainregression.log", replace

/*******************************************************************************
#1. Loop through datasets
*******************************************************************************/

local datasets labconfirmsarscov2 clinreportcovid19
local outcome mace acs mi /*angina*/ stroke hf arrhythmia

foreach suffix of local datasets {

	use "$datadir/analysis_mace`suffix'", clear
	merge m:1 newid using "$datadir/eventtype_`suffix'", keep(master match) nogen
	
	unique newid
	gen mace=1
/*******************************************************************************
#2. Create csv for output
*******************************************************************************/
foreach cond of local outcome {
	
	cap file close textfile 
	file open textfile using "$outputdir/`cond'`suffix'.csv", write replace
	file write textfile "sep=;" _n
	file write textfile "Risk period" ";" "All" ";" ";" "QRISK3" ";" ";" ";" ";" "Hypertension" _n
	file write textfile ";" ";" ";" "Raised risk" ";"  ";" "Low risk" ";" ";" "Raised risk" ";"  ";" "Low risk" _n
	file write textfile ";" "N events" ";" "IR (95% CI)" ";" "N events" ";" "IR (95% CI)" ";" "N events" ";" "IR (95% CI)" ";" "N events" ";" "IR (95% CI)" ";" "N events" ";" "IR (95% CI)" _n

	
/*******************************************************************************
#3. Number of events
*******************************************************************************/
	

	*Number of events
	foreach group in cohort qrisk1 qrisk0 hrisk1 hrisk0 {
		forvalues x=0/4 {
			if "`group'"=="cohort" unique newid if exgr==`x' & nevents==1 & `cond'==1
			if "`group'"=="qrisk1" unique newid if exgr==`x' & nevents==1 & qrisk3==1 & `cond'==1
			if "`group'"=="qrisk0" unique newid if exgr==`x' & nevents==1 & qrisk3==0 & `cond'==1
			if "`group'"=="hrisk1" unique newid if exgr==`x' & nevents==1 & hrisk==1 & `cond'==1
			if "`group'"=="hrisk0" unique newid if exgr==`x' & nevents==1 & hrisk==0 & `cond'==1
			
			local exgr_`group'_`x'=r(sum) 
			}  // end forvalues of exgr	

		if "`group'"=="cohort" unique newid if exgr==21 & nevents==1 & `cond'==1
			
/*******************************************************************************
#4. Season adjusted model
*******************************************************************************/		
		
		if "`group'"=="cohort" display in red "*******************Whole cohort `cond'*******************"	
		if "`group'"=="cohort" xi: xtpoisson nevents i.exgr i.season if `cond'==1, fe i(newid) offset(loginterval) irr
		if "`group'"=="qrisk1" display in red "*******************High QRISK3 `cond'********************"
		if "`group'"=="qrisk1" xi: xtpoisson nevents i.exgr i.season if qrisk3==1 & `cond'==1, fe i(newid) offset(loginterval) irr
		if "`group'"=="qrisk0" display in red "*******************Low QRISK3 `cond'*********************"
		if "`group'"=="qrisk0" xi: xtpoisson nevents i.exgr i.season if qrisk3==0 & `cond'==1, fe i(newid) offset(loginterval) irr
		if "`group'"=="hrisk1" display in red "*******************Hypertension `cond'*******************"	
		if "`group'"=="hrisk1" xi: xtpoisson nevents i.exgr i.season if hrisk==1 & `cond'==1, fe i(newid) offset(loginterval) irr
		if "`group'"=="hrisk0" display in red "*******************No hypertension `cond'****************"
		if "`group'"=="hrisk0" xi: xtpoisson nevents i.exgr i.season if hrisk==0 & `cond'==1, fe i(newid) offset(loginterval) irr
		forvalues x=1/4 {	
			local irr_`group'_`x'=exp(_b[_Iexgr_`x'])
			local ci1_`group'_`x'=exp(_b[_Iexgr_`x']-1.96*_se[_Iexgr_`x'])
			local ci2_`group'_`x'=exp(_b[_Iexgr_`x']+1.96*_se[_Iexgr_`x'])	
			} // end forvalues of exgr
		} // end foreach group

		
/*******************************************************************************
#5. Number of events by wave
*******************************************************************************/
	
	*Number of events
	foreach group in cohort qrisk1 qrisk0 hrisk1 hrisk0 {
	forvalues w=0/1 {
		forvalues x=0/4 {
			if "`group'"=="cohort" unique newid if exgr==`x' & nevents==1 & wave==`w' & `cond'==1
			if "`group'"=="qrisk1" unique newid if exgr==`x' & nevents==1 & qrisk3==1 & wave==`w' & `cond'==1
			if "`group'"=="qrisk0" unique newid if exgr==`x' & nevents==1 & qrisk3==0 & wave==`w' & `cond'==1
			if "`group'"=="hrisk1" unique newid if exgr==`x' & nevents==1 & hrisk==1 & wave==`w' & `cond'==1
			if "`group'"=="hrisk0" unique newid if exgr==`x' & nevents==1 & hrisk==0 & wave==`w' & `cond'==1
			
			local exgr_`group'_`x'_`w'=r(sum) 
			}  // end forvalues of exgr	
		} // end forvalues of wave	
			
/*******************************************************************************
#6. Season adjusted model by wave
*******************************************************************************/		
		
		forvalues w=0/1 {
		if "`group'"=="cohort" display in red "*******************Whole cohort `cond' wave `w'*******************"
		if "`group'"=="cohort" xi: xtpoisson nevents i.exgr i.season if wave==`w' & `cond'==1, fe i(newid) offset(loginterval) irr
		if "`group'"=="qrisk1" display in red "*******************High QRISK3 `cond' wave `w'********************"
		if "`group'"=="qrisk1" xi: xtpoisson nevents i.exgr i.season if qrisk3==1 & wave==`w' & `cond'==1, fe i(newid) offset(loginterval) irr
		if "`group'"=="qrisk0" display in red "*******************Low QRISK3 `cond' wave `w'*********************"
		if "`group'"=="qrisk0" xi: xtpoisson nevents i.exgr i.season if qrisk3==0 & wave==`w' & `cond'==1, fe i(newid) offset(loginterval) irr
		if "`group'"=="hrisk1" display in red "*******************Hypertension `cond' wave `w'*******************"
		if "`group'"=="hrisk1" xi: xtpoisson nevents i.exgr i.season if hrisk==1 & wave==`w' & `cond'==1, fe i(newid) offset(loginterval) irr
		if "`group'"=="hrisk0" display in red "*******************No hypertension `cond' wave `w'****************"
		if "`group'"=="hrisk0" xi: xtpoisson nevents i.exgr i.season if hrisk==0 & wave==`w' & `cond'==1, fe i(newid) offset(loginterval) irr
		
		forvalues x=1/4 {	
			local irr_`group'_`x'_`w'=exp(_b[_Iexgr_`x'])
			local ci1_`group'_`x'_`w'=exp(_b[_Iexgr_`x']-1.96*_se[_Iexgr_`x'])
			local ci2_`group'_`x'_`w'=exp(_b[_Iexgr_`x']+1.96*_se[_Iexgr_`x'])	
			} // end forvalues of exgr
		} // end forvalues of wave	
		} // end foreach group
		
			
/*******************************************************************************
#7. Output results to csv
*******************************************************************************/
		
	*Create label for risk period to add to output
	local label1 "1-7 days"
	local label2 "8-14 days"
	local label3 "15-28 days"
	local label4 "29-91 days"
	
	*Output results
	forvalues x=1/4 {
		file write textfile "`label`x''" ";" (`exgr_cohort_`x'') ";" %5.2f (`irr_cohort_`x'') " (" %4.2f (`ci1_cohort_`x'') "-" %4.2f (`ci2_cohort_`x'') ")" ";" (`exgr_qrisk1_`x'') ";" %5.2f (`irr_qrisk1_`x'') " (" %4.2f (`ci1_qrisk1_`x'') "-" %4.2f (`ci2_qrisk1_`x'') ")" ";" (`exgr_qrisk0_`x'') ";" %5.2f (`irr_qrisk0_`x'') " (" %4.2f (`ci1_qrisk0_`x'') "-" %4.2f (`ci2_qrisk0_`x'') ")" ";" (`exgr_hrisk1_`x'') ";" %5.2f (`irr_hrisk1_`x'') " (" %4.2f (`ci1_hrisk1_`x'') "-" %4.2f (`ci2_hrisk1_`x'') ")" ";" (`exgr_hrisk0_`x'') ";" %5.2f (`irr_hrisk0_`x'') " (" %4.2f (`ci1_hrisk0_`x'') "-" %4.2f (`ci2_hrisk0_`x'') ")" _n 
	} // end forvalues 
	file write textfile "Baseline" ";" (`exgr_cohort_0') ";" "ref" ";" (`exgr_qrisk1_0') ";" "ref" ";" (`exgr_qrisk0_0') ";" "ref" ";" (`exgr_hrisk1_0') ";" "ref" ";" (`exgr_hrisk0_0') ";" "ref" _n /*baseline period*/
	
	file write textfile "Wave 1" _n
	forvalues x=1/4 {
		file write textfile "`label`x''" ";" (`exgr_cohort_`x'_0') ";" %5.2f (`irr_cohort_`x'_0') " (" %4.2f (`ci1_cohort_`x'_0') "-" %4.2f (`ci2_cohort_`x'_0') ")" ";" (`exgr_qrisk1_`x'_0') ";" %5.2f (`irr_qrisk1_`x'_0') " (" %4.2f (`ci1_qrisk1_`x'_0') "-" %4.2f (`ci2_qrisk1_`x'_0') ")" ";" (`exgr_qrisk0_`x'_0') ";" %5.2f (`irr_qrisk0_`x'_0') " (" %4.2f (`ci1_qrisk0_`x'_0') "-" %4.2f (`ci2_qrisk0_`x'_0') ")" ";" (`exgr_hrisk1_`x'_0') ";" %5.2f (`irr_hrisk1_`x'_0') " (" %4.2f (`ci1_hrisk1_`x'_0') "-" %4.2f (`ci2_hrisk1_`x'_0') ")" ";" (`exgr_hrisk0_`x'_0') ";" %5.2f (`irr_hrisk0_`x'_0') " (" %4.2f (`ci1_hrisk0_`x'_0') "-" %4.2f (`ci2_hrisk0_`x'_0') ")" _n 
	} // end forvalues 
	file write textfile "Baseline" ";" (`exgr_cohort_0_0') ";" "ref" ";" (`exgr_qrisk1_0_0') ";" "ref" ";" (`exgr_qrisk0_0_0') ";" "ref" ";" (`exgr_hrisk1_0_0') ";" "ref" ";" (`exgr_hrisk0_0_0') ";" "ref" _n
	
	file write textfile "Wave 2" _n
	forvalues x=1/4 {
		file write textfile "`label`x''" ";" (`exgr_cohort_`x'_1') ";" %5.2f (`irr_cohort_`x'_1') " (" %4.2f (`ci1_cohort_`x'_1') "-" %4.2f (`ci2_cohort_`x'_1') ")" ";" (`exgr_qrisk1_`x'_1') ";" %5.2f (`irr_qrisk1_`x'_1') " (" %4.2f (`ci1_qrisk1_`x'_1') "-" %4.2f (`ci2_qrisk1_`x'_1') ")" ";" (`exgr_qrisk0_`x'_1') ";" %5.2f (`irr_qrisk0_`x'_1') " (" %4.2f (`ci1_qrisk0_`x'_1') "-" %4.2f (`ci2_qrisk0_`x'_1') ")" ";" (`exgr_hrisk1_`x'_1') ";" %5.2f (`irr_hrisk1_`x'_1') " (" %4.2f (`ci1_hrisk1_`x'_1') "-" %4.2f (`ci2_hrisk1_`x'_1') ")" ";" (`exgr_hrisk0_`x'_1') ";" %5.2f (`irr_hrisk0_`x'_1') " (" %4.2f (`ci1_hrisk0_`x'_1') "-" %4.2f (`ci2_hrisk0_`x'_1') ")" _n 
	} // end forvalues 
	file write textfile "Baseline" ";" (`exgr_cohort_0_1') ";" "ref" ";" (`exgr_qrisk1_0_1') ";" "ref" ";" (`exgr_qrisk0_0_1') ";" "ref" ";" (`exgr_hrisk1_0_1') ";" "ref" ";" (`exgr_hrisk0_0_1') ";" "ref" _n
	
	
	capture file close textfile 	

/*******************************************************************************
>> end loops	
*******************************************************************************/	
	
} /*end foreach suffix in datasets*/
} /*end foreach cond in outcome*/




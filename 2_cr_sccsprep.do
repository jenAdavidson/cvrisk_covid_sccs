/*=========================================================================
DO FILE NAME:		3_cr_sccsprep

AUTHOR:				Jennifer Davidson (edited do file provided by Kate 
					Mansfield)		
DATE:				23/06/2021
			
PURPOSE:			Create datasets needed for SCCS analysis with:
					interval variable, season, wave and severity
			
ADDITIONAL INFO:	File goes through each bit of info (cutpoint) needed 
					and creates a dataset containing: patid; season; ageband; 
					exposure interval; sex; eventdate; cv risk 
			
DATASET USED: 		final`suffix' (where suffix is lab confirmed or 
					clinically reported COVID-19)
													
DATASET CREATED:	analysis`cond'`suffix' (where file is season 
					breakdown and cond is outcome)
					
*=========================================================================*/

clear all
capture log close
log using "$logdir\sccsprep.log", replace

/*******************************************************************************
#1. loop through datasets i.e. main analysis and sensitivity analysis
*******************************************************************************/

* set up local macro containing suffixes for datasets to open
* and also to use for name of saved dataset
local datasets labconfirmsarscov2 clinreportcovid19
local outcome mace acs mi angina stroke hf arrhythmia

* loop through each dataset
foreach suffix of local datasets {
foreach cond of local outcome {
use "$datadir/`cond'`suffix'", clear		
gen suffix = "`suffix'"

/*******************************************************************************
*******************************************************************************
#A. MAKE dates DATASET - date at change in state of each var
*******************************************************************************
*******************************************************************************/	

/* 
Go through each variable and create a dataset containing the variable and date. 
Then append each of these datasets to each other. 
Using the completed dataset containing all dates:
	 - check for records that happen on the same day
	 - time-update all records using [_n+/-1]
	 - calculate age in days at date (used to create interval)
	 - calculate interval in days between changes in state (records)
	 - drop any records before indexdate and after end of FU
	 - identify whether MACE event happens within interval
*/
	
	
/*******************************************************************************
#2. Find start and end dates and save to file
*******************************************************************************/	
preserve	
keep newid indexdate enddate // only keep useful vars
		
* 2.1 date1 = start FU
gen long date1=indexdate
label var date1 "indexdate"
drop indexdate

* 2.2 date2 = end FU
gen long date2=enddate
label var date2 "end FU date"
drop enddate	
				
* 2.3 reshape
reshape long date, i(newid) j(type)
label var date "date at cutpoint"
format date %td
label var type "type of cutpoint, use with label var"
		
* 2.4 label as start/end dates
gen label=1
label var label "label to identify type of cutpt"
* create value labels for all values of the label var
label define lablelbl 1 "startend" 2 "period" 3 "season" 
label values label lablelbl
		
* 2.5 save
label data "dates and cutpoint types"
notes: cutpts and cutpoint types
save "$datadir/dates_`cond'`suffix'", replace
restore
	
	
/*******************************************************************************
#3. Generate exposure group cutpoints and reshape and save to dates file
*******************************************************************************/
preserve 
keep newid exposdate enddate // only keep useful vars

* 3.1 create cutpoints for each risk period
local counter=1 // counter var used to number cutpoints
gen long date`counter' = exposdate-14 // sets cutpt for pre-exposed time
label var date`counter' "end of unexposed pre-covid period"
local counter=`counter'+1
gen long date`counter' = exposdate // sets cut pt exposure date
label var date`counter' "end of pre-exposure period"
local counter=`counter'+1
gen long date`counter' = exposdate + 7 // sets end of risk period 1
label var date`counter' "end of risk period 1 (1 to 7 days)"
local counter=`counter'+1
gen long date`counter' = exposdate + 14 // sets end of risk period 2 
label var date`counter' "end of risk period 2 (8 to 14 days)"
local counter=`counter'+1
gen long date`counter' = exposdate + 28 // sets end of risk period 3 
label var date`counter' "end of risk period 3 (15 to 28 days)"
local counter=`counter'+1
gen long date`counter' = exposdate + 91 // sets end of risk period 4 
label var date`counter' "end of risk period 4 (29 to 91 days)"
local counter=`counter'+1
gen long date`counter' = enddate // sets cutpt for baseline time
label var date`counter' "end of unexposed post-covid period"

drop exposdate enddate
		
* 3.2 reshape
reshape long date, i(newid) j(type)
format date %td
label var date "date at cutpoint"
label var type "type of cutpoint, use with label var"
		
* 3.3 label as risk periods
gen label=2
label var label "label to identify type of cutpt"
label define lablelbl 1 "startend" 2 "period" 3 "season" 
label values label lablelbl

/*
1=end unexposed (prior infection)
2=end pre-exposed
3=end risk1
4=end risk2
5=end risk3
6=end risk4
7=end unexposed (post infection)
*/
	
* create exgr var
recode type 1=0 7=0 2=21 3=1 4=2 5=3 6=4 , gen(exgr)
		
label define exgr 0 "unexp" 21 "pre-exp" 1 "risk1" 2 "risk2" 3 "risk3" 4 "risk4" 
label values exgr exgr
label var exgr "exposure: baseline, pre-exp, risk periods"

* append to dataset containing start and end cutpoints
sort newid date
order newid date type exgr
append using "$datadir/dates_`cond'`suffix'"
save "$datadir/dates_`cond'`suffix'", replace
restore


/******************************************************************************
#6. Generate season
******************************************************************************/	
preserve
keep newid // only keep useful vars	
	/*
* 6.1 season in main analysis - first day of each season: warmer (Apr-Sep) 
*& colder (Oct-Mar)
local counter=1
forvalues yr = 2019/2021 { // do years around time period to ensure there are no empty values
foreach month in 4 10 {
	gen str = "1/`month'/`yr'"
	gen date`counter' = date(str, "DMY")
	label var date`counter' "date at start of `month'"
	format date`counter' %td
	drop str
	local counter=`counter'+1
	} // end foreach month
	} // end forvalues yr 
*/
gen date1=d(30sep2020)
format date1 %td
gen date2=d(31dec2020)
format date2 %td
	
* 6.2 reshape
reshape long date, i(newid) j(type)
label var date "date at cutpoint"
label var type "type of cutpoint, use with label var"
gen season = mod(type,2) 
label var season "season of year"
label define season 0 "colder" 1 "warmer"  
label values season season
		

* 6.3 label as risk periods
gen label=3
label var label "label to identify type of cutpt"
label define lablelbl 1 "startend" 2 "period" 3 "season" 
label values label lablelbl
			
* 6.4 append to dataset containing start and end cutpoints
order newid date type label 
append using "$datadir/dates_`cond'`suffix'"
save "$datadir/dates_`cond'`suffix'", replace
restore


/*******************************************************************************
*******************************************************************************
#B. MAKE static variables dataset 
*******************************************************************************
*******************************************************************************/	

/*******************************************************************************
#4. Age group
*******************************************************************************/	
preserve
keep newid dob indexdate // only keep useful vars
	
gen age=int(indexdate-dob)/365.25
replace age=round(age)
	
* create exposure var
gen ageband=1 if age <55
replace ageband=2 if age <65 & ageband==.
replace ageband=3 if age <75 & ageband==.
replace ageband=4 if ageband==.
label var ageband "ageband"
label define agebands 1 "40-54" 2 "55-64" 3 "65-74" 4 "75-84"
label values ageband agebands
save "$datadir/vars_`cond'`suffix'", replace
restore


/******************************************************************************
#5. Generate wave
******************************************************************************/	
preserve
keep newid exposdate // only keep useful vars	
	
* 5.1 last day of each wave: wave 1 (12 Mar-16 Aug) & wave 2 (17 Aug-31 Dec)
gen wave = 0 if exposdate<= d(16/08/2020)
replace wave = 1 if wave==.

label define wave 0 "Wave 1" 1 "Wave 2"
label values wave wave
		
merge 1:1 newid using "$datadir/vars_`cond'`suffix'", nogen
save "$datadir/vars_`cond'`suffix'", replace
restore

/******************************************************************************
#7. CV risk
******************************************************************************/	

preserve
keep newid 
merge 1:1 newid using "$datadirobj1/qrisk3", nogen keep(master match) keepusing(qrisk3 qrisk3_sens)
merge 1:1 newid using "$datadirobj1/qrisk2", nogen keep(master match) keepusing(qrisk2)
merge 1:1 newid using "$datadirobj1/hypertens", nogen keep(master match) keepusing(hrisk)
replace hrisk=0 if hrisk!=1
merge 1:1 newid using "$datadir/vars_`cond'`suffix'", nogen
save "$datadir/vars_`cond'`suffix'", replace
restore


/******************************************************************************
#8. Severity - estimated by hospitalisation
******************************************************************************/	
/*
preserve
import delimited "J:\EHR-Working\Jennifer\RaisedCVDRisk_COVID19\analysis\objective1_incidence\datafiles\linkeddata\hes_diagnosis_hosp_20_000135_request2_DM.txt", stringcols(1 2) 
merge m:1 patid using "J:\EHR-Working\Jennifer\RaisedCVDRisk_COVID19\analysis\objective1_incidence\datafiles\patid_master.dta"
keep if _merge==3
gen admidate2=date(admidate, "DMY")
format admidate2 %td
drop admidate
rename admidate2 admidate
gen discharged2=date(discharged, "DMY")
format discharged2 %td
drop discharged
rename discharged2 discharged
drop if discharged<d(31dec2019)
keep newid admidate discharged
save "J:\EHR-Working\Jennifer\RaisedCVDRisk_COVID19\analysis\objective3_sccs\datafiles\2020hospadmit.dta"
restore
*/

preserve
keep newid exposdate suffix // only keep useful vars

if suffix=="labconfirmsarscov2" { 
merge 1:1 newid using "$datadirobj1/chosp_labconfirmsarscov2", keep (master match) nogen
rename hosplabcov hosp
rename hosplabcovdate hospdate
}

else {
merge 1:1 newid using "$datadirobj1/chosp_clinreportcovid19", keep (master match) nogen
rename hospclincov hosp
rename hospclincovdate hospdate
}
		
* 8.1 add in COVID-19 hospitalisation records
merge 1:m newid using "$datadir/2020hospadmit", keep(master match) nogen	
replace admidate=. if hosp==.
replace discharged=. if hosp==.
duplicates drop
gen days=discharge-hospdate if hospdate<discharged
sort newid days
duplicates drop newid, force
gen hosplen=discharged-admidate

drop hospdate admidate discharged

* 8.2 merge to dataset containing cutpoint
merge 1:1 newid using "$datadir/vars_`cond'`suffix'", nogen
save "$datadir/vars_`cond'`suffix'", replace
restore


/*******************************************************************************
*******************************************************************************
#C. Make SCCS dataset
*******************************************************************************
*******************************************************************************/	

/*******************************************************************************
#9. Open cutpoint dates file 
*******************************************************************************/

use "$datadir/dates_`cond'`suffix'", clear
sort newid date	
	
	
/*******************************************************************************
#10. deal with records on the same day 
*******************************************************************************/
collapse (firstnm) label type season exgr, by(newid date)
by newid: replace exgr=exgr[_n+1] if exgr==0 & exgr[_n+1]>0 & exgr[_n+1]!=. & type==7
sort newid date
label var label "label to identify type of cutpt"
label var type "type of cutpoint, use with label var"
label var exgr "infection risk periods"
label var season "season of year"
label var date "date of changes in state"
	
label values label lablelbl
label values exgr exgr
label values season season
	
* identify if there are any further records on the same day
duplicates report newid date // all good


/*******************************************************************************
#11. add individual based data 
*******************************************************************************/	
* add main dataset and further variables created
merge m:1 newid using "$datadir/`cond'`suffix'", nogen
merge m:1 newid using "$datadir/vars_`cond'`suffix'", nogen
	
gen long cutp=(date-dob)


/*******************************************************************************
#12. propagate agebands through patient's record in ageband var 
*******************************************************************************/	
* need to update age for intervals where data is missing
sort newid date	
count if ageband>=.
local nmiss = r(N)
local nchange = 1 // counter var
while `nchange'>0{
	bysort newid: replace ageband = ageband[_n+1] if ageband>=.
	count if ageband>=.
	local nchange = `nmiss'-r(N)
	local nmiss = r(N)
	} /* end while `nchange'>0*/
	
	
/*******************************************************************************
#13. propagate wave through patient's record in season var 
*******************************************************************************/	
* need to update season for intervals where data is missing
sort newid date	
count if wave>=.
local nmiss = r(N)
local nchange = 1 // counter var
while `nchange'>0{
	bysort newid: replace wave = wave[_n+1] if wave>=.
	count if wave>=.
	local nchange = `nmiss'-r(N)
	local nmiss = r(N)
	} /* end while `nchange'>0*/	
	
	
/*******************************************************************************
#14. propagate season through patient's record in season var 
*******************************************************************************/	
* need to update season for intervals where data is missing
sort newid date	
count if season>=.
local nmiss = r(N)
local nchange = 1 // counter var
while `nchange'>0{
	bysort newid: replace season = season[_n+1] if season>=.
	count if season>=.
	local nchange = `nmiss'-r(N)
	local nmiss = r(N)
	} /* end while `nchange'>0*/

	
/*******************************************************************************
#15. propagate exgr through patient's record in exgr var 
*******************************************************************************/	
* need to update any exposure values for intervals that happen between different
* risk periods - for example, season could have changed during a risk period 
* therefore need to make sure that the indvidual is appropriately exposed
* deal with missing values of exposure
sort newid date
count if exgr>=.
local nmiss = r(N)
local nchange = 1 // counter var
while `nchange'>0{
	bysort newid: replace exgr = exgr[_n+1] if exgr>=.
	count if exgr>=.
	local nchange = `nmiss'-r(N)
	local nmiss = r(N)
	}/* end while `nchange'>0*/
replace exgr = 0 if exgr==.	
	
	
/*******************************************************************************
#16. Drop records before startdate and after enddate
*******************************************************************************/
drop if date<indexdate
replace exgr=exgr[_n+1] if exgr==0 & (exgr[_n+1]>=1 & exgr[_n+1]<=4)
drop if date>enddate
	
	
/*******************************************************************************
#17. Identify number of MACE events in each interval
*******************************************************************************/
* eventdate = age at MACE
gen long eventday=eventdate-dob
label var eventday "age in days at MACE"

* count the number of events in each interval
bysort newid: gen long nevents = 1 if eventday > cutp[_n-1]+0.5 & eventday <= cutp[_n]+0.5
collapse (sum) nevents, by(newid cutp date label type ageband wave season exgr sex indexdate exposdate eventdate enddate died diedfu hrisk qrisk3 qrisk3_sens qrisk2 hosp hosplen) // create one record /px / day / cut point type
replace nevents = 1 if indexdate==eventdate
label var nevents "number of MACE events in interval"
sort newid date

	
/*******************************************************************************
#18. Count number of days within each interval
*******************************************************************************/
*intervals
bysort newid: gen long interval = cutp[_n] - cutp[_n-1]
label var interval "duration of interval in days"
	

/*******************************************************************************
#19. Clean up:
- No longer need type or cutp
- can't do anything with zero intervals so drop
- missing intervals represent cut points for start of FU
- used to calculate first interval, no longer needed so drop
*******************************************************************************/
*fit issue of those with exactly 14 days between index and exposdate not having an interval
by newid: replace interval=1 if exgr[_n+1]==21 & cutp[_n+1]-cutp[_n]==14 & _n==1
*fit issue of those with exposure date on index
by newid: replace interval=1 if exgr==21 & interval==. & indexdate==exposdate & _n==1
drop type cutp
count if interval==0
count if interval==.
drop if interval==0 | interval==.

*drop patients who have no baseline time
unique newid
gen baseline=1 if exgr==0 
by newid: egen baselinecount=max(baseline)
drop if baselinecount==.

drop baseline baselinecount
/*
*below changes are because date used for cut is the start of the following season
gen month=month(date)
replace season=1 if label==4 & season==0 & month==10
replace season=0 if label==4 & season==1 & month==4
drop month
*/

replace sex=1 if sex==.

/*******************************************************************************
#20. Generate loginterval and save dataset
*******************************************************************************/	
* create loginterval for analysis
gen loginterval = log(interval)
label var loginterval "log of interval duration"
	
label data "SCCS `cond '`suffix' analysis dataset"
notes: SCCS `cond' `suffix' analysis dataset
	
* save
compress
save "$datadir/analysis_`cond'`suffix'", replace
	}/* end foreach cond in `outcome' */
	} /* end foreach suffix in `datasets' */
log close
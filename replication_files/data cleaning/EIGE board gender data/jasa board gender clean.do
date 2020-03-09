****************************************************
*Staggered Synthetic Control*
*Import and prepare board gender data*
*Data: EIGE board female ratio*
*Last updated 03/08/2020*
****************************************************

/*
Source: EIGE https://eige.europa.eu/gender-statistics/dgs/indicator/wmidm_bus_bus__wmid_comp_compbm/datatable
Purpose: clean data and format for matlab
*/

set more off

**step 1
clear
import excel "wmidm_bus_bus__wmid_comp_compbm.xlsx", sheet("stata") firstrow
rename GeographicregionTime country
*keep only B2, which is year end
keep country *B2

**step 2
*reshape
reshape long Y, i(country) j(temp) string
rename Y fratio
gen year=substr(temp,1,4)
destring year, replace
drop temp

**step 3
*bring in policy
merge m:1 country using policy
keep if _merge==3
drop _merge

**step 4
*check which countries to drop due to missing data
tab country if missing(fratio)
tab year
*after browsing, no missing data: 0 is truly 0%

**step 5
*clean up
*all drop 4 too early adopt
drop if inlist(country,"Norway","Spain","Finland","Iceland")

**step 6
*create treatment variables
gen t=year-Year
gen t_quota=year-Year if Policy=="Quota"
gen t_disclosure=year-Year if Policy=="Disclosure"

gen treat=year>=Year 
gen treatq=treat if Policy=="Quota"
replace treatq=0 if missing(treatq)
gen treatd=treat if Policy=="Disclosure"
replace treatd=0 if missing(treatd)

egen id=group(country)
save jasa_boardgendereige, replace

**step 7
use jasa_boardgendereige, clear
preserve
gen time=year
gen unit=country
sort id time
keep  id unit time treat treatq treatd fratio
order  id unit time treat treatq treatd fratio
export delimited using "data_boardgendereige.csv", replace
restore


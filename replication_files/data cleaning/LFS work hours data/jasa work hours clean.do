****************************************************
*Staggered Synthetic Control*
*Eurostat weekly work hours data*
*Last updated 03/08/2020*
****************************************************

/*
Source: Eurostat lfsq_ewhuis https://ec.europa.eu/eurostat/databrowser/view/lfsq_ewhuis/default/table
Purpose: clean data and format for matlab
*/

/*note on data variables
isco08 is occupation: 1 managers 2 professionals 3 technicians and associate pro 4 clerical support 5 service and sales
6 skilled agri forestry fishery 7 craft and related trades 8 plants and machine 9 elementary 0 armed forces
worktime is PT FT or TOTAL
*/

**set the working directory to the folder containing this do file
cap cd "~/replication_files/data cleaning/LFS work hours data"
set more off


**step 1 import
import excel "lfsq_ewhuis.xlsx", sheet("for stata") firstrow clear

**step 2
*turn to panel
reshape long Y, i(geo sex isco08 ) j(ym) string
rename Y value
drop worktime wstatus

**step 3
merge m:1 geo using policy
keep if _merge==3
drop _merge

*step 4 create treatment variables
gen year=substr(ym,1,4)
gen quarter=substr(ym,6,1)
destring year quarter, replace
gen time=yq(year, quarter)
gen Time=yq(Year, Quarter)
gen treat=time>=Time

**step 5
*reshape into value by gender and occupation
reshape wide value, i(country time isco08) j(sex) string
reshape wide value*, i(country time ) j(isco08) string
*for each value, create fratio
foreach i of numlist 1/9 { 
gen  valuefratioOC`i' = (valueFOC`i'/( valueFOC`i'+ valueMOC`i'))*100
}

**step 6
*check missing data
tab time if !missing( valuefratioOC2)
keep if time>=172
drop if time==239
tab country if missing( valuefratioOC2)
*drop countries with any missing data
drop if inlist(country,"Cyprus","Germany","Ireland","Switzerland")
*all drop 4 too early adopt
drop if inlist(country,"Norway","Spain","Finland","Iceland")
*tab country if missing( valuefratioOC1)
tab country if missing( valuefratioOC2)
*tab country if missing( valuefratioOC3)
*tab country if missing( valuefratioOC4)
*tab country if missing( valuefratioOC5)
*tab country if missing( valuefratioOC6)
*tab country if missing( valuefratioOC7)
*tab country if missing( valuefratioOC8)
tab country if missing( valuefratioOC9)

*
**step 7
*create treatment variables
gen t=time-Time
gen t_quota=time-Time if Policy=="Quota"
gen t_disclosure=time-Time if Policy=="Disclosure"

gen treatq=treat if Policy=="Quota"
replace treatq=0 if missing(treatq)
gen treatd=treat if Policy=="Disclosure"
replace treatd=0 if missing(treatd)

egen id=group(country)
save jasa_hr, replace

**step 8
*export to csv for matlab
use jasa_hr, clear
preserve
gen unit=country
sort id time
keep  id unit time treat treatq treatd *OC2 *OC9
order  id unit time treat treatq treatd *OC2 *OC9
export delimited using "data_hr.csv", replace
restore







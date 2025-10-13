****************************************************
*Staggered Synthetic Control*
*Eurostat FT employment data*
*Last updated 03/08/2020*

/*
Source: Eurostat lfsq_epgais https://ec.europa.eu/eurostat/databrowser/view/LFSQ_EPGAIS__custom_6036320/default/table
Purpose: clean data and format for matlab
*/

/*note on data variables
isco08 is occupation: 1 managers 2 professionals 3 technicians and associate pro 4 clerical support 5 service and sales
6 skilled agri forestry fishery 7 craft and related trades 8 plants and machine 9 elementary 0 armed forces
worktime is PT FT or TOTAL
*/

**set the working directory to the folder containing this do file
cap cd "~/replication_files/data cleaning/LFS full time employment data"
set more off

**step 1
clear
import excel "lfsq_epgais.xlsx", sheet("for stata") firstrow

**step 2
keep if age=="Y_GE15"
keep if worktime=="FT"
drop age
*turn to panel
reshape long Y, i(geo sex isco08 worktime ) j(ym) string
rename Y value

**step 3
*bring in policy
merge m:1 geo using policy
keep if _merge==3
drop _merge

**step 4 
*create treatment variables
gen year=substr(ym,1,4)
gen quarter=substr(ym,6,1)
destring year quarter, replace
gen time=yq(year, quarter)
gen Time=yq(Year, Quarter)
gen treat=time>=Time

**step 5
*reshape by gender and occupation
reshape wide value, i(country time isco08) j(sex) string
reshape wide value*, i(country time ) j(isco08) string
*for each value, create fratio
foreach i of numlist 1/9 { 
gen  valuefratioOC`i' = (valueFOC`i'/( valueFOC`i'+ valueMOC`i'))*100
}

**step 6
*check for missing data
tab time if !missing( valuefratioOC2)
keep if time>=172
tab country if missing( valuefratioOC2)
*5 countries with missing data
drop if inlist(country,"Croatia", "Cyprus","Germany","Ireland","Switzerland")
*drop 4 early adoptors
drop if inlist(country,"Norway","Spain","Finland","Iceland")
*check missing values
*tab country if missing( valuefratioOC1)
tab country if missing( valuefratioOC2)
*tab country if missing( valuefratioOC3)
*tab country if missing( valuefratioOC4)
*tab country if missing( valuefratioOC5)
*tab country if missing( valuefratioOC6)
*tab country if missing( valuefratioOC7)
*tab country if missing( valuefratioOC8)
tab country if missing( valuefratioOC9)

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
save jasa_ft, replace

**step 8
*export to csv for matlab
use jasa_ft, clear
preserve
sort id time
drop unit
gen unit=country
keep  id unit time treat treatq treatd *OC2 *OC9
order  id unit time treat treatq treatd *OC2 *OC9
export delimited using "data_ft.csv", replace
restore


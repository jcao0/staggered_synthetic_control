****************************************************
*Staggered Synthetic Control*
*Plot DID for board gender result*
*Data: EIGE board female ratio*
*Last updated 03/08/2020*
****************************************************
/*
plot dynamics treatment effects using generalized-diff-in-diff
*/
set more off
*figure formatting
*graph
grstyle init
grstyle set plain, horizontal grid
grstyle set legend 10, inside nobox
grstyle set symbol

use jasa_boardgendereige, clear

preserve
gen time=t
replace time=-6 if time<=-6 & !missing(time)
replace time=6 if time>=6 & !missing(time) 
*create time relative to treat indicator
*let time_9 be t=-1, time_10 is first year of treatment
levelsof time, local(time)
foreach time in `time' {
	local count=`time'+10
	gen treat_time_`count'=0 
	replace treat_time_`count'=1 if time==`time'
}
drop treat_time_9
global treat_time treat_time_*

*main regression
reghdfe fratio $treat_time, absorb(id year) vce(cluster id) nocon

*quick plot of the coefficients 
coefplot, vert recast(line)

*formally plotting the time graph
*save betas and SEs
foreach v of varlist $treat_time {
      qui gen b_`v' = _b[`v']
      qui gen se_`v' = _se[`v']
}

gen b_treat_time_9=0
gen se_treat_time_9=0

*keep only one row
keep in 1
keep b_* se_*

*reshape
qui gen i=.
reshape long b_ se_, i(i) j(x) s
drop i
gen time=substr(x,12,.)
destring time, replace
replace time = time-10

*generate SE upper and lower bounds
qui gen u = b_ + 1.96*se_
qui gen l = b_ - 1.96*se_

*plot results
twoway rcap u l time || scatter b_ time, xline(-0.5) yline(0) ttitle("Event Time") ylabel(-10 (5) 25) xlabel(-6(1)6, angle(vertical)) legend(off) 
*title("Board Gender Policy Effect Time Plot" "country FE and year FE") 
graph export figureA1.png, replace
restore

****************************************************
*Staggered Synthetic Control*
*Table 1 Summary Stats*
*Last updated 03/08/2020*
****************************************************
/*
Purpose: to create summary stats in table 1
*/
set matsize 11000

*board female ratio data
use jasa_boardgender, clear
label var fratio "Board female ratio (in percentage)"
outreg2 using boardgendereige , replace sum(detail) keep(fratio) eqkeep(N mean sd p25 p50 p75)  tex(fragment) dec(2) label

*full time employment data
use jasa_ft, clear
label var valuefratioOC2 "Professional female ratio (in percentage)"
label var valueFOC2 "Professional female (thousands)"
label var valueMOC2 "Professional male (thousands)"
label var valuefratioOC9 "Non-professional female ratio (in percentage)"
label var valueFOC9 "Non-professional female (thousands)"
label var valueMOC9 "Non-professional male (thousands)"
order valuefratioOC2 valueFOC2 valueMOC2   valuefratioOC9  valueFOC9 valueMOC9  
local var valuefratioOC2 valueFOC2 valueMOC2   valuefratioOC9  valueFOC9 valueMOC9  
outreg2 using ft , replace sum(detail) keep(`var') eqkeep(N mean sd p25 p50 p75)  tex(fragment) dec(2) label

*work hours data
use jasa_hr, clear
label var valuefratioOC2 "Professional female ratio (in percentage)"
label var valueFOC2 "Professional female"
label var valueMOC2 "Professional male"
label var valuefratioOC9 "Non-professional female ratio (in percentage)"
label var valueFOC9 "Non-professional female "
label var valueMOC9 "Non-professional male"
order valuefratioOC2 valueFOC2 valueMOC2   valuefratioOC9  valueFOC9 valueMOC9  
local var valuefratioOC2 valueFOC2 valueMOC2   valuefratioOC9  valueFOC9 valueMOC9  
outreg2 using hr , replace sum(detail) keep(`var') eqkeep(N mean sd p25 p50 p75)  tex(fragment) dec(2) label

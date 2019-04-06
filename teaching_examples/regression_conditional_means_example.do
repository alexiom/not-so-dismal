clear all
log using "/Users/alexander/Desktop/example_log.smcl", replace

*===================================================================
* 				Generate Example Data
*===================================================================

/*
IMPORTANT:
This part of the code simply generates random data that we can use
in this example. Feel free to skip this as it will mean nothing to
you.
*/


qui{

global N = 40000
global mean_income_1 = 35000
global mean_income_2 = 37000
global mean_income_3 = 32000
global mean_income_4 = 36500

set seed 123
set obs $N

gen id = _n
gen obs_score = runiform()
sort obs_score 

egen group_id = cut(obs_score), at(0 0.25 0.5 0.75 1)
sort id
drop obs_score

replace group_id = 1 if group_id==0
replace group_id = 2 if group_id==0.25
replace group_id = 3 if group_id==0.5
replace group_id = 4 if group_id==0.75

gen income_1 = rnormal($mean_income_1 ,10000)
gen income_2 = rnormal($mean_income_2 ,10000)
gen income_3 = rnormal($mean_income_3 ,10000)
gen income_4 = rnormal($mean_income_4 ,10000)

gen income = income_1
replace income = income_2 if group_id==2
replace income = income_3 if group_id==3
replace income = income_4 if group_id==4

drop income_*

rename group_id wbho
}

*===================================================================
* Method 1a: simple sample mean by group using "if" statements
*===================================================================

summarize income if wbho == 1
summarize income if wbho == 2
summarize income if wbho == 3
summarize income if wbho == 4

*===================================================================
* Method 1b: simple sample mean by group using bysort syntax
*===================================================================

bysort wbho: sum income

*===================================================================
* Method 2: regression
*===================================================================

reg income i.wbho


*Show our coefficient estimates
matrix list e(b)

/* 
  Notice how 1b.wbho is 0, this is because group 1 is the base group
  group (this is why it has a b after the 1) and so it's average  
  income is the value of the constant.
*/

*Calculate actual means from the coefficients

*Average Income Group 1
display _b[_cons] 
*Average Income Group 2
display _b[_cons]+_b[2.wbho]
*Average Income Group 3
display _b[_cons]+_b[3.wbho]
*Average Income Group 4
display _b[_cons]+_b[4.wbho]

/*
You can eyeball these and see that they are identical from the
estimates in methods 1a and 1b. However, to see both estimates side 
by side, see the following code.

IMPORTANT:
You don't have to worry about the code as you likely won't follow
it. However, know that the first estimates to be displayed will be 
the means by group when like in methods 1a/1b and the second 
estimates will be from the regression method.
*/


*===================================================================
* 			Side-by-Side Comparison
*===================================================================

*make matrix MM2 from coefficient matrix
qui{
matrix MM2 = (0,0,0,0)
matrix MM2[1,1] = _b[_cons] 
matrix MM2[1,2] = _b[_cons] + _b[2.wbho] 
matrix MM2[1,3] = _b[_cons] + _b[3.wbho] 
matrix MM2[1,4] = _b[_cons] + _b[4.wbho] 

matrix rownames MM2 = "average income"
matrix colnames MM2 = "group 1" "group 2" "group 3" "group 4"
}

qui{
by wbho: eststo: estpost summarize income
}
esttab, cells("mean(fmt(0))") noobs nonumbers /// 
	mtitles("group 1" "group 2" "group3" "group4") 

matlist  MM2, format(%9.2g)

log close

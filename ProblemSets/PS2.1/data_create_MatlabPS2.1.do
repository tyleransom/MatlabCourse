clear all
set mem 2g
set maxvar 32000
version 11.0
capture cd "/afs/econ.duke.edu/home/t/tmr17/Teaching/PS2.1"
capture log close
set more off
log using data_create_MatlabPS2.1.log, replace

* Clean data from master.dta to class example
* variables desired:
* male
* AFQT
* Mhgc
* hgc
* exper
* Diploma
* AA
* BA
* wage
* activity

use "/afs/econ.duke.edu/home/t/tmr17/RA_Peter/NLSY97_Data/master.dta", clear

generat activity = 1 if Primary_activity<2
replace activity = 2 if Primary_activity>2 & Primary_activity<5
replace activity = 3 if activity >=.

replace Highest_Grade_Completed = . if Highest_Grade_Completed>25
bys ID: gen t = _n
bys ID: ipolate Highest_Grade_Completed t, gen(hgc) epolate

replace Bio_mother_highest_educ = . if Bio_mother_highest_educ>25
rename  Bio_mother_highest_educ Mhgc

drop if ASVAB>=.
drop if hgc>=.
drop if Mhgc>=.

drop BA AA Diploma GED

gen Diploma = year>=Diploma_year
gen AA      = year>=AA_year
gen BA      = year>=BA_year

replace BA = 0 if Diploma==0
replace AA = 0 if Diploma==0

zscore ASVAB
ren    z_ASVAB AFQT

gen in_work = activity==2
bys ID: gen exper = sum(in_work)

merge 1:1 ID year using /afs/econ.duke.edu/data/vjh3/WageReturns/Data/y97/y97_master, keepusing(activityOct) nogen

keep if age_now>=20 & age_now<=24

keep ID age_now male t AFQT Mhgc hgc exper Diploma AA BA log_wage activity activityOct
ren age_now age
drop t
bys ID: gen t = _n
xtset ID t
drop if ID==3082 | ID==3234 | ID==4530 | ID==6574 | ID==7353
xtset

*** Analysis quickly
mlogit activity male AFQT Mhgc hgc exper Diploma AA BA, base(3)
mlogit activityOct male AFQT Mhgc hgc exper Diploma AA BA, base(6)
reg    log_wage male hgc exper Diploma AA BA
reg    log_wage male AFQT Mhgc hgc exper Diploma AA BA
xtreg  log_wage male AFQT Mhgc hgc exper Diploma AA BA, fe
xtreg  log_wage male AFQT Mhgc hgc exper Diploma AA BA, re

gen schwrk = (activityOct==2)
gen workPT = (activityOct==3)

program normal
version 11.0
args lnf Xb sigma
quietly replace `lnf'=ln(normalden($ML_y1, `Xb', `sigma'))
end

ml model lf normal (log_wage=male AFQT Mhgc hgc exper Diploma AA BA schwrk workPT) /sigma
ml max

reg log_wage male AFQT Mhgc hgc exper Diploma AA BA schwrk workPT

drop schwrk workPT
reshape wide age log_wage activity activityOct hgc exper Diploma AA BA, i(ID) j(t)
outsheet ID male AFQT Mhgc hgc? exper? Diploma? AA? BA? log_wage? activity? using nlsy97.csv, comma nol replace
outsheet ID male AFQT Mhgc hgc? exper? Diploma? AA? BA? log_wage? activityOct? using nlsy97m.csv, comma nol replace

log close
exit

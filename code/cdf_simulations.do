net install disco, from("https://raw.githubusercontent.com/Davidvandijcke/DiSCos_stata/dev/src/") replace

global maindir = "/Users/davidvandijcke/University of Michigan Dropbox/David Van Dijcke/Flo_GSRA/sj_replication"
global figs = "${maindir}/results/figs"

************************************************************
* 1. Set up the dataset
************************************************************
clear all
set seed 12349   // For reproducibility; you can set any seed you like

* We want 4 distributions, each with 100 draws -> total = 400
set obs 500

************************************************************
* 2. Create an identifier for distribution "id" 
*    and a uniform random variable "u"
************************************************************
gen id = .
gen u = runiform()


forvalues i = 1/5 {
    // Lower bound
    local first = ((`i' - 1) * 100) + 1
    // Upper bound
    local last  = `i' * 100

    // Now replace id = i from observation `first' to `last'
    replace id = `i' in `first'/`last'
}


************************************************************
* 3. Generate the outcome variable "y" 
*    according to different probability distributions
************************************************************
gen y = .

* ---- Distribution #1: More mass on 1 ----
* Suppose p(1)=0.4, p(2)=0.2, p(3)=0.2, p(4)=0.2
replace y = 1 if id==1 & u<0.4
replace y = 2 if id==1 & u>=0.4 & u<0.6
replace y = 3 if id==1 & u>=0.6 & u<0.8
replace y = 4 if id==1 & u>=0.8

* ---- Distribution #2: More mass on 2 ----
* Suppose p(1)=0.2, p(2)=0.4, p(3)=0.2, p(4)=0.2
replace y = 1 if id==2 & u<0.2
replace y = 2 if id==2 & u>=0.2 & u<0.6
replace y = 3 if id==2 & u>=0.6 & u<0.8
replace y = 4 if id==2 & u>=0.8

* ---- Distribution #3: More mass on 3 ----
* Suppose p(1)=0.2, p(2)=0.2, p(3)=0.4, p(4)=0.2
replace y = 1 if id==3 & u<0.2
replace y = 2 if id==3 & u>=0.2 & u<0.4
replace y = 3 if id==3 & u>=0.4 & u<0.8
replace y = 4 if id==3 & u>=0.8

* ---- Distribution #4: More mass on 4 ----
* Suppose p(1)=0.2, p(2)=0.2, p(3)=0.2, p(4)=0.4
replace y = 1 if id==4 & u<0.2
replace y = 2 if id==4 & u>=0.2 & u<0.4
replace y = 3 if id==4 & u>=0.4 & u<0.6
replace y = 4 if id==4 & u>=0.6

* ---- Target distribution: discrete uniform ----
replace y = 1 if id==5 & u<0.25
replace y = 2 if id==5 & u>=0.25 & u<0.5
replace y = 3 if id==5 & u>=0.5 & u<0.75
replace y = 4 if id==5 & u>=0.75

************************************************************
* 4. Clean up and show a sample of the data
************************************************************
drop u
sort id
list in 1/20



gen t = 1

***************************************************************
* (A) Barycenter approach (default)
***************************************************************
disco y id t, idtarget(5) t0(2)  // Donors are id=1..4, target is id=5

* Store the relevant quantile/cdf matrices from e()
matrix Bary_Qtarget  = e(quantile_t)     // True (target) QF
matrix Bary_Qsynth   = e(quantile_synth) // Barycenter synthetic QF
matrix Bary_Ctarget  = e(cdf_t)          // True (target) CDF
matrix Bary_Csynth   = e(cdf_synth)      // Barycenter synthetic CDF

***************************************************************
* (B) Mixture approach
***************************************************************
disco y id t, idtarget(5) t0(2) mixture

matrix Mix_Qtarget  = e(quantile_t)      // True (target) QF
matrix Mix_Qsynth   = e(quantile_synth)  // Mixture synthetic QF
matrix Mix_Ctarget  = e(cdf_t)           // True (target) CDF
matrix Mix_Csynth   = e(cdf_synth)       // Mixture synthetic CDF




***************************************************************
* Plot quantile: True vs. Barycenter vs. Mixture
***************************************************************
preserve
clear

* 1) 100 grid points for quantile
set obs 100
gen tau = (_n - 1)/(100 - 1)

* 2) Populate with the stored matrices
svmat Bary_Qtarget, name(QT_true)
svmat Bary_Qsynth,  name(QT_bary)
svmat Mix_Qsynth,   name(QT_mix)

rename QT_true1  q_true
rename QT_bary1  q_bary
rename QT_mix1   q_mix

* 3) Overlay lines + custom legend, axis labels, etc.
twoway ///
    (line q_true tau,  lcolor(black) lwidth(medium)) ///
    (line q_bary tau,  lcolor(blue)  lpattern(dash)      lwidth(medium)) ///
    (line q_mix  tau,  lcolor(red)   lpattern(dash_dot)  lwidth(medium)), ///
    /// Legend formatting
    legend(order(1 "True (Target)" 2 "Qtile-Based" 3 "CDF-Based") ///
           ring(0) position(5) cols(1) size(large)) ///
    /// Axis tick labels
    xlabel(0(.2)1, grid labsize(vlarge)) ///
    ylabel(1(1)4, angle(horiz) grid labsize(vlarge)) ///
    /// Axis titles
    xtitle("Quantile", size(vlarge)) ///
    ytitle("Y", size(vlarge)) ///
    name(qf_compare, replace)

graph export "${figs}/mixture_quantile.pdf", replace

restore




***************************************************************
* Plot CDF: True vs. Barycenter vs. Mixture
***************************************************************
preserve
clear

* 1) 100 grid points for CDF
set obs 100
gen grid_val = 1 + (_n - 1)*(4 - 1)/(100 - 1)

* 2) Populate with stored matrices
svmat Bary_Ctarget, name(CT_true)
svmat Bary_Csynth,  name(CT_bary)
svmat Mix_Csynth,   name(CT_mix)

rename CT_true1   c_true
rename CT_bary1   c_bary
rename CT_mix1    c_mix

* 3) Overlay lines + custom legend, axis labels, etc.
twoway ///
    (line c_true  grid_val, lcolor(black) lwidth(medium)) ///
    (line c_bary  grid_val, lcolor(blue)  lpattern(dash)     lwidth(medium)) ///
    (line c_mix   grid_val, lcolor(red)   lpattern(dash_dot) lwidth(medium)), ///
    /// Legend formatting
    legend(order(1 "True (Target)" 2 "Qtile-Based" 3 "CDF-based") ///
           ring(0) position(5) cols(1) size(large)) ///
    /// Axis tick labels
    xlabel(1(1)4, grid labsize(vlarge)) ///
    ylabel(0(.2)1, angle(horiz) grid labsize(vlarge)) ///
    /// Axis titles
    xtitle("y", size(vlarge)) ///
    ytitle("Pr(Y ≤ y)", size(vlarge)) ///
    name(cdf_compare, replace)

graph export "${figs}/mixture_cdf.pdf", replace

restore



*******************************************************************************
* Manually compute each donor's quantile function (id=1..4) plus target (id=5),
* then overlay them in one graph.
*******************************************************************************

save temp, replace


********************************************************************************
* 2) Build a 100 x 5 matrix "donorQ"
********************************************************************************
matrix donorQ = J(100, 5, .)

forvalues d = 1/5 {
    // Reload original data each iteration
    use temp, clear
    keep if id == `d'
    
    // Stata immediate command: _pctile
    // which puts the quantiles in r(c_1) ... r(c_100)
    _pctile y, nquantiles(101)

    forvalues j = 1/100 {
        matrix donorQ[`j', `d'] = r(r`j')
    }
}

********************************************************************************
* 3) Turn that matrix into a new dataset of 100 rows
********************************************************************************
clear
set obs 100
gen tau = (_n - 1)/(100 - 1)   // in [0,1]

svmat donorQ, name(Q_)
rename Q_1 Q1
rename Q_2 Q2
rename Q_3 Q3
rename Q_4 Q4
rename Q_5 Qtarget

********************************************************************************
* 4) Plot them on a single graph
********************************************************************************
twoway ///
    (line Q1 tau,      lcolor(red)    lwidth(medium)) ///
    (line Q2 tau,      lcolor(blue)   lwidth(medium)) ///
    (line Q3 tau,      lcolor(green)  lwidth(medium)) ///
    (line Q4 tau,      lcolor(orange) lwidth(medium)) ///
    (line Qtarget tau, lcolor(black)  lwidth(thick)), ///
    legend(order(1 "Donor #1" 2 "Donor #2" 3 "Donor #3" 4 "Donor #4" 5 "Target") ///
           ring(0) pos(5) cols(1)) ///
    xlabel(0(.2)1, grid labsize(large)) ///
    ylabel(1(1)4, angle(horiz) grid labsize(large)) ///
    xtitle("Quantile", size(large)) ytitle("Y", size(large)) ///
    name(donors_qf, replace)
	
graph export "${figs}/mixture_raw_quantiles.pdf", replace

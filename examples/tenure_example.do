* tenure_example.do
* Replicates the tenure (quantile-based) results in
* Gunsilius and Van Dijcke, "disco: Distributional Synthetic Controls".
*
* Run this file from the examples/ directory of the replication package:
*   . cd examples
*   . do tenure_example.do
* The log and figures are written to the working directory.

version 18.0
clear all
set more off

* install the bundled disco package (skip these two lines if disco is
* already installed, e.g., from the SJ archives)
net install disco, from("`c(pwd)'/../src") replace
mata: mata mlib index

capture log close _all
log using "tenure_example.log", replace text

* load and inspect the anonymized tenure data
use "tenure_anonymized.dta", clear
list in 1/5, ab(20)

* main disco command
disco y_col id_col time_col, idtarget(2) t0(3) agg("quantileDiff") ///
	seed(12143) g(10) m(100) ci boots(300)

* top 5 weights by company name
disco_weight id_col company_name, n(5)

* summary table of quantile effects
disco_estat summary

* plot quantile effects
disco_plot, title(" ") ytitle("Difference in Tenure (Days)") hline(0) ///
	color1(black) cicolor(gs12) scheme(sj)
graph export "tenure_quantileDiff.pdf", replace

* plot observed vs. synthetic quantile functions
disco_plot, title(" ") ytitle("Tenure (Days)") agg("quantile") ///
	yrange(0 3000) color1(black) color2(black) scheme(sj)
graph export "tenure_quantile.pdf", replace

log close

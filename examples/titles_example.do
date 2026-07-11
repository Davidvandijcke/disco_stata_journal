* titles_example.do
* Replicates the title (CDF-based, mixture) results in
* Gunsilius and Van Dijcke, "disco: Distributional Synthetic Controls".
*
* Run this file from the examples/ directory of the replication package:
*   . cd examples
*   . do titles_example.do
* The log and figures are written to the working directory.

version 18.0
clear all
set more off

* install the bundled disco package (skip these two lines if disco is
* already installed, e.g., from the SJ archives)
net install disco, from("`c(pwd)'/../src") replace
mata: mata mlib index

capture log close _all
log using "titles_example.log", replace text

* load and inspect the anonymized titles data
use "titles_anonymized.dta", clear
list in 1/5, ab(20)

* main disco command: CDF-based (mixture) approach for the categorical outcome
disco y_col id_col time_col, idtarget(2) t0(3) agg("cdfDiff") seed(12143) ///
	mixture g(10) m(10) ci boots(300)

* top 5 weights by company name
disco_weight id_col company_name, n(5)

* summary table of CDF effects
disco_estat summary

* plot CDF effects
disco_plot, title(" ") ytitle("Change in CDF") hline(0) categorical ///
	color1(gs10) cicolor(gs4) scheme(sj)
graph export "title_cdfDiff.pdf", replace

* plot synthetic vs. treated CDF
disco_plot, title(" ") ytitle("CDF") agg("cdf") color1(black) color2(black) scheme(sj)
graph export "title_cdf.pdf", replace

log close

clear

net install disco, from("https://raw.githubusercontent.com/Davidvandijcke/DiSCos_stata/dev/src/") replace
net install gzimport, from(https://raw.githubusercontent.com/mdroste/stata-gzimport/master/) replace

net install disco, from("/Users/davidvandijcke/University of Michigan Dropbox/David Van Dijcke/Flo_GSRA/stata_repo/src") replace



//**************************
// Set Paths
//**************************
global maindir = "/Users/davidvandijcke/University of Michigan Dropbox/David Van Dijcke/Flo_GSRA/sj_replication"
global figs = "${maindir}/results/figs"
global dataIn = "${maindir}/data/in"
global dataOut = "${maindir}/data/out"



//**************************
// Load data and anonymize
//**************************

local process_data = 0

if `process_data' {
	// Tenure

	clear all
	gzimport delimited using "${dataIn}/tenure_stata.csv.gz", clear


	// Define the range for uniform noise
	quietly: summarize y_col
	gen noise = runiform() * 0.1 * r(mean) - 0.05 * r(mean) // Noise in [-5%, +5%] of mean

	// Add noise to y_col
	gen y_col_anonymized = y_col + noise

	// Ensure the minimum remains at least 30 (approximately one month)
	replace y_col_anonymized = max(y_col_anonymized, 30)

	// Summary statistics to check the effect
	summarize y_col y_col_anonymized

	replace y_col = y_col_anonymized

	keep id_col company_name y_col time_col
	
	
	* First, create a percentile rank within each id-time group
	bysort id_col time_col: egen pctile_rank = pctile(y_col), p(90)

	* Keep only observations that are below or equal to the 90th percentile
	keep if y_col <= pctile_rank

	* Clean up
	drop pctile_rank
	
	gen random = runiform()
	sort random
	drop random

	// Save the anonymized dataset as a Stata file
	save "${dataOut}/tenure_anonymized.dta", replace


	// Title
	clear all
	gzimport delimited using "${dataIn}/titles_stata.csv.gz", clear

	// Generate noise with position-specific probabilities
	gen noise = 0
	// Reduced noise (5% chance each way) for lowest and highest positions
	replace noise = -1 if runiform() < 0.05 & (y_col == 1 | y_col == 10)
	replace noise = 1 if runiform() < 0.05 & noise == 0 & (y_col == 1 | y_col == 10)

	// Moderate noise (15% chance each way) for middle positions
	replace noise = -1 if runiform() < 0.1 & y_col > 1 & y_col < 10 & noise == 0
	replace noise = 1 if runiform() < 0.1 & y_col > 1 & y_col < 10 & noise == 0

	// Apply noise and enforce bounds
	gen y_col_anonymized = y_col + noise
	replace y_col_anonymized = 1 if y_col_anonymized < 1
	replace y_col_anonymized = 10 if y_col_anonymized > 10


	// Check the effect
	tabulate y_col y_col_anonymized

	// Update and keep relevant columns
	replace y_col = y_col_anonymized
	keep id_col company_name y_col time_col


	// Save the anonymized dataset as a Stata file
	save "${dataOut}/titles_anonymized.dta", replace

}



//**************************
// Reproduce tenure results
//**************************

// disco <- DiSCo(grouped, id_col.target = id_col.target, t0 = t0, q_max=0.9, G = G, M=M, num.cores=20,
//                cl=0.95,uniform=TRUE, permutation = TRUE, CI = TRUE, boots = 1000, simplex=TRUE, seed=30, qtype=7,
//                qmethod=NULL) # seed 5


cd "$dataOut"

capture sjlog close
sjlog using "tenure_analysis", replace 

// load and inspect data
use "tenure_anonymized.dta", clear 
list in 1/5, ab(20) 


// run disco commands
disco y_col id_col time_col, idtarget(2) t0(3) agg("quantileDiff") ///
	seed(12143) g(10) m(100) ci boots(300) 
disco_weight id_col company_name, n(5)
disco_estat summary
disco_plot, title(" ") ytitle("Difference in Tenure (Days)") hline(0) scheme("stsj")
graph export "${figs}/tenure_quantileDiff.pdf", replace

// plot quantile functions separately
disco_plot, title(" ") ytitle("Tenure (Days)") agg("quantile") /// 
	yrange(0 3000) scheme("stsj") // vline(0)
graph export "${figs}/tenure_quantile.pdf", replace


capture sjlog close





//**************************
// Reproduce title results
//**************************

// disco <- DiSCo(grouped, id_col.target = id_col.target, t0 = t0, q_max=0.9, G = G, M=M, num.cores=20,
//                cl=0.95,uniform=TRUE, permutation = TRUE, CI = TRUE, boots = 1000, simplex=TRUE, seed=30, qtype=7,
//                qmethod=NULL) # seed 5

cd "$dataOut"

sjlog using "titles_analysis.log", replace 

use "titles_anonymized.dta", clear 
list in 1/5, ab(20) 

// main disco command
disco y_col id_col time_col, idtarget(2) t0(3) agg("cdfDiff") seed(12143) /// 
	mixture g(10) m(10) ci boots(300)
	 
// plot top 5 weights
disco_weight id_col company_name, n(5)

// plot summary table
disco_estat summary

// plot CDF effects
disco_plot, title(" ") ytitle("Change in CDF") hline(0) categorical /// 
	scheme("stsj") color("bluishgray") 
graph export "${figs}/title_cdfDiff.pdf", replace

// plot synthetic vs. treated CDF
disco_plot, title(" ") ytitle("Tenure (Days)") agg("cdf") scheme("stsj")

graph export "${figs}/title_cdf.pdf", replace

capture sjlog close, replace





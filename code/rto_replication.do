clear

net install disco, from("https://raw.githubusercontent.com/Davidvandijcke/DiSCos_stata/dev/src/") replace
net install gzimport, from(https://raw.githubusercontent.com/mdroste/stata-gzimport/master/) replace


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

	// Save the anonymized dataset as a Stata file
	save "${dataOut}/tenure_anonymized.dta", replace


	// Title
	clear all
	gzimport delimited using "${dataIn}/titles_stata.csv.gz", clear

	// Add random perturbation to y_col
	gen noise = round(runiform() * 2 - 1) // Random noise: -1, 0, or 1
	gen y_col_anonymized = y_col + noise

	// Ensure y_col_anonymized stays within the valid range [1, 10]
	replace y_col_anonymized = cond(y_col_anonymized < 1, 1, y_col_anonymized) // Cap minimum at 1
	replace y_col_anonymized = cond(y_col_anonymized > 10, 10, y_col_anonymized) // Cap maximum at 10

	// Summary statistics to check the effect
	tabulate y_col y_col_anonymized
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


disco y_col id_col time_col, idtarget(2) t0(3) agg("quantileDiff") qmax(0.9) seed(1242)    // ci boots(500)




* Create a temporary file to store the matches
tempname memhold

* Store the id-company name pairs we need
postfile `memhold' str32 company_name double weight using weights_matched.dta, replace


* Loop through the cids and weights to match with company names
forvalues i = 1/`=colsof(e(cids))' {
    local id = e(cids)[1,`i']
    local w = e(weights)[1,`i']
    * Get company name for this id
	qui levelsof company_name if id_col == `id', local(company) clean
    post `memhold' ("`company'") (`w')
}
postclose `memhold'

preserve
* Load the matched data
use weights_matched.dta, clear


gsort -weight

keep in 1/5
replace weight = round(weight, 0.0001)

texsave * using top_weights.tex, replace 

restore 

* tenure_data_creation.do
* Creates the releasable file examples/tenure_anonymized.dta from the
* confidential raw tenure data (Van Dijcke, Gunsilius, and Wright 2026).
*
* FOR PROVENANCE ONLY. This script cannot be run without the confidential
* input file tenure_stata.csv.gz, which is not distributed with this
* package. It is documented here so that the construction of the anonymized
* dataset is transparent; it is not part of the SJ online distribution.

version 18.0
clear all
set more off

* ---- EDIT: directory holding the confidential raw inputs ----
local dataIn "PATH-TO-CONFIDENTIAL-INPUTS"

gzimport delimited using "`dataIn'/tenure_stata.csv.gz", clear

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

* Scramble the row order
gen random = runiform()
sort random
drop random

// Save the anonymized dataset next to the example do-files
save "../examples/tenure_anonymized.dta", replace

* titles_data_creation.do
* Creates the releasable file examples/titles_anonymized.dta from the
* confidential raw titles data (Van Dijcke, Gunsilius, and Wright 2026).
*
* FOR PROVENANCE ONLY. This script cannot be run without the confidential
* input file titles_stata.csv.gz, which is not distributed with this
* package. It is documented here so that the construction of the anonymized
* dataset is transparent; it is not part of the SJ online distribution.

version 18.0
clear all
set more off

* ---- EDIT: directory holding the confidential raw inputs ----
local dataIn "PATH-TO-CONFIDENTIAL-INPUTS"

gzimport delimited using "`dataIn'/titles_stata.csv.gz", clear

// Generate noise with position-specific probabilities
gen noise = 0
// Reduced noise (5% chance each way) for lowest and highest positions
replace noise = -1 if runiform() < 0.05 & (y_col == 1 | y_col == 10)
replace noise = 1 if runiform() < 0.05 & noise == 0 & (y_col == 1 | y_col == 10)

// Moderate noise (10% chance each way) for middle positions
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

// Save the anonymized dataset next to the example do-files
save "../examples/titles_anonymized.dta", replace

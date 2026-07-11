clear all

capture log using build_disco_mlib, replace

mata:
mata clear
mata set matastrict off
end

// Load all Mata functions into memory
do disco_utils.mata
do quadprog.mata

// Now create and populate the Mata library
mata:
mata mlib create ldisco, replace

// Add each function defined in disco_utils.mata
mata mlib add ldisco CI_out()
mata mlib add ldisco disco_out()
mata mlib add ldisco get_unique()
mata mlib add ldisco disco_prob_grid()
mata mlib add ldisco disco_quantile_sorted()
mata mlib add ldisco disco_quantile_points()
mata mlib add ldisco disco_quantile()
mata mlib add ldisco disco_cdf_sorted()
mata mlib add ldisco cdf_builder()
mata mlib add ldisco cdf_at_points()
mata mlib add ldisco disco_boot_sorted()
mata mlib add ldisco disco_build_cells()
mata mlib add ldisco disco_precompute()
mata mlib add ldisco disco_run_from_store()
mata mlib add ldisco disco_solve_weights()
mata mlib add ldisco disco_mixture_weights()
mata mlib add ldisco disco_full_run()
mata mlib add ldisco disco_ratio_from_diff()
mata mlib add ldisco disco_compute_ratio()
mata mlib add ldisco disco_permutation_test()
mata mlib add ldisco disco_bootstrap_CI()
mata mlib add ldisco disco_wrapper()
mata mlib add ldisco disco_ci_wrapper()
mata mlib add ldisco compute_summary_stats()
mata mlib add ldisco solve_quadprog()
end


capture log close

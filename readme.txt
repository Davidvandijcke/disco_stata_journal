Package name: disco stata journal replication package

DOI:

Title: Stata Journal Replication Package for disco

Author 1 name: Florian Gunsilius
Author 1 from: Emory University
Author 1 email: fgunsil@emory.edu

Author 2 name: David Van Dijcke
Author 2 from: University of Virginia
Author 2 email: eju5az@virginia.edu

Author 3 name:
Author 3 from:
Author 3 email:

Author 4 name:
Author 4 from:
Author 4 email:

Author 5 name:
Author 5 from:
Author 5 email:

Help keywords: DiSCo, Stata, cdf, quantile, synthetic control, replication

File list:
- src: the disco Stata package. The do-files install it automatically with
  net install. It contains the command and help files (disco.ado, disco_estat.ado,
  disco_plot.ado, disco_weight.ado, quadprog.ado and their .sthlp files), the
  compiled Mata library (ldisco.mlib), the platform plugins
  (quadprog_mata_{mac,mac_intel,linux,win}), the Mata source (disco_utils.mata,
  quadprog.mata, build_disco_mlib.do), and the package index (disco.pkg, stata.toc).
- code:
    * code/cdf_simulations.do: generates the figures that illustrate the mixture
      (CDF-based) approach. Self-contained; creates its own simulated data.
    * code/rto_replication.do: replicates the tenure and title results, including the
      disco_estat summary table and the figures, using the anonymized data in data/out.
- data:
    * data/out/tenure_anonymized.dta: anonymized, perturbed tenure data from
      Van Dijcke, Gunsilius, and Wright (2024).
    * data/out/titles_anonymized.dta: anonymized, perturbed title data from the same
      paper.
- results:
    * results/cdf_simulations.log: Stata log produced by code/cdf_simulations.do.
    * results/rto_replication.log: Stata log produced by code/rto_replication.do.
    * results/figs/*.pdf: the figures that appear in the article (monochrome).
    * results/paper/: the compiled article PDF.
- readme.txt: this file.

How to reproduce:
1. Install Stata 18 or later.
2. Open code/cdf_simulations.do and code/rto_replication.do and set the global
   "root" at the top of each file to the directory where you unpacked this package.
3. Run code/cdf_simulations.do, then code/rto_replication.do. Each file installs the
   disco package from src/, writes its log to results/, and writes its figures to
   results/figs/.

Notes:
The original employment records are confidential, so the data shipped here are a
perturbed version: uniform noise is added to the outcomes, which preserves the shape
of the distributions while keeping the dataset shareable. Estimated weights and point
estimates match those in the article; bootstrap standard errors can differ slightly
across Stata versions because the bootstrap draws depend on the random-number stream.

References:
Van Dijcke, David, Florian Gunsilius, and Austin Wright. 2024. "Return to Office and
the Tenure Distribution." arXiv preprint arXiv:2405.04352.

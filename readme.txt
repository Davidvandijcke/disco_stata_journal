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
- src: the disco Stata package. The example do-files install it automatically with
  net install. It contains the command and help files (disco.ado, disco_estat.ado,
  disco_plot.ado, disco_weight.ado, quadprog.ado and their .sthlp files), the
  compiled Mata library (ldisco.mlib), the platform plugins
  (quadprog_mata_{mac,mac_intel,linux,win}), the Mata source (disco_utils.mata,
  quadprog.mata, build_disco_mlib.do), and the package index (disco.pkg, stata.toc).
- examples: self-contained example do-files and the datasets they use. Each file
  can be run out of the box from the examples/ directory and writes its log and
  figures to the working directory.
    * examples/tenure_example.do: replicates the tenure (quantile-based) results in
      the article, including the disco_estat summary table and the figures. Starts
      from "use tenure_anonymized.dta, clear".
    * examples/titles_example.do: replicates the title (CDF-based, mixture) results.
      Starts from "use titles_anonymized.dta, clear".
    * examples/cdf_simulations.do: generates the figures that illustrate the mixture
      (CDF-based) approach. Creates its own simulated data.
    * examples/tenure_anonymized.dta, examples/titles_anonymized.dta: anonymized,
      perturbed employment data from Van Dijcke, Gunsilius, and Wright (2026).
- data_creation (not part of the SJ archive): for provenance only. The scripts
  that construct the anonymized datasets from the confidential raw data. They
  cannot be run without the confidential inputs, which are not distributed. They
  are available in the public repository accompanying this package at
  https://github.com/Davidvandijcke/disco_stata_journal.
- results: the logs and figures produced by running the three example files, and
  the compiled article (results/paper/).
- readme.txt: this file.

How to reproduce:
1. Install Stata 18 or later.
2. Change to the examples/ directory of this package.
3. Run, in any order:
     do cdf_simulations.do
     do tenure_example.do
     do titles_example.do
   Each file is self-contained: it installs the bundled disco package from ../src
   (skip those lines if disco is already installed), loads its dataset from the
   working directory, and writes its log and figures to the working directory.
4. The results/ folder contains the logs and figures we obtained from these runs.

Notes:
The original employment records are confidential, so the data shipped here are a
perturbed version: uniform noise is added to the outcomes, which preserves the shape
of the distributions while keeping the dataset shareable. Estimated weights and point
estimates match those in the article; the bootstrap confidence intervals are
reproducible via the seed() option used in the example files.

References:
Van Dijcke, David, Florian Gunsilius, and Austin Wright. 2026. "Return to Office and
the Tenure Distribution." Review of Economics & Statistics, conditionally accepted.

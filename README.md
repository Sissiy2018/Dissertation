# Dissertation
Aim: develop statistical methods using genealogical inference techniques to estimate recombination rates (and detect variations in these rates over time)

1. **Simulation of Genealogies**: Using the ARGON2 software, simulate genealogies under different recombination rates.

2. **Estimation of Recombination Rate Maps**: Develop methods to build estimators of recombination rate maps based on these simulated genealogies.

3. **Application to Simulated Data**: Apply the developed estimators to genealogies inferred from simulated data. Compared with LDhat and Pyrho.

5. **Potential Application to Real Data**: If time permits, extend the analysis to genealogies inferred from real datasets, such as those from the 1000 Genomes Project.

First stage:
- simulate using msprime (https://tskit.dev/msprime/docs/stable/intro.html) to obtain both simulated genotypes (in VCF) and the tree sequence object (in tskit format).
  e.g. a constant population size, say 20,000 and 1000 diploid samples.
  tools useful for checking and manipulating VCF files later (bcftools: https://samtools.github.io/bcftools/howtos/index.html).

- simulate with a varying recombination rate (https://tskit.dev/msprime/docs/stable/rate_maps.html).
  the real human recombination rate (can be found here: https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/). should use the GRCh38 version

- set up and run pyrho (https://github.com/popgenmethods/pyrho?tab=readme-ov-file) using your simulated data.
  have a plot like Fig.1a in https://www.science.org/doi/10.1126/sciadv.aaw9206, obtaining the spearman correlation between true recombination map and inferred recombination map and compare our results to that of pyrho.

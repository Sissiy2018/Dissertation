# Dissertation
Aim: develop statistical methods using genealogical inference techniques to estimate recombination rates (and detect variations in these rates over time)

1. **Simulation of Genealogies**: Using the ARGON2 software, simulate genealogies under different recombination rates.

2. **Estimation of Recombination Rate Maps**: Develop methods to build estimators of recombination rate maps based on these simulated genealogies.

3. **Application to Simulated Data**: Apply the developed estimators to genealogies inferred from simulated data. Compared with LDhat and Pyrho.

5. **Potential Application to Real Data**: If time permits, extend the analysis to genealogies inferred from real datasets, such as those from the 1000 Genomes Project.

- simulate using msprime (https://tskit.dev/msprime/docs/stable/intro.html) to obtain both simulated genotypes (in VCF) and the tree sequence object (in tskit format).
  - e.g. a constant population size, say 20,000 and 1000 diploid samples.
  - tools useful for checking and manipulating VCF files later (bcftools: https://samtools.github.io/bcftools/howtos/index.html).

- simulate with a varying recombination rate (https://tskit.dev/msprime/docs/stable/rate_maps.html).
  - the real human recombination rate (can be found here: https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/). -- the GRCh38 version

- set up and run pyrho (https://github.com/popgenmethods/pyrho?tab=readme-ov-file) using the simulated data.
  - obtain a plot like Fig.1a in https://www.science.org/doi/10.1126/sciadv.aaw9206, obtaining the spearman correlation between true recombination map and inferred recombination map and compare ARG results to that of pyrho.

- reconstruct ARGs from the simulated data with arg_needle and get recombination rates along the geneome. calculate the there spearman correlation with the truth (with different window sizes).


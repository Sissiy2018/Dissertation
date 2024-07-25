# This example script is a walk through of a typical use for pyrho
# We will infer a recombination map using 10 arbitrary diploid individuals from
# ACB for the first 200 SNPs of chromosome 1.

# This script assumes that you've successfully installed pyrho.

# All of the pyrho commands have a number of options.  If anything
# is unclear, running pyrho <command> -h will return a list
# of options and their meaning for that command.

# Here all of the options are chosen such that this example should run in
# a few minutes on a single core on a laptop computer.
# I will point out places where in practice you would want to
# change the options to be a bit more thorough.

# If you encounter any issues or bugs with this script or pyrho in general
# please do not hesitate to open an issue at 
# https://github.com/popgenmethods/pyrho/issues



# Before we do anything else, we need to pre-compute a lookup table for our
# 10 diploid individuals.  For large samples this step is the most expensive step both
# computationally and memory-wise.  We recommend always using the --approx
# flag for sample sizes greater than say 20.  Also note that the memory
# usage can be substantial for sample sizes in the hundreds (~100G of RAM
# are required).

pyrho make_table -n 20 -N 25 --mu 1.25e-8 --logfile . --outfile ACB_n_20_N_40_lookuptable.hdf \
	--approx --smcpp_file ACB_pop_sizes.csv --decimate_rel_tol 0.1

pyrho make_table -n 20 -N 30 --mu 1e-8 --logfile . --outfile constant_n_20_N_30_lookuptable.hdf \
	--approx --smcpp_file constant_pop_sizes.csv

pyrho make_table --samplesize <n> --approx  --moran_pop_size <N> \
--numthreads <par> --mu <mu> --outfile <outfile> \
--popsizes <size1>,<size2>,<size3> --epochtimes <breakpoint1>,<breakpoint2>

pyrho make_table -n 20 -N 30 --mu 1e-8 --logfile constant_n_20_N_30_ldpop_run.log \
    --outfile constant_n_20_N_30_lookuptable.hdf --approx --numthreads 8 --popsizes 200,200 --epochtimes 0

pyrho make_table -n 50 -N 75 --mu 1e-8 --logfile constant_n50_N75_Ne500_ldpop_run.log \
    --outfile constant_n50_N75_Ne500_lookuptable.hdf --approx --numthreads 8 --popsizes 500,500 --epochtimes 0

pyrho make_table -n 100 -N 125 --mu 1e-8 --logfile constant_n100_N125_Ne1000_ldpop_run.log \
    --outfile constant_n100_N125_Ne1000_lookuptable.hdf \
	--approx --numthreads 8 --popsizes 1000,1000 --epochtimes 0

# in the pyrho article, n ∈ {20,40,60,80,100,120}

pyrho make_table -n 100 -N 125 --mu 1e-8 --logfile constant_n100_N125_Ne10000_ldpop_run.log \
    --outfile constant_n100_N125_Ne10000_lookuptable.hdf \
	--approx --numthreads 8 --popsizes 10000,10000 --epochtimes 0

pyrho make_table -n 50 -N 75 --mu 1e-8 --logfile constant_n50_N75_Ne10000_ldpop_run.log \
    --outfile constant_n50_N75_Ne10000_lookuptable.hdf \
	--approx --numthreads 8 --popsizes 10000,10000 --epochtimes 0

pyrho make_table -n 200 -N 256 --mu 1e-8 --logfile constant_n250_N256_Ne2000_ldpop_run.log \
    --outfile constant_n200_N256_Ne2000_lookuptable.hdf \
	--approx --numthreads 8 --popsizes 2000,2000 --epochtimes 0


# This command will output a lot of information while it is computing the lookup table.
# If you wanted to store that output in a different file, you could set logfile to be
# some name.

# The result will be a lookup table stored in ACB_n_20_N_40_lookuptable.hdf
# This lookup table will take the demography stored in ACB_pop_sizes.csv into account.
# I generated ACB_pop_sizes.csv by running smc++ to make a model file
# and then running smc++ plot --csv ACB_pop_sizes.pdf model.final.json

# Note that we set n=20 because we have 10 diploids --> 20 haploids
# Using the approx flag, we want to set N as high as computationally reasonable
# but to make this example run faster we set it to only 25. It is usually good
# to shoot for about N being 25-50% larger than n.

# decimate_rel_tol will smooth out small changes in the smc++ population
# sizes.  This will generally speed things up a bit, but the larger
# decimate_rel_tol is, the less closely the lookup table will match the
# inferred demography.




# Now that we have a table, it would be good to find hyper-parameter settings
# that work well for this demography.  We do this using pyrho hyperparam:

pyrho hyperparam -n 20 --mu 1.25e-8 --blockpenalty 50,100 \
	--windowsize 25,50 --logfile . --tablefile ACB_n_20_N_40_lookuptable.hdf \
	--num_sims 3 \
	--smcpp_file ACB_pop_sizes.csv --outfile ACB_hyperparam_results.txt

pyrho hyperparam -n 20 --mu 1e-8 --blockpenalty 15,20,25,30,35,40,45,50 \
	--windowsize 30,40,50,60,70,80,90 --logfile . --tablefile constant_n_20_N_30_lookuptable.hdf \
	--num_sims 5 \
	--smcpp_file constant_pop_sizes.csv --outfile constant_n20_hyperparam_results.txt

pyrho hyperparam -n 20 --mu 1e-8 --blockpenalty 15,20,25,30,35,40,45,50 \
	--windowsize 30,40,50,60,70,80,90 --logfile constant_n20_N30_Ne200_hyperparam.log --tablefile constant_n20_N30_Ne200_lookuptable.hdf \
	--num_sims 5 --ploidy 1 --numthreads 8 \
	--popsizes 200,200 --epochtimes 0 --outfile constant_n20_hyperparam_results.txt
# n20_N30_Ne200: 15, 30

pyrho hyperparam -n 50 --mu 1e-8 --blockpenalty 15,20,25,30,35,40,45,50 \
	--windowsize 30,40,50,60,70,80,90 --logfile constant_n50_N75_Ne500_hyperparam.log --tablefile constant_n50_N75_Ne500_lookuptable.hdf \
	--num_sims 5 --ploidy 1 --numthreads 8 \
	--popsizes 500,500 --epochtimes 0 --outfile constant_n50_N75_Ne500_hyperparam_results.txt
# n50_N75_Ne500: 50, 90

pyrho hyperparam -n 100 --mu 1e-8 --blockpenalty 15,20,25,30,35,40,45,50 \
	--windowsize 30,40,50,60,70,80,90 --logfile constant_n100_N125_Ne1000_hyperparam.log --tablefile constant_n100_N125_Ne1000_lookuptable.hdf \
	--num_sims 5 --ploidy 1 --numthreads 8 \
	--popsizes 1000,1000 --epochtimes 0 --outfile constant_n100_N125_Ne1000_hyperparam_results.txt
# n100_N125_Ne1000: 15, 90 or 50, 50

# (w, l) ∈ {30,40,50,60,70,80,90} × {15,20,25,30,35,40,45,50}

pyrho hyperparam -n 100 --mu 1e-8 --blockpenalty 15,20,25,30,35,40,45,50 \
	--windowsize 30,40,50,60,70,80,90 --logfile constant_n100_N125_Ne10000_hyperparam.log --tablefile constant_n100_N125_Ne10000_lookuptable.hdf \
	--num_sims 5 --ploidy 1 --numthreads 8 \
	--popsizes 10000,10000 --epochtimes 0 --outfile constant_n100_N125_Ne10000_hyperparam_results.txt
# n100_N125_Ne10000: 15, 40

pyrho hyperparam -n 50 --mu 1e-8 --blockpenalty 15,20,25,30,35,40,45,50 \
	--windowsize 30,40,50,60,70,80,90 --logfile constant_n50_N75_Ne10000_hyperparam.log --tablefile constant_n50_N75_Ne10000_lookuptable.hdf \
	--num_sims 5 --ploidy 1 --numthreads 8 \
	--popsizes 10000,10000 --epochtimes 0 --outfile constant_n50_N75_Ne10000_hyperparam_results.txt
# n50_N75_Ne10000: 15, 40

# This will search over different settings of the hyperparameters (blockpenalty and
# windowsize).  It will simulate <num_sims> 1Mb chunks and compute some statistics
# of how well the optimization works on those chunks.
# It is obviously better to set <num_sims> to be as high as feasible to have
# better estimates of the accuracy for the different hyperparameter settings.
# It is also obviously better to search over a large, fine grid of blockpenalties
# and windosizes, but note that doing so incurs a computational cost.
# Note, however, that this step is optional.  You may certainly jump straight
# to pyrho optimize using the default value of the hyperparameters.  We found
# that tuning the hyperparameters can result in slightly better performance
# but there is a wide range of hyperparameter values that work well.

# There is some randomness involved in the simulator, and so your results may differ.
# When I ran this command and looked at ACB_hyperparam_results.txt
# it was clear that windowsize = 50 and blockpenalty = 50 was the best
# setting for the hyperparameters.




bgzip your_file.vcf
zcat  xxx.vcf|cut -f1-10|less

# Now that we have some sensible hyper-parameters, we are finally ready to
# infer some fine-scale recombination rates for out dataset.  We do this
# using pyrho optimize:

pyrho optimize --tablefile ACB_n_20_N_40_lookuptable.hdf \
	--vcffile ACB_chr_1_subset.vcf.gz \
	--outfile ACB_chr_1_subset.rmap \
	--blockpenalty 50 --windowsize 50 \
	--logfile .

pyrho optimize --tablefile constant_n_20_N_30_lookuptable.hdf \
	--vcffile varmutGRCh38all_200pop_10dip.vcf.gz \
	--outfile constant_chr_1_first1M_subset.rmap \
	--blockpenalty 50 --windowsize 50 \
	--logfile .

pyrho optimize --tablefile constant_n20_N30_Ne200_lookuptable.hdf \
	--vcffile chr22_20Mb_10_Ne200.vcf \
	--outfile constant_n20_N30_Ne200_map.rmap \
	--blockpenalty 15 --windowsize 30 --ploidy 1 \
	--logfile constant_n20_N30_Ne200_optmize.log

pyrho optimize --tablefile constant_n50_N75_Ne500_lookuptable.hdf \
	--vcffile chr22_25_Ne500.vcf \
	--outfile constant_n50_N75_Ne500_map.rmap \
	--blockpenalty 50 --windowsize 90 --ploidy 1 \
	--logfile constant_n50_N75_Ne500_optmize.log

pyrho optimize --tablefile constant_n100_N125_Ne1000_lookuptable.hdf \
	--vcffile chr22_50_Ne1000.vcf \
	--outfile constant_n100_N125_Ne1000_map.rmap \
	--blockpenalty 50 --windowsize 50 --ploidy 1 \
	--logfile constant_n100_N125_Ne1000_optmize.log

pyrho optimize --tablefile constant_n100_N125_Ne10000_lookuptable.hdf \
	--vcffile chr22_50_Ne10000.vcf \
	--outfile constant_n100_N125_Ne10000_map.rmap \
	--blockpenalty 50 --windowsize 60 --ploidy 1 \
	--logfile constant_n100_N125_Ne10000_optmize.log

pyrho optimize --tablefile constant_n50_N75_Ne10000_lookuptable.hdf \
	--vcffile chr22_25_Ne10000.vcf \
	--outfile constant_n50_N75_Ne10000_map.rmap \
	--blockpenalty 15 --windowsize 40 --ploidy 1 \
	--logfile constant_n50_N75_Ne10000_optmize.log

# Note that this VCF only contains about 80kb and so the optimization
# should take about 0.1 seconds, so don't be alarmed if
# it terminates extremely quickly.

# The resulting recombination map is stored in ACB_chr_1_subset.rmap.
# Running the above should result in every interval inferred to have
# a recombination rate of about 8.8968e-09



# pyrho also has a utility to compute the distribution of r^2, a measure
# of linkage-disequillibrium from the lookup tables. This can be of interest
# to see if inferred demographies or recombination maps fit the data well.
# Note that the distribution of r^2 depends heavily on the demographic
# history.

pyrho compute_r2 --quantiles .25,.5,.75 --compute_mean --samplesize 20 \
	--tablefile ACB_n_20_N_40_lookuptable.hdf \
	--outfile ACB_r2.txt

# This command will compute the 25th, 50th, and 75th percentiles of the
# distribution of r^2, as well as its mean and store them in ACB_r2.txt.

# The distribution of r^2 is very sensitive to low-frequency alleles,
# which are also the most difficult to genotype accurately.  Therefore
# it is sometimes desirable to compute the distribution of r^2
# conditioned on the minor allele frequency being above some cutoff.
# This can be achieved by using the --MAFcut <thresh> option,
# where only pairs of sites where both sites have a minor allele frequency
# above thresh are considered.

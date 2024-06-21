import msprime
import numpy as np
import pandas as pd
import argparse
import gzip
import arg_needle_lib
import timeancestry as tac


def export_output(ts, name):
    a = arg_needle_lib.tskit_to_arg(ts)
    arg_needle_lib.serialize_arg(a, name + ".argn")

def sim_one_const(Ne, sample_size, length):
    # Create a demographic model for a constant-size population
    demographic_model = msprime.Demography()
    demographic_model.add_population(initial_size=Ne)
    # Simulate the tree sequence
    ts = msprime.sim_ancestry(samples=sample_size,
                                        demography=demographic_model,
                                        sequence_length = length,
					recombination_rate = 1e-8,
                                        ploidy=2)
    # Simulate mutations
    return ts

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate trees')
    parser.add_argument("-s", '--sample', dest="nsample", required=True, help="Number of haploid samples, which is not the same as the diploid sample size required by msprime")
    parser.add_argument("-l", '--len', dest="len", required=True)
    parser.add_argument("-Ne", '--Ne', dest="ne", required=True)
    parser.add_argument("-ndeme", "--ndeme", dest="ndeme", default="1")
    parser.add_argument("-o", '--output', dest="output", required=True)
    
    args = vars(parser.parse_args())
    print(args)
    ne = int(args['ne'])
    nsample = int(args['nsample'])
    length = int(args['len'])

    ts = sim_one_const(ne, nsample, length)
    export_output(ts, args['output'])


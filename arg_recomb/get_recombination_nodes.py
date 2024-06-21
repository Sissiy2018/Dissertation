import arg_needle_lib
import argparse
import numpy as np
import pandas as pd


def get_recombination(arg_file):
   arg = arg_needle_lib.deserialize_arg(arg_file)
   arg.populate_children_and_roots()

   child_id = []
   parent_id = []
   parent_start = []
   child_height = []
   parent_height = []

   for i in arg.node_ids():
       # Find all non-root nodes
       if len(arg.node(i).parent_starts()) != 0:
           parent_start.extend(arg.node(i).parent_starts())
           for parent_edge in arg.node(i).parent_edges():
               child_id.append(parent_edge.child.ID)
               parent_id.append(parent_edge.parent.ID)
               child_height.append(parent_edge.child.height)
               parent_height.append(parent_edge.parent.height)
    
   df = pd.DataFrame({
          'child_id': np.array(child_id), 
          'parent_id': np.array(parent_id),
          'parent_height': np.array(parent_height),
          'child_height': np.array(child_height),
          'parent_start': np.array(parent_start)
          })
   return(df)
   



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate trees')
    parser.add_argument("-a", '--arg', dest="arg_file", required=True, help="Number of haploid samples, which is not the same as the diploid sample size required by msprime")
    parser.add_argument("-o", '--output', dest="output", required=True)
    
    args = vars(parser.parse_args())
    arg_file = args['arg_file']
    output = args['output']

    nrecomb = get_recombination(arg_file)
    nrecomb.to_csv(output + ".csv", index=False)



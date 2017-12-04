import cooler
import hic_basics
import argparse

parser = argparse.ArgumentParser(description='Run intrachromosomal clustering algorithm.')
parser.add_argument('-cool',required=True,help="Cooler file containing Hi-C data.")
parser.add_argument('-outdir',required=True,help="Output dir, where 'clusters_chr*.txt' files will be generated.")
parser.add_argument('-k',required=True,help="Number of clusters.")
params = parser.parse_args()

cf = cooler.Cooler(params.cool)

hic_basics.cluster_compartments(cf=cf,k=params.k,chrlist=cf.chromnames,outdir=params.outdir)

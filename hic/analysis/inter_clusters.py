import cooler
import hic_basics
import argparse

parser = argparse.ArgumentParser(description='Run intrachromosomal clustering algorithm.')
parser.add_argument('-cool',required=True,help="Cooler file containing Hi-C data.")
parser.add_argument('-input',required=True,help="Intrachromosomal clusters file (all chromosomes in a single file).")
parser.add_argument('-out',required=True,help="Output file with global clusters.")
parser.add_argument('-k',required=True,help="Number of clusters.")
params = parser.parse_args()

cf = cooler.Cooler(params.cool)

hic_basics.interchromosomal_clusters(cf=cf,k=params.k,cluster_file=params.input,out_file=params.out)

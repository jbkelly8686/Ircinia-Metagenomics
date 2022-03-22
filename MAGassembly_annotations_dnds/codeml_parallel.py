#!/usr/bin/env python

#author: David E. Carlson david.carlson@stonybrook.edu

import argparse
import sys
import glob
from Bio.Phylo.PAML import codeml
import os
import multiprocessing
import shutil

parser = argparse.ArgumentParser(description='Take tree as input and run codeml in parallel over several alignments/trees')

parser.add_argument('--phylo_dir', required=True, dest='phylo_dir', action='store')
parser.add_argument('--begin', required=True, help='Beginning of range of files to work with', action='store')
parser.add_argument('--end', required=True, help='End of range of files to work with', action='store')
parser.add_argument('--threads', type=int, required=False, default=1, help='Number of threads to use. Default is 1', action='store')
args=parser.parse_args()

indir = args.phylo_dir
begin = int(args.begin)
end = int(args.end)
threads = int(args.threads)
# find genes trees that are finished running
seqs = glob.glob(indir + '/*.iqtree')

#print(seqs)

seq_ids = [seq.split('.')[0].split('/')[-1]for seq in seqs]
#seq_ids = [seq.split('.')[0] for seq in seqs]


#print(seq_ids)

def run_codeml(id):

	aln = "phylip_aln/" + id + ".phy"
	tree = 'phylo_out/' + id + '.treefile'
	outfile = 'codeml_out/' + id + '.cml.model0.out'
	working = 'temp/' + id + '_working'
	
	#print(aln, tree, outfile)
	
	# create blank outfile if it doesn't exist
	if not os.path.exists(outfile):
		open(outfile, 'a').close()
	
	with open(outfile) as f:
		if 'Time used' in f.read():
			print(f"codeml analysis already finished for gene {id}")
			sys.stdout.flush()
		else:
			print(f"Running codeml model 0 for gene {id}")
			sys.stdout.flush()
			cml = codeml.Codeml(alignment = aln, tree = tree,
				out_file = outfile, working_dir = working)
			
			cml.set_options(seqtype = 1)
			cml.set_options(runmode = 0)
			cml.set_options(model = 0)
			cml.set_options(kappa = 2)
			cml.set_options(fix_kappa = 0)
			cml.set_options(fix_omega = 0)
			cml.set_options(omega = 0.1)
			cml.set_options(fix_blength = 1)
			cml.set_options(CodonFreq = 2)
			cml.set_options(ndata = 1)
			cml.set_options(NSsites = [0])
			cml.set_options(icode = 2)
			cml.set_options(alpha = 0)
			cml.set_options(fix_alpha = 1)
			cml.set_options(ncatG = 8)
			cml.set_options(cleandata = 0)
			cml.run(verbose = True)
			
#run_codeml('K00801')

# paralelize over inputs

p = multiprocessing.Pool(processes=threads)
for id in seq_ids[begin:end+1]:    # launch a process for each gene (up to the number of threads specified)
	p.apply_async(run_codeml, [id])


p.close()
p.join() # Wait for all child processes to close.

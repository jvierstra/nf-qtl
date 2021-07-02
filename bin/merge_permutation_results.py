#!/usr/bin/env python3
# Author: Francois Aguet (Modified by: Jeff Vierstra)
import argparse
import os
import numpy as np
import subprocess
import gzip
import contextlib
from datetime import datetime
import tempfile

script_dir = os.path.abspath(os.path.dirname(__file__))

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


def merge_chunks(chunk_list_file, header, output_dir, prefix):
    """Merge FastQTL output chunks and add header"""
    with tempfile.NamedTemporaryFile(mode='w+b', dir=output_dir) as header_file:
        with open(header_file.name, mode='wt') as file:
            file.write('\t'.join(header)+'\n')
            file.flush()
        cmd = "xargs -I $ sh -c 'zcat $ | tail -n +2 ' < {} | cat {} - | gzip -c -1 > {}".format(
            chunk_list_file, header_file.name, os.path.join(output_dir, prefix+'.txt.gz'))
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Merge tensorQTL permutation results and calculate significance')
    parser.add_argument('chunk_list', help='List of chunks')
    parser.add_argument('prefix', help='Prefix for output file name')
    parser.add_argument('--fdr', default=0.05, type=np.double)
    parser.add_argument('--qvalue_lambda', default=None, help='lambda parameter for pi0est in qvalue.')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    args = parser.parse_args()
    fastqtl_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Merging chunks', flush=True)
   
    # header = [
    #     'gene_id', 'num_var', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df',
    #     'variant_id', 'tss_distance', 'ma_samples', 'ma_count', 'maf', 'ref_factor',
    #     'pval_nominal', 'slope', 'slope_se', 'pval_perm', 'pval_beta'
    # ]

    with open(args.chunk_list) as chunk_list_file:
        with gzip.open(chunk_list_file.readline().strip(), 'rt') as chunk_file:
            header = chunk_file.readline().strip().split('\t')
    
    with cd(args.output_dir):
        merge_chunks(args.chunk_list, header, args.output_dir, args.prefix)
     
        print('Calculating q-values', flush=True)
        cmd = os.path.join(script_dir, 'compute_significance.R')+' '+args.prefix+'.txt.gz '+str(args.fdr)+' '+args.prefix+'.phenotypes.txt.gz'
        if args.qvalue_lambda is not None:
            cmd += ' --lambda '+args.qvalue_lambda
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        os.remove(args.prefix+'.txt.gz')

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Done', flush=True)
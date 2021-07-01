#!/usr/bin/env python3

import sys

import torch

import numpy as np
import pandas as pd

import tensorqtl
from tensorqtl import cis

from pandas_plink import read_plink

plink_prefix_path = sys.argv[1]
phenotype_bed_path = sys.argv[2]
covariates_path = sys.argv[3]
include_chrs = [sys.argv[4]]

prefix = 'all'
output_dir = '.'


window = 1000000
maf_threshold = 0.01
permuations = 10000


print("Loading covariates...")

covariates_df = pd.read_table(covariates_path, header=0)
covariates_df = covariates_df.drop(columns=covariates_df.columns[0]).rename(columns = {"IID": "id"}).set_index("id")

print("Done!")

# load phenotypes

print("Loading phenotype data...")

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_path)
phenotype_sample_ids = phenotype_df.columns.tolist()

select_phenotypes = phenotype_pos_df['chr'].isin(include_chrs).values
phenotype_df = phenotype_df[select_phenotypes]
phenotype_pos_df = phenotype_pos_df[select_phenotypes]

print("Done!")

# plink genotypes

print("Scanning genotype data...")

bim, fam, bed = read_plink(plink_prefix_path, verbose=True)

bed[np.isnan(bed)] = -1  # convert missing (NaN) to -1 for int8
bed = bed.astype(np.int8, copy=False)

# Filter chromosomes
m = bim['chrom'].isin(include_chrs).values
bed = bed[m,:]
bim = bim[m]
bim.reset_index(drop=True, inplace=True)
bim['i'] = bim.index

# Filter samples
fam_sample_ids = fam['iid'].tolist()

select_samples = list(set(phenotype_sample_ids) & set(fam_sample_ids))
ix = [fam_sample_ids.index(i) for i in select_samples]

fam = fam.loc[ix]
bed = bed[:, ix]

covariates_df = covariates_df.loc[select_samples]

print("Found {} samples in both genotype and phenotype files!".format(len(ix)))

print("Loading genotype data into memory...")

genotype_df = pd.DataFrame(bed.compute(), index=bim['snp'], columns=fam['iid'])
variant_df = bim.set_index('snp')[['chrom', 'pos']]

print("Done!")

print("Running emperical permutations...")


cis_empirical = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                     covariates_df=covariates_df, window=window, 
                     maf_threshold=maf_threshold, nperm=permuations)

out_file = os.path.join(output_dir, prefix +'.cis_qtl.' + include_chrs[0] + '.txt.gz')

print("Writing results to {}...".format(out_file))

cis_empirical.to_csv(out_file, header=True, index=True, sep="\t")

print("Done!")

print("Running nominal permutations...")

cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, 
                covariates_df=covariates_df, window=window, 
                maf_threshold=maf_threshold, write_stats=True,
                output_dir=output_dir)

print("Done!")

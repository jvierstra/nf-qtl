#!/usr/bin/env python3

import sys
import os

import torch

import numpy as np
import pandas as pd

import tensorqtl
from tensorqtl import cis

from pandas_plink import read_plink

def unpack_region(s):
    (chrom, coords) = s.split(":")
    (start, end) = coords.split("-")
    return chrom, int(start), int(end)

plink_prefix_path = sys.argv[1]
phenotype_bed_path = sys.argv[2]
covariates_path = sys.argv[3]
region = sys.argv[4]

(region_chrom, region_start, region_end) = unpack_region(region)

prefix = region
output_dir = '.'


window = 500000
maf_threshold = 0.01
permuations = 10000


print("Loading covariates...")

covariates_df = pd.read_table(covariates_path, header=0)
covariates_df = covariates_df.drop(columns=covariates_df.columns[0]).rename(columns = {"IID": "id"}).set_index("id")

print("Done!")

# load phenotypes

print("Loading phenotype data...")

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_path)
#phenotype_sample_ids = phenotype_df.columns.tolist()

phenotype_sample_ids = phenotype_df.columns
phenotype_indiv_ids = phenotype_df.columns.str.split('>').str[0]
phenotype_celltypes =  phenotype_df.columns.str.split('>').str[1]
phenotype_sample_df = pd.DataFrame({'sample_id': phenotype_sample_ids, 'indiv_id': phenotype_indiv_ids, 'celltype': phenotype_celltypes})

select_phenotypes = (phenotype_pos_df['chr']==region_chrom) & (phenotype_pos_df['tss']>region_start) & (phenotype_pos_df['tss']<=region_end)
phenotype_df = phenotype_df[select_phenotypes]
phenotype_pos_df = phenotype_pos_df[select_phenotypes]

print("Done!")

# plink genotypes

print("Scanning genotype data...")

bim, fam, bed = read_plink(plink_prefix_path, verbose=True)

bed = 2 - bed # flip alleles!!!!
bed[np.isnan(bed)] = -1  # convert missing (NaN) to -1 for int8
bed = bed.astype(np.int8, copy=False)

# Filter chromosomes
#m = bim['chrom'].isin(include_chrs).values

windowed_region_start = (region_start - window)
windowed_region_end = (region_end + window)

m = (bim['chrom'].isin([region_chrom]).values) & (bim['pos']>windowed_region_start) & (bim['pos']<=windowed_region_end)

bed = bed[m,:]
bim = bim[m]
bim.reset_index(drop=True, inplace=True)
bim['i'] = bim.index

# Filter samples
fam_sample_ids = fam['iid'].tolist()

#select_samples_id = sorted(list(set(phenotype_sample_ids) & set(fam_sample_ids)))
#
#covariates_df = covariates_df.loc[select_samples]
#phenotype_df = phenotype_df[select_samples]
#
#ix = [fam_sample_ids.index(i) for i in select_samples]
#fam = fam.loc[ix]
#bed = bed[:, ix]

i = phenotype_sample_df['indiv_id'].isin(fam_sample_ids)


select_phenotypes_id = phenotype_sample_df[i]['sample_id']
select_indivs_id = phenotype_sample_df[i]['indiv_id']

interaction_df = phenotype_sample_df[i].pivot_table(columns = 'celltype', index='sample_id', aggfunc=lambda x: 1, fill_value=0)
interaction_df = interaction_df.loc[phenotype_df.columns]

if len(interaction_df.columns) < 2:
    interaction_df = None

assert len(select_phenotypes_id) == len(select_indivs_id)

covariates_df = covariates_df.loc[select_indivs_id]
phenotype_df = phenotype_df[select_phenotypes_id]

ix = [fam_sample_ids.index(i) for i in select_indivs_id]
fam = fam.loc[ix]
bed = bed[:, ix]

print("Found {} samples in both genotype and phenotype files!".format(len(ix)))

print("Loading genotype data into memory...")

# genotype_df = pd.DataFrame(bed.compute(), index=bim['snp'], columns=fam['iid'])
covariates_df.index = select_phenotypes_id
genotype_df = pd.DataFrame(bed.compute(), index=bim['snp'], columns=select_phenotypes_id)
variant_df = bim.set_index('snp')[['chrom', 'pos']]

print("Done!")

print("Running emperical permutations...")

cis_empirical = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                     covariates_df=covariates_df, window=window, 
                     maf_threshold=maf_threshold, nperm=permuations)

out_file = os.path.join(output_dir, prefix +'.cis_qtl.' + region_chrom + '.txt.gz')

print("Writing results to {}...".format(out_file))

cis_empirical.to_csv(out_file, header=True, index=True, sep="\t")

print("Done!")

print("Running nominal permutations...")

cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix,
                covariates_df=covariates_df, interaction_df=interaction_df, window=window, 
                maf_threshold=maf_threshold, write_stats=True,
                output_dir=output_dir)

print("Done!")

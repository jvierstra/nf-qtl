#!/usr/bin/env python3

import sys

import pandas as pd
import numpy as np
import scipy as sp

from sklearn.preprocessing import QuantileTransformer


regions_annotation_filepath = sys.argv[1] #'/net/seq/data/projects/regulotyping/dnase/by_celltype_donor/h.CD3+/index/masterlist_DHSs_h.CD3+_nonovl_core_chunkIDs.gc-content.K76mappable.bed'
raw_tag_counts_filepath =  sys.argv[2] #'/net/seq/data/projects/regulotyping/dnase/by_celltype_donor/h.CD3+/index/tag_counts/matrix_tagcounts.txt.gz'
samples_filepath = sys.argv[3] #'/tmp/samples.txt'
exclude_chrs = ['chrX', 'chrY', 'chrM']

raw_tag_counts = pd.read_csv(raw_tag_counts_filepath, header=0, index_col=0, delimiter="\t")

regions_annotations = pd.read_csv(regions_annotation_filepath, header=None, delimiter="\t")
regions_annotations.columns = ["chr", "start", "end", "n_bases", "n_gc", "percent_gc", "n_mappable", "region_id", "mid"]
regions_annotations.set_index("region_id", inplace=True)

with open(samples_filepath) as f:
	samples = [line.rstrip() for line in f]

# Check whether files match
assert raw_tag_counts.index.equals(regions_annotations.index), "Counts and annotation files do not match!"

# Remove unwanted chromosomes
filtered_regions = regions_annotations.index[~regions_annotations["chr"].isin(exclude_chrs)]
raw_tag_counts = raw_tag_counts.loc[filtered_regions]
regions_annotations = regions_annotations.loc[filtered_regions]

# Get relevant samples (e.g., to match VCF file)
sample_cols = raw_tag_counts.columns.intersection(samples)
raw_tag_counts = raw_tag_counts[sample_cols]

# Normalize by total counts & # of mappable base in element

normalized_tag_counts = raw_tag_counts.div(raw_tag_counts.sum(axis=0)) * 1e6
normalized_tag_counts = normalized_tag_counts.div(regions_annotations["n_mappable"].values, axis=0)

# GC content normalization

normalized_tag_counts["gc_bin"] = pd.qcut(regions_annotations["percent_gc"], q=50, duplicates='drop').values

gc_bins = normalized_tag_counts.groupby('gc_bin')
gc_medians = normalized_tag_counts.groupby('gc_bin').transform(np.median)

normalized_tag_counts = normalized_tag_counts[gc_medians.columns] - gc_medians

# Mean and variance scaling

row_means = normalized_tag_counts.mean(axis=1)
row_sigmas = normalized_tag_counts.std(axis=1)
normalized_tag_counts = normalized_tag_counts.subtract(row_means, axis=0).div(row_sigmas, axis=0)


# Quantile normalization

qt = QuantileTransformer(n_quantiles=1000, random_state=0, output_distribution='normal')

normalized_tag_counts = pd.DataFrame(qt.fit_transform(normalized_tag_counts), 
                                     index=normalized_tag_counts.index, 
                                     columns=normalized_tag_counts.columns)

# Drop rows that are all NAs
normalized_tag_counts.dropna(axis='rows', inplace=True)

df = regions_annotations.reset_index()[["chr", "mid", "end", "region_id"]].join(normalized_tag_counts, on="region_id", how="right")
df.sort_values(by = ["chr", "mid"], inplace=True)

# Rename columns for compatibility with TensorQTL
df.rename(columns={"chr": "#chr", "mid": "start", "region_id": "phenotype_id"}, inplace=True)
df["end"] = df["start"] + 1


df.to_csv(sys.stdout, header=True, index=False, sep="\t", float_format="%0.4f")


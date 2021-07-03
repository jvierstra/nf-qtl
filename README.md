# nf-qtl
Nextflow QTL mapping pipeline using tensorQTL

## Overview

## Requirements

**Python:**
- tensorQTL
- PyTorch

**R:** 
- qvalue

## Usage

```
module load jdk/1.8.0_92
module load nextflow/0.29.1

# use an up to date version of nextflow instead of our creaky old decrepit one
export NXF_VER=21.04.1
```

## Pipeline overview

<detail><summary>Step 1: Normalize count matrix</summary>
<p></p>
</detail>

<detail><summary>Step 2: Convert VCF to plink & generate population genetic structure covariates</summary>
<p></p>
</detail>

<detail><summary>Step 3: Chunk genome</summary>
<p></p>
</detail>

<detail><summary>Step 4: Run QTL analysis for each chunk</summary>
<p></p>
</detail>

<detail><summary>Step 5: Compute FDR for permutation pass</summary>
<p></p>
</detail>

<detail><summary>Step 6: Apply FDR to a variant-phenotye pairs</summary>
<p></p>
</detail>
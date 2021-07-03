#!/usr/bin/env nextflow

params.vcf_file='/net/seq/data/projects/regulotyping/genotypes/genotype_panel/imputed_genotypes/chroms1-22.phaseI+II.annotated.ancestral.vcf.gz'
params.count_matrix_file='/net/seq/data/projects/regulotyping/dnase/by_celltype_donor/h.CD3+/index/tag_counts/matrix_tagcounts.txt.gz'
params.regions_file='/net/seq/data/projects/regulotyping/dnase/by_celltype_donor/h.CD3+/index/masterlist_DHSs_h.CD3+_nonovl_core_chunkIDs.bed'

params.genome='/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts'
params.chunksize=25000000

params.outdir='output'

//DO NOT EDIT BELOW

genome_fasta_file="$params.genome"  + ".fa"
genome_chrom_sizes_file="$params.genome"  + ".chrom_sizes"
genome_mappable_file="$params.genome" + ".K76.mappable_only.bed"


process normalize_count_matrix {
	executor 'local'

	publishDir params.outdir, mode: 'copy'

	input:
	file 'matrix_counts.txt.gz' from file(params.count_matrix_file)
	file 'genotypes.vcf.gz' from file(params.vcf_file)
	file 'regions.bed' from file(params.regions_file)
	
	file 'genome.fa' from file(genome_fasta_file)
	file 'genome.fa.fai' from file("${genome_fasta_file}.fai")
	file 'mappable.bed' from file(genome_mappable_file)

	output:
	file 'regions_annotated.bed.gz*'
	set file('matrix_counts.norm.bed.gz'), file('matrix_counts.norm.bed.gz.tbi') into NORMED_COUNTS_FILES, NORMED_COUNTS_REGIONS

	script:
	"""
	faidx -i nucleotide -b regions.bed genome.fa \
		| awk -v OFS="\\t" 'NR>1 { total =\$4+\$5+\$6+\$7+\$8; cg=\$6+\$7; print \$1, \$2-1, \$3,total, cg, cg/total;  }' \
		| bedmap --delim "\t" --echo --bases-uniq - mappable.bed \
		| paste - <(cut -f4,9 regions.bed) \
	| bgzip -c > regions_annotated.bed.gz

	bcftools query -l genotypes.vcf.gz > samples.txt

	normalize_counts.py regions_annotated.bed.gz matrix_counts.txt.gz samples.txt \
	| bgzip -c > matrix_counts.norm.bed.gz
	tabix -p bed matrix_counts.norm.bed.gz
	"""
}

// Select bi-allelic SNVs and make plink files; PCA on genotypes for covariates
process make_plink {
	executor 'local'

	publishDir params.outdir + '/plink' , mode: 'copy'

	input:
	file 'genotypes.vcf.gz' from file(params.vcf_file)

	output:
	file "plink.*" into PLINK_FILES

	script:
	"""
	plink2 --make-bed \
    	--output-chr chrM \
    	--vcf genotypes.vcf.gz \
    	--snps-only \
    	--out plink

    plink2 \
    	--bfile plink \
    	--pca \
    	--out plink
	"""
}

// Chunk genome up only look at regions with in the phenotype matrix
process create_genome_chunks {
	executor 'local'

	input:
	file(count_matrix) from NORMED_COUNTS_REGIONS.map{ it[0] }
	file 'chrom_sizes.txt' from file(genome_chrom_sizes_file)
	val chunksize from params.chunksize

	output:
	stdout into GENOME_CHUNKS

	script:
	"""
	zcat ${count_matrix} | cut -f1-3 | sort-bed - > regions.bed

	cat chrom_sizes.txt \
  	| awk -v step=${chunksize} -v OFS="\\t" \
		'{ \
			for(i=step; i<=\$2; i+=step) { \
				print \$1, i-step+1, i; \
			} \
			print \$1, i-step+1, \$2; \
		}' \
	| sort-bed - > chunks.bed

	bedops -e 1 chunks.bed regions.bed | awk -v OFS="\\t" '{ print \$1":"\$2"-"\$3; }'
	"""
} 

process qtl_by_region {
	tag "${region}"

	//publishDir params.outdir + '/qtl', mode: 'copy'

	label 'gpu'

	input: 
	set file(count_matrix), file(count_matrix_index) from NORMED_COUNTS_FILES
	each region from GENOME_CHUNKS.flatMap{ it.split() }
	file '*' from PLINK_FILES.collect()

	output:
	file "*.txt.gz" into QTL_EMPIRICAL
	file "*.parquet" into QTL_PAIRS_NOMINAL

	script:
	"""
	source /home/jvierstra/.local/share/venv/gpu-python3.8/bin/activate

	qtl.py plink ${count_matrix} plink.eigenvec ${region}
	"""
}

process merge_permutations {
	executor 'local'
	
	publishDir params.outdir + '/qtl', mode: 'copy'

	module "R/4.0.5:python/3.8.10"

	input:
	file '*' from QTL_EMPIRICAL.collect()

	output:
	file 'all.phenotypes.txt.gz' into QTL_EMPIRICAL_SIGNIF

	script:
	"""
	find \$PWD -name "chr*.txt.gz" > filelist.txt

	merge_permutation_results.py filelist.txt all
	"""
}

QTL_PAIRS_NOMINAL
	.map{ it -> 
		def names = (it.name.split(":")) 
		tuple(names[0], it)
	}
	.groupTuple(by: 0)
	.set{ QTL_PAIRS_NOMINAL_BY_CHR }

process filter_nominal_pairs {
	tag "${chr}"

	publishDir params.outdir + '/qtl', mode: 'copy'

	executor 'local'
	module 'python/3.6.4'

	input:
	set val(chr), file('*') from QTL_PAIRS_NOMINAL_BY_CHR 
	file phenotypes_file from QTL_EMPIRICAL_SIGNIF

	output:
	file "${chr}.signifpairs.txt.gz" into QTL_PAIRS_SIGNIF_BY_CHR

	script:
	"""
	ls *.parquet > filelist.txt

	merge_nominal_results.py --fdr 0.05 ${phenotypes_file} filelist.txt ${chr}.signifpairs.txt.gz
	"""
}

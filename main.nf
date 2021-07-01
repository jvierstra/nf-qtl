#!/usr/bin/env nextflow

params.outdir='output'

params.vcf_filepath = '/net/seq/data/projects/regulotyping/genotypes/genotype_panel/imputed_genotypes/chroms1-22.phaseI+II.annotated.ancestral.vcf.gz'
params.regions_filepath = '/net/seq/data/projects/regulotyping/dnase/by_celltype_donor/h.CD3+/index/masterlist_DHSs_h.CD3+_nonovl_core_chunkIDs.bed'
params.count_matrix = '/net/seq/data/projects/regulotyping/dnase/by_celltype_donor/h.CD3+/index/tag_counts/matrix_tagcounts.txt.gz'
params.genome='/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts'


fasta_reference_filepath="$params.genome" + ".fa"
mappable="$params.genome" + ".K76.mappable_only.bed"


process annotate_regions {
	executor 'local'

	publishDir params.outdir, mode: 'copy'

	input:
	file 'regions.bed' from file(params.regions_filepath)
	file 'genome.fa' from file(fasta_reference_filepath)
	file 'genome.fa.fai' from file("${fasta_reference_filepath}.fai")
	file 'mappable.bed' from file(mappable)

	file 'matrix_counts.txt.gz' from file(params.count_matrix)
	file 'genotypes.vcf.gz' from file(params.vcf_filepath)

	output:
	file 'regions_annotated.bed.gz*'
	set file('matrix_counts.norm.bed.gz'), file('matrix_counts.norm.bed.gz.tbi') into NORMED_COUNTS

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

NORMED_COUNTS.into{NORMED_COUNTS_FILES;NORMED_COUNTS_REGIONS}

process collect_chrs {
	executor 'local'

	input:
	file(count_matrix) from NORMED_COUNTS_REGIONS.map{ it[0] }.first()

	output:
	stdout into CHRS

	script:
	"""
	zcat ${count_matrix} | cut -f1 | tail -n +2 | sort | uniq
	"""
}

// Select bi-allelic SNVs and make plink files
process make_plink {
	executor 'local'

	publishDir params.outdir + '/plink' , mode: 'copy'

	input:
	file 'genotypes.vcf.gz' from file(params.vcf_filepath)

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

process qtl_by_chr {
	tag "${chr}"

	publishDir params.outdir + '/qtl', mode: 'copy'

	label 'gpu'

	input: 
	set file(count_matrix), file(count_matrix_index) from NORMED_COUNTS_FILES
	each chr from CHRS.flatMap{ it.split() }
	file '*' from PLINK_FILES.collect()

	output:
	file "*.txt.gz" into QTL_PAIRS_NOMINAL
	file "*.parquet" into QTL_EMPIRICAL

	script:
	"""
	source /home/jvierstra/.local/share/venv/gpu/bin/activate

	qtl.py plink ${count_matrix} plink.eigenvec ${chr}
	"""
}
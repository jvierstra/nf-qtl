#!/usr/bin/env Rscript


# Author: Francois Aguet (modified by Jeff Vierstra for tensorQTL)
# Description: This is a post-processing function for tensorQTL. It calculates 
# q-values (Storey) and p-value thresholds for all genes in the permutation results file.

suppressMessages(library(qvalue))
suppressMessages(library(tools))
suppressMessages(library(argparser))

# parse inputs
p <- arg_parser("Annotates TensorQTL permutation output and runs qvalue")
p <- add_argument(p, "tensorqtl_output", help="tensorQTL permutation ouptut")
p <- add_argument(p, "fdr", type="numeric", help="FDR cutoff")
p <- add_argument(p, "outfile", help="Output file")
p <- add_argument(p, "--lambda", type="numeric", help="Use lambda to calculate q-values", default=NULL)
args <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))

cat("Processing tensorQTL output (", args$tensorqtl_output, ") with FDR=", args$fdr, "\n", sep="")
tensorqtl.df <- read.table(args$tensorqtl_output, header=TRUE, stringsAsFactors=FALSE)
stopifnot("pval_beta" %in% colnames(tensorqtl.df))

# remove genes w/o variants
nanrows <- is.na(tensorqtl.df[, 'pval_beta'])
tensorqtl.df <- tensorqtl.df[!nanrows, ]
cat("  * Number of phenotypes tested: ", nrow(tensorqtl.df), " (excluding ",
    sum(nanrows), " phenotypes w/o variants)\n", sep="")
cat("  * Correlation between Beta-approximated and empirical p-values: ",
    round(cor(tensorqtl.df[, 'pval_perm'], tensorqtl.df[, 'pval_beta']), 4), "\n", sep="")

# calculate q-values
if (is.null(args$lambda) || is.na(args$lambda)) {
    Q <- qvalue(tensorqtl.df[, 'pval_beta'])
} else {
    cat("  * Calculating q-values with lambda = ", args$lambda, "\n", sep="")
    Q <- qvalue(tensorqtl.df[, 'pval_beta'], lambda=args$lambda)
}

tensorqtl.df$qval <- signif(Q$qvalues, 6)
cat("  * Proportion of significant phenotypes (1-pi0): " , round((1 - Q$pi0), 2), "\n", sep="")
cat("  * ePhenotypes @ FDR ", args$fdr, ":   ", sum(tensorqtl.df[, 'qval']<args$fdr), "\n", sep="")

# determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
ub <- sort(tensorqtl.df[tensorqtl.df$qval > args$fdr, 'pval_beta'])[1]  # smallest p-value above FDR
lb <- -sort(-tensorqtl.df[tensorqtl.df$qval <= args$fdr, 'pval_beta'])[1]  # largest p-value below FDR
pthreshold <- (lb+ub)/2
cat("  * min p-value threshold @ FDR ", args$fdr, ": ", pthreshold, "\n", sep="")
tensorqtl.df[, 'pval_nominal_threshold'] <- signif(qbeta(pthreshold,
    tensorqtl.df[, 'beta_shape1'], tensorqtl.df[, 'beta_shape2'], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)

write.table(tensorqtl.df, gzfile(args$outfile), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
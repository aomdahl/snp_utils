# snp_utils
6/7/2021: Tools for analyzing list of SNPs (mapping to genes, checking for enrichment, rank by given metric, etc.)
* *SNPsToBed.R* : takes a list of SNPs (either RSID or chr:pos) and converts into a bed file. If provide a list of values associated with the SNPs (scores, GWAS sum stats, etc.) can also return a list with just a certain percentage of the SNPs. Initially optimized for use with `gwas_decomp_ldsc`, but moved here for wider use.
* *SNPsToTSV.R* : takes a list of SNPs and (either RSID or chr:pos) and converts to a TSV file for querying on (OpenCravat) [https://run.opencravat.org/submit/nocache/index.html]

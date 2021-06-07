# snp_utils
6/7/2021: Tools for analyzing list of SNPs (mapping to genes, checking for enrichment, rank by given metric, etc.)
* **SNPsToBed.R** : takes a list of SNPs (either RSID or chr:pos) and converts into a bed file. If provide a list of values associated with the SNPs (scores, GWAS sum stats, etc.) can also return a list with just a certain percentage of the SNPs. Initially optimized for use with `gwas_decomp_ldsc`, but moved here for wider use.
 ```USAGE:
 Rscript src/SNPsToBed.R --mapping gwas_extracts/seed2_thresh0.9_h2-0.1_vars1e-5/seed2_thresh0.9_h2-0.1_vars1e-5.pruned_rsids.txt --snp_list results/seed2_thresh0.9_h2-0.1_vars1e-5/factorization/snp.ids.txt --outdir results/seed2_thresh0.9_h2-0.1_vars1e-5/gene_set_enrichment/full_list.bed
#yields issue with order
  bedtools closest -a results/seed2_thresh0.9_h2-0.1_vars1e-5/gene_set_enrichment/full_list.bed -b ../../reference_data/protein_co
 ```
* *SNPsToTSV.R* : takes a list of SNPs and (either RSID or chr:pos) and converts to a TSV file for querying on [OpenCravat](https://run.opencravat.org/submit/nocache/index.html)
* *SNPsToGenes.sh* : takes a BED file of SNPs and matches them to their closest genes. (hg19)

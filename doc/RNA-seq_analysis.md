## RNAseq

This data is focused on the processing the RNA-Seq data. The workflow for this process is outlined in a 'snakemake' file.

`snakemake -s expression_workflow.snakefile -p --configfile config.yml`

This workflow includes the following steps:

1. Make an HTSeq matrix from the individual raw htseq files from each tissue:
 	  * `/home/hirschc1/shared/projects/fractionation/cache/rnaseq/B73_mRNA_htseq-count-matrix.txt`
  	* `/home/hirschc1/shared/projects/fractionation/cache/rnaseq/PH207_mRNA_htseq-count-matrix.txt`  
2. Append the B73 and PH207 matrices together into a single file:
  	* `/home/hirschc1/shared/projects/fractionation/cache/rnaseq/B73_PH207_CDS_appended-count-matrix.txt`
3. Correct for differences in transcript length between B73 and PH207 for these raw counts:
  	* `/home/hirschc1/shared/projects/fractionation/cache/rnaseq/B73_PH207_mRNA_htseq-count-matrix-len-corrected.txt`
4. Import length-corrected matrix to DESeq to call differentially expressed genes and export normalized counts
 	* `B73_PH207_mRNA_htseq-count-matrix-len-corrected_htseq-normalized.txt` - normalized DE values before separating by tissue
5. Create file for manipulation in R of the average DESeq counts for each de gene
	* `/home/hirschc1/shared/projects/fractionation/cache/rnaseq/b73_ph207_CDS_htseq-count-matrix-len-corrected_htseq-normalized-avr.txt` - this the norm counts table

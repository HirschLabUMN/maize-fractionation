CODE=config["dir"]["munge"]
CACHE=config["dir"]["cache"] + "/rnaseq/"
GENOS = ['b73', 'ph207']

rule all:
	input:
		expand('{dir}/{gen}_CDS_htseq-count-matrix.txt', gen=GENOS,dir=CACHE),
		CACHE+GENOS[0] + "_" + GENOS[1] + "_CDS_appended-count-matrix.txt",
		CACHE+GENOS[0] + "_" + GENOS[1] + "_CDS_htseq-count-matrix-len-corrected.txt",
		CACHE+GENOS[0] + "_" + GENOS[1] + "_CDS_htseq-count-matrix-len-corrected_htseq-normalized.txt"
		#CACHE+'normalized_rpkm_count_matrix.txt'

# Create a basic count matrix from htseq counts by looking into each tophat directory
rule matrix:
	input:
		workdir = '/scratch.global/broha006/projects/frac/reads/rnaseq/{gen}'
	output:
		"{dir}/{gen}_CDS_htseq-count-matrix.txt"
	shell:
		'perl {CODE}/make_htseq_matrix.pl --workdir {input.workdir} --matrix_out {output}'
# Combine b73 and ph207 matrices together
rule combine:
	input:
		in1 = CACHE+GENOS[0] + '_CDS_htseq-count-matrix.txt',
		in2 = CACHE+GENOS[1] + '_CDS_htseq-count-matrix.txt'
	output:
		CACHE+GENOS[0] + "_" + GENOS[1] + "_CDS_appended-count-matrix.txt"
	run:
		"""
		touch {output}
		perl {CODE}/combine_b73andph207_count_matrices.pl -i {input.in1} -p {input.in2} -o {output}
		perl -i -pe 's/\.v1\.1//g' {output}
		"""
# Correct for differing transcript lengths between B73 and PH207 syntelogs and reformat matrix
rule len_correct:
	params:
		masterFile = '/home/hirschc1/shared/projects/fractionation/cache/coge/b73-ph207-synteny-bt2corrected-nah-coord.txt'
	input:
		rules.combine.output
	output:
		CACHE+GENOS[0] + "_" + GENOS[1] + "_CDS_htseq-count-matrix-len-corrected.txt"
	shell:
		"""
		touch {output}
		perl {CODE}/normalize_count_matrix_by_transcript_len.pl -c {params.masterFile} -i {input} -o {output}
		"""
# Run DESeq2 to find differentially expressed genes between syntelogs
rule deseq:
	message: "Running DESeq2 now"
	input: rules.len_correct.output
	output: CACHE+GENOS[0] + "_" + GENOS[1] + "_CDS_htseq-count-matrix-len-corrected_htseq-normalized.txt"
	shell:
		"""
		Rscript {CODE}/run_DESeq2.Rscript {input} {output}
		perl -i -pe 's/"//g' {output}
		"""

#perl make-de-gene-matrix.pl -d /home/hirschc1/shared/projects/fractionation/cache/rnaseq/ -o out.test

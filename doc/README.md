# Fractionation Project  

*Alex Brohammer*  
*University of Minnesota, Twin-Cities*

## Setup directory structure

`mkdir -p ./{src,graphs,reports,lib,logs/std_out,logs/std_error,data,doc,diagnostics,munge,config,cache}`

## Download CoGe files

[Sorghum vs. B73](https://genomevolution.org/r/o6rf)  
[Sorghum vs. PH207](https://genomevolution.org/r/o6rg)  
[Rice vs. B73](https://genomevolution.org/r/o6rd)  
[Rice vs. PH207](https://genomevolution.org/r/o6re)  

## Configuration Set-up

### Create config file

These are hardcoded paths to default directories, genomes, gffs, etc.

Add these things to [yaml](https://github.com/darvid/trine/wiki/YAML-Primer) file for easy retrieval within many of the scripts.

See config `src/config.yml` for this and `src/cluster.yaml`.

## Parse tandem duplicate output analysis

See `Compare_SynMap_to_Orthofinder.py` script for filtering process.

- False tandem duplicates:  `B73_False_Tandem_Consensus.txt` and `PH207_False_Tandem_Consensus.txt`  
	- These files have some meta-data on how they were found:
		- BOTH: SynMap and Orthofinder both find the same maize-sorghum ortholog  
		- BLAST: SynMap finds as tandem but one has higher sequence similarity (ds)  
		- BLAST_Leftmost: Synmap finds as tandem, but about equal ds so left-most is chosen by default  
		- SynMap_Chosen: Synmap finds as tandem and orthoFinder finds different maize-sorghum ortholog so choose SynMap  
		- Other_Orthologue: False tandem with different orthologs - additional relationships from orthoFinder.   
		- Orthofinder_Only: OrthoFinder finds assignment but SynMap does not.  
- True tandem duplicates: `B73_True_Tandems_WithID.txt` and `PH207_True_Tandems_WithID.txt`  

`B73_All_Tandem_Fractionation_Representative_wStatus.txt` and `PH207_All_Tandem_Fractionation_Representative_wStatus.txt` contain both False tandem duplicates and true tandem duplicates. Since some of these tandem duplicates were found by orthoFinder and not SynMap they aren't all necessarily collinear, so the next script will determine this.  

B73:
`perl /home/hirschc1/shared/projects/fractionation/munge/tandem-passfail-assign.pl -i /home/hirschc1/shared/projects/fractionation/data/tandem_assignments/B73_All_Tandem_Fractionation_Representative_wStatus.txt -m /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_no-duplicates.txt -q b73 -s sb -o /home/hirschc1/shared/projects/fractionation/data/tandem_assignments/B73_All_Tandem_Fractionation_Representative_wStatus-passfail.txt`  

PH207:
`perl /home/hirschc1/shared/projects/fractionation/munge/prefilter_toms_file.pl -i /home/hirschc1/shared/projects/fractionation/data/tandem_assignments/PH207_All_Tandem_Fractionation_Representative_wStatus.txt -m /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_no-duplicates.txt -q ph207 -s sb -o /home/hirschc1/shared/projects/fractionation/data/tandem_assignments/PH207_All_Tandem_Fractionation_Representative_wStatus-passfail.txt`  

Combine B73 and PH207 parsed tandems:  
`cat *passfail.txt > B73_PH207_All_Tandem_Fractionation_Representative_wStatus-passfail.txt`  

### Set up storable database

This will store the following information in a hash stored on disk and retrieved at various points in the workflow. I am using the Perl module called [Storable](https://perldoc.perl.org/Storable.html) for this. See [here](http://www.dcs.ed.ac.uk/home/perl5/Storable.html) for a tutorial.

This is basically a mini-database that contains the following:

* gene (key) 
	* cds_length - length of representative CDS
	* chr - chromosome
	* cognate - the homolog from the same subgenome in the other genotype
	* homeolog - the homolog from the opposite subgenome
	* length - the gene length
	* order - the order of the gene along the chromosome
	* ortholog - the sorghum/rice ortholog
	* reciprocal - The homolog from the opposite genotype
	* rep_cds - the representative cds Zm00001d..._T003
	* rep_trans - the reprentative transcript Zm00001d..._T002
	* start - gene start coordinates
	* stop  - gene stop coordinates
	* subgenome - maize1 or maize2
	* synteny - 1 if syntenic


Create the storable dbs using this script:

```
perl /home/hirschc1/shared/projects/fractionation/munge/create_storable_dbs.pl -i b73 -m Gramene
perl /home/hirschc1/shared/projects/fractionation/munge/create_storable_dbs.pl -i ph207 -m Phytozome
perl /home/hirschc1/shared/projects/fractionation/munge/create_storable_dbs.pl -i sb -m Phytozome
perl /home/hirschc1/shared/projects/fractionation/munge/create_storable_dbs.pl -i os -m Phytozome
```

*Note:* These may need to be updated with additional information throughout the workflow. 

### Parse SynMap Results

Ok, now we are ready to parse the SynMap results. Most of the analysis is semi-automated in what I refer to as 'modules'. Each module is a `snakemake` file. [Snakemake](https://snakemake.readthedocs.io/en/latest/) is a python implementation of `GNU Makefiles` with some additional extensions.

Run this first module which just parses the raw SynMap output and the second which builds the preliminary list of syntenic assignments.

*Note:* make sure directory structure and file locations match what is in the modules. 

```
snakemake -s Snakefile -F -p --configfile config.yml
snakemake -s Snakefile2.snakemake --configfile config.yml -p -F
```

These script will do the following:

* Parse the Quota Align output to a more readable form
* Report some basic stats
* Match up homeologs in both B73 and PH207 outputs
* Parse the raw SynMap tandem duplicate file
* Make a master file with information
* Run tblastx to get rid of non-tandem duplicates in preliminary master file

Ok, that's it! We now have our first master file. Now we can verify and try to find more syntelogs using BLAST.

### Blast Runs

Thankfully, there is another makefile that largely automates this process.

But we have some work to do setting-up a few things before-hand:

1. Make sure you have blast database created
	`perl make_fasta_index.pl -i /home/maize/shared/databases/blast/Zea_mays/B73/Zea_mays.AGPv4.cds.all.fa -o /home/hirschc1/shared/projects/fractionation/data/assests/B73CdsFaIndex`  
  `perl make_fasta_index.pl -i /home/maize/shared/databases/blast/Zea_mays/PH207/ZmaysPH207_443_v1.1.cds.fa -o /home/hirschc1/shared/projects/fractionation/data/assests/PH207CdsFaIndex`  
2. Create fasta index for bioperl
3. Load required modules: `module load ncbi_blast+/2.6.0`

### Pairwise Blast

**NOTE:** Make sure to update the config file with the latest version of the master file!!

Once that is complete we are ready to perform the pairwise blast runs. This will be done for every maize1 and maize2 gene within a genotype and between every maize(1,2)/maize(1,2) gene across genotypes.

`snakemake -s pairwise-blast.snakemake -p --configfile config.yml --cluster-config cluster.yaml --cluster "qsub -l {cluster.l} -m {cluster.m} -N {cluster.N} -r {cluster.r} -V" --jobs 5`

### Genomewide Blast

This step consists of two blast runs:

1. Blast to genome space but require mapping around gene region. This eliminates issues with local annotation.
2. Blast to recover NA genes. This mitigates issues with a gene being missed by `SynMap` or the entire gene not being annotated.

`snakemake -s genomewide-blast.snakemake -p -F --configfile config.yml --cluster-config cluster.yaml --cluster "qsub -l {cluster.l} -m {cluster.m} -N {cluster.N} -r {cluster.r} -V" --jobs 6`

### Post-blast parsing

*Note:* There are a lot of scripts for this post-processing, not all of which are always necessary. There is some description of the process below and also see the wiki page for additional summaries of individual scripts. 

Now we need to combine these files (*remember to delete headers*)

`cat B73vB73_blastNA_output-scored.txt PH207vPH207_blastNA_output-scored.txt PH207vB73_blastNA_output-scored.txt B73vPH207_blastNA_output-scored.txt >> B73_PH207_blastNA-output.txt`

`cat B73vB73_blastn-parsed-scored.txt PH207vPH207_blastn-parsed-scored.txt PH207vB73_blastn-parsed-scored.txt B73vPH207_blastn-parsed-scored.txt >> B73_PH207_blastn-output.txt`

`cat B73vB73_blast-parsed-scored.txt PH207vPH207_blast-parsed-scored.txt B73vPH207_blast-parsed-scored.txt >> B73_PH207_pairwise-blast-output.txt`

```  
perl /home/hirschc1/shared/projects/fractionation/munge/add_pairwise_status_to_master_file.pl \
  --master /scratch.global/broha006/projects/frac/masterFile_no-duplicates.txt \
  --blast2seq /scratch.global/broha006/projects/frac/B73_PH207_pairwise-blast-output.txt \
  --blastn /scratch.global/broha006/projects/frac/B73_PH207_blastn-output.txt \
  --blastNA /scratch.global/broha006/projects/frac/B73_PH207_blastNA-output.txt \
  --output /scratch.global/broha006/projects/frac/masterFile_parsed.txt  
```  

To deal with possible cases of identifying the wrong maize gene when two different maize genes overlap, I insert the mapping coordinates into the master file when they fail the pairwise blast step but pass in the blastn step. Here is the full rundown:

1.) Remove duplicates using tblastx (again)  
2.) After incorporating new genes, do a new round of NA blast  
3.) Score this alignment and incorporate new genes  
4.) Do a third round of duplicate removal to get rid of any leftover duplicates  

`snakemake -p -F -s post-blast.snakemake --configfile config.yml`
`snakemake -p -F -s post-blast2.snakemake --configfile config.yml`
`snakemake -p -F -s post-blast3.snakemake --configfile config.yml`

The third module will just run another iteration of blast for the sake of completeness, but doesn't necessarily need to be included for all uses.

Ok, this last module should always be ran. It will do one final targeted blast to a window where the gene would be expected to be located and also fill out differentially fractionated genes that are likely due to inconsistent annotation (specifically gene fusions where there are two gene models that represent a single gene in PH207).

`snakemake -p -F -s post-blast4.snakemake --configfile config.yml`

## Master-file QC

At this point, you have a finalized master file and should just be checking for anything that doesn't make sense.

Add a note to discriminate between probable gene fusions:  

```
perl /home/hirschc1/shared/projects/fractionation/munge/add-notes-to-master-file.pl -i probable_gene_fusions.txt > annotated-fusions.txt
awk '{print $0"\t","."}' masterFile_no-duplicates_r5.txt > b73-ph207-synteny.txt
cat annotated-fusions.txt >> b73-ph207-synteny.txt
```

You can also add on CoGe links: 

See `generate-gevo-links.pl`

## SQLite Database Setup  

I use a SQLite database to hold various datasets. This is especially useful for the expression analysis. 

#### Set up SQLite database  

You can get an idea of the various tables and schema from thhis script: 

`perl /home/hirschc1/shared/projects/fractionation/munge/make_pav_db.pl`  

This script is important for generating the assignments table in the database. 

*Note:* - Make sure the input file is latest version of syntenic gene list. 

```
perl /home/hirschc1/shared/projects/fractionation/munge/get_masterFile_summary.pl \
	-i /home/hirschc1/shared/projects/fractionation/cache/coge/b73-ph207-synteny.txt \
	-o /home/hirschc1/shared/projects/fractionation/data/sqlite/inputs/masterFile_summary.txt
```  

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

Further expression analysis was done using `make_expression_plots.Rscript`. This script is probably most easily ran interactively. 

## Non-allelic Homolog Workflow

There is another module for this analysis 

`snakemake -s bt2-realign-tracking.snakemake --configfile config.yml -p`  

We can actually use this method to fill in a handful of genes that we previously missed, so the above two scripts can be ran twice: 1. just looking for any missed syntenic genes, and 2. looking for non-allelic homologs. The `incorporate_bt2_fixes.pl` script at the end of this code block will insert the missed genes into the master file.  

See `analyze_nonSynDFGs.Rscript` for further downstream analysis of non-allelic homologs. 

## PAV analysis

To get coverage of CDS sequences run this on filtered BAM files: 

`qsub processDNAseqForExonCov.sh`  

To re-do analysis use:  
`pav-plot.Rscript`. 


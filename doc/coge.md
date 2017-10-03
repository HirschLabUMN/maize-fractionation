#CoGe  

## Filtering of GFF  

#!/usr/bin/perl -w
use strict;

```perl
while ( my $line = <$in> ) {
  chomp $line;
  if ( $line =~ /^#/ ) {
    print "$line\n";
    next;
  }
  my ( $chr, undef, $feat, undef ) = split '\t', $line;
  if ( $chr !~/(^Mt|^Pt)/ ) {
    if ( $feat =~ /(^gene$|CDS|chromosome|exon|mRNA|contig)/ ) {
      print "$line\n";
    }
  }
}
```

```bash  
perl -i -pe 's/;Name.*;biotype/;biotype/g' Zea_mays.AGPv4.33.gff5  
perl -i -pe 's/CDS://g' Zea_mays.AGPv4.33.gff5
perl -i -pe 's/gene://g' Zea_mays.AGPv4.33.gff5
perl -i -pe 's/transcript://g' Zea_mays.AGPv4.33.gff5
```  


961123 CDS
    10 chromosome
   255 contig
1250081 exon
 39324 gene
     1 #!genebuild-last-updated 2015-12
     1 #!genome-build wareLab AGPv4
     1 #!genome-date 2015-12
     1 #!genome-version AGPv4
     1 ##gff-version 3
133968 mRNA  

## Create on-disk hashes  

```
perl munge/create_storable_dbs.pl B73
perl munge/create_storable_dbs.pl PH207
perl munge/create_storable_dbs.pl Sb
perl munge/create_storable_dbs.pl Os
```

## CoGe re-run links  

iput Zea_mays.AGPv4.33.gff5 coge_data
Link: https://genomevolution.org/r/mhco

These are all in `\data\coge\` and the `README.txt` file.

Everything was ran with the following options:
* Blast Algorithm - LAST
* DAGChainer - Nucleotide Distance
* QuotaAlign Merge / QuotaAlign - 2:1 Zm:Sb / Zm:Os  
* Tandem duplication distance - 15
* C-score - 0.1

## Filter quota alignment output, Get block designations and stats  

`snakemake src/Snakefile`

Makes calls to:  
* `munge/filter_quota_align.pl`
* `enhanced_block_statistics.pl`  
* `filter_condensed_syntelog.pl`  
* `find_sb_duplicates.pl`  
* `find_zm_duplicates.pl`
* `munge/combine_condensed_syntelog_files.pl`  
* `munge/add_OS_genes_to_combined_syntenlog_file.pl`  

All outputs are in `cache`

## Misc.

Number of Block Summary:

|        | Num of Blocks |
|--------|---------------|  
|B73vSb  | 1031 |
|B73vOs  | 1113 |
|PH207vSb| 276 |
|PH207vOs| 375 |  

## Get percent of genome in block space  

This is done by counting the total block span of each block and correcting for overlaps between blocks. So each block will initially be the start of the first CDS and the end of the block will be the last CDS. Everything in between is counted as syntenic.  

perl ../../munge/merge-blocks.pl B73vSb_block_list.txt > B73vSb_merged_blocks.txt
perl ../../munge/merge-blocks.pl B73vOs_block_list.txt > B73vOs_merged_blocks.txt
cat B73vSb_merged_blocks.txt B73vOs_merged_blocks.txt >> B73vSbandOs_merged_blocks.txt
perl ../../munge/merge-blocks.pl B73vSbandOs_merged_blocks.txt > B73vSbandOs_merged_blocks_final.txt

perl ../../munge/merge-blocks.pl PH207vSb_block_list.txt > PH207vSb_merged_blocks.txt
perl ../../munge/merge-blocks.pl PH207vOs_block_list.txt > PH207vOs_merged_blocks.txt
cat PH207vSb_merged_blocks.txt PH207vOs_merged_blocks.txt >> PH207vSbandOs_merged_blocks.txt
perl ../../munge/merge-blocks.pl PH207vSbandOs_merged_blocks.txt > PHvSbandOs_merged_blocks_final.txt

awk 'BEGIN {sum=0} ; {if ($4 == "maize1") {sum+=$3-$2} } END {print sum}'
PH207vSbandOs_merged_blocks_final.txt

awk 'BEGIN {sum=0} ; {if ($4 == "maize2") {sum+=$3-$2} } END {print sum}' PH207vSbandOs_merged_blocks_final.txt

## Get fractionation class summary

Get a list of times that each fractionation class is observed in the master file

*note:* prints to screen

`perl munge/get_master_file_stats.pl -i Zm.vs.Sb_designated_wOsgenes.txt`  

### Make Master File with Pairwise Status Column  

```  
cd /munge
perl add_pairwise_status_to_master_file.pl -i ../cache/coge/Zm.vs.Sb_designated_wOsgenes.txt -b ../data/blast_results/concat_bl2seq_scored.txt -p ../data/blast_results/concat_blsnRealn_scored.txt -n ../data/blast_results/concat_blsnNA_scored.txt -o ../cache/coge/masterFile_wStatus.txt2  
```  

### Filter out poor syntelog assignments  

A value of '9' in the status column indicates that this gene failed in the pairwise blast step and failed to be recovered when realigning to gene space. After contemplating some complicated ways of deciding how 'bad' an assignment had to be to be filtered out of the master file, I decided that a simple solution is best for now. Any status column that contains multiple '9s' is filtered out. The following script will perform this action.

```
perl munge/filter_on_status_column.pl -i /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus.txt2 -o /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filtered.txt
```
From this only 45 lines were filtered from the file.

### Use bowtie2 alignments to fill in missing genes  

More info on this pipeline is in the `QC_master_file.md` file. Genes filled in using this method will be marked 'b'.

```  
perl /home/hirschc1/shared/projects/fractionation/munge/incorporate_bt2_fixes.pl -i /scratch.global/broha006/projects/frac/qc -m /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filtered.txt -o /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC.txt
```

## Incorporate Coordinate fills positions overlapping genes into master file.  

```  
perl /home/hirschc1/shared/projects/fractionation/munge/incorporate_coord_fill_qc.pl -i /home/hirschc1/shared/projects/fractionation/data/coord_fill -m /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC.txt -o /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC2.txt
```  

## Get a list of duplicate genes to filter out of this list  

I think at some point it would be helpful to manually curate these gene lists. There are often clear winners; time is the issue.

```
cut -f 2 /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC2.txt | sort | uniq -d >> duplicates.txt
cut -f 3 /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC2.txt | sort | uniq -d >> duplicates.txt
cut -f 4 /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC2.txt | sort | uniq -d >> duplicates.txt
cut -f 5 /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC2.txt | sort | uniq -d >> duplicates.txt
grep ",Zm" /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC2.txt >> duplicates.txt  
```

## Run QC steps  

See `QC_master_file.md`.

```
perl /home/hirschc1/shared/projects/fractionation/munge/incorporate_bt2_fixes.pl -i /home/hirschc1/shared/projects/fractionation/data/alns/Reseq -m /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filtered.txt -o /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC.txt
```

## Get Summary for SQLite database  

```  
perl /home/hirschc1/shared/projects/fractionation/munge/incorporate_bt2_fixes.pl -i /home/hirschc1/shared/projects/fractionation/data/alns/Reseq -m /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filtered.txt -o /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC.txt  
```  

### Realign genomic short reads

```
qsub /home/hirschc1/shared/projects/fractionation/munge/runBowtie2 -F "B73"
qsub /home/hirschc1/shared/projects/fractionation/munge/runBowtie2 -F "PH207"  
```  

### Make Karyotype Plots  

```  
 Rscript ../../munge/ColoredSorghumChr_Karyotype.Rscript B73vSb_quota_aln_lite.txt PH207vSb_quota_aln_lite.txt  
```  

### Verfiy Fractionation events

This script is used as part of the bowtie2 realignments. These alignments serve two purposes: 1.) For re-calling PAVs, 2.) For QC'ing differentially fractionated genes. This script takes the raw BAM outputs and sorts them, gets per base coverage of exon sequence, and parses the `bedtools` output to a summarized from. See `alns/Reseq/Bt2_DFG_QC/*exon_coverage.txt`.

```  
qsub processDNAseqForExonCov.sh
```

### Get PAV lists  

```
cd /home/hirschc1/shared/projects/fractionation/data/alns/Reseq/Bt2_DFG_QC

perl /home/hirschc1/shared/projects/fractionation/munge/get_coverage_lists_for_pav.pl -i B73vPH207_exon_coverage.txt -p PH207vPH207_exon_coverage_fromCH.txt -o PH207_pav_list.txt

perl /home/hirschc1/shared/projects/fractionation/munge/get_coverage_lists_for_pav.pl -i PH207vB73_exon_coverage.txt -p B73vB73_exon_coverage.txt -o B73_pav_list.txt
```  

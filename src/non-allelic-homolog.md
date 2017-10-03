This document outlines steps taken to examine non-allelic homologs using W22. I want to know whether these non-allelic homologs are real *or* are simply the result of misassembly. W22 will be useful for this because the assembly is not based off of B73 like PH207, CML247, F7, etc. 

Working directory: 
cd /home/hirschc1/shared/projects/fractionation/cache/non-allelic

### Get B73 and PH207 genes 
```  
awk '$1 ~/Zm00001d/' b73-ph207-regions-mappedBack-sorted-merged-overlap-hits-wcov20-uniq-homeologs-filtered.txt | cut -f 1 > b73-genes.txt
awk '$1 ~/Zm00008a/' b73-ph207-regions-mappedBack-sorted-merged-overlap-hits-wcov20-uniq-homeologs-filtered.txt | cut -f 1 > ph207-genes.txt
```  
### Get fasta sequences for each genotype  

```  
module load ncbi_blast+/2.6.0 
blastdbcmd -entry_batch b73-genes.txt -db /home/maize/shared/databases/blast/Zea_mays/B73/Zea_mays.AGPv4.cds.all.fa -out b73-non-allelic-cds.fasta
blastdbcmd -entry_batch ph207-genes.txt -db /home/maize/shared/databases/blast/Zea_mays/PH207/ZmaysPH207_443_v1.1.cds.fa -out ph207-non-allelic-cds.fasta
```  

### Use last to align B73 and Ph207 genes against W22
`qsub run-last.sh`  

### Parse the last output to get top hits  
```
perl get-last-hits.pl -i b73-alns.maf -o b73-alns-parsed.txt
perl get-last-hits.pl -i ph207-alns.maf -o ph207-alns-parsed.txt
```  

### Summarize whether these genes hit to the unexpected or expected chromosome in W22  

*note:* the input file is created by using awk to get only the genotype you want from the `b73-ph207-regions-mappedBack-sorted-merged-overlap-hits-wcov20-uniq-homeologs-filtered.txt` file. 

```  
perl get-aln-results.pl -i b73_genes.2.txt -r b73-alns-parsed.txt -o b73-parsed-results.txt
perl get-aln-results.pl -i ph207_genes.2.txt -r ph207-alns-parsed.txt -o ph207-parsed-results.txt
```  


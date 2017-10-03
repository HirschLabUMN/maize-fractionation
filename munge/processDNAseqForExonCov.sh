#!/bin/bash
#PBS -l walltime=14:00:00
#PBS -l nodes=1:ppn=9
#PBS -l mem=6gb
#PBS -V
#PBS -r n


INPUT="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
GENOTYPE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
#cd /panfs/roc/scratch/tk_maize/genome_content_var/Merged/Against_B73v4
cd /panfs/roc/scratch/tk_maize/genome_content_var/Merged/Against_PH207v1
module load samtools
module load bedtools
echo "..sorting ${INPUT}"
#samtools view -b -@ 9 -q 20 -o ${GENOTYPE}_B73v4_Merged.unique.bam ${INPUT}
#bedtools coverage -d -sorted -d -a /home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.33.rep-CDS-wCoords.gff3 -b ${GENOTYPE}_B73v4_Merged.unique.bam > ${GENOTYPE}_B73v4_Merged_unique_cds-cov.txt
#perl /home/hirschc1/shared/projects/fractionation/munge/get_exon_coverage.pl -i ${GENOTYPE}_B73v4_Merged_unique_cds-cov.txt -g b73 -d 0 > ${GENOTYPE}_B73v4_Merged_unique_cds-cov-summary.txt
#rm ${GENOTYPE}_B73v4_Merged_unique_cds-cov.txt

samtools view -b -@ 9 -q 20 -o ${GENOTYPE}_PH207v1_Merged.unique.bam ${INPUT}
#bedtools coverage -d -sorted -d -a /home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.gene_exons_cds-only.gff3 -b ${GENOTYPE}_PH207v1_Merged.unique2.bam > ${GENOTYPE}_PH207v1_Merged_unique_cds-cov.txt
bedtools coverage -d -sorted -g /home/hirschc1/shared/projects/fractionation/munge/ph207.genome -a /home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.gene_exons_cds-only.gff4 -b ${GENOTYPE}_PH207v1_Merged.unique.bam > ${GENOTYPE}_PH207v1_Merged_unique_cds-cov2.txt
perl /home/hirschc1/shared/projects/fractionation/munge/get_exon_coverage.pl -i ${GENOTYPE}_PH207v1_Merged_unique_cds-cov2.txt -g ph207 -d 0 > ${GENOTYPE}_PH207v1_Merged_unique_cds-cov-summary.txt
rm ${GENOTYPE}_PH207v1_Merged_unique_cds-cov.txt



#echo "${GENOTYPE}, ${INPUT}"
#
#samtools sort -O bam -@ 9 -T /panfs/roc/scratch/brohamm1/sambamba_tmp ${INPUT} | samtools view -b -@ 9 -q 20 - > "${GENOTYPE}_B73v4_Merged.sorted.unique.bam"
#echo "getting coverage"
#bedtools coverage -sorted -d -a /home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.33.exons-only.gff4 -b W606S_B73v4_Merged.sorted.unique.bam > W606S_B73v4_Merged.sorted.unique.exon-coverage.txt
#bedtools coverage -sorted -d -a /home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.33.exons-only.gff4 -b ${GENOTYPE}_B73v4_Merged.sorted.unique.bam > ${GENOTPE}_B73v4_Merged.sorted.unique.exon-coverage.txt


##################
##  B73 vs B73  ##
##################
#

# cd /scratch.global/broha006/projects/frac/reads/reseq/B73
#samtools sort B73_B73_bt2_unique.bam -o B73_B73_bt2_unique.sorted.bam
#samtools sort B73_PH207_bt2_unique.bam -o B73_PH207_bt2_unique.sorted.bam
#samtools sort PH207_PH207_bt2_unique.bam PH207_PH207_bt2_unique.sorted.bam

# echo "getting BvB coverage"
# bedtools coverage -sorted -d -a /home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.33.exons.sorted.bed -b B73_B73_bt2_unique.sorted.bam > B73_B73_bt2_unique.sorted.coverage.txt
# echo "parsing BvB coverage"
# perl /home/hirschc1/shared/projects/fractionation/munge/get_exon_coverage.pl -i B73_B73_bt2_unique.sorted.coverage.txt -g B73 > B73vsB73_exon_coverage.txt

# echo "getting BvP coverage"
# bedtools coverage -sorted -d -a /home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.gene_exons.exons.sorted.bed -b B73_PH207_bt2_unique.sorted.bam > B73_PH207_bt2_unique.sorted.coverage.txt
# echo "parsing BvsP coverage"
# perl /home/hirschc1/shared/projects/fractionation/munge/get_exon_coverage.pl -i B73_PH207_bt2_unique.sorted.coverage.txt -g PH207 > B73vsPH207_exon_coverage.txt

# cd /scratch.global/broha006/projects/frac/reads/reseq/PH207
# bedtools coverage -sorted -d -a /home/maize/shared/databases/genomes/Zea_mays/B73/Zea_mays.AGPv4.33.exons.sorted.bed -b PH207_B73_bt2_unique.sorted.bam > PH207_B73_bt2_unique.sorted.coverage.txt
#
# bedtools coverage -sorted -d -a /home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.gene_exons.exons.sorted.bed -b PH207_PH207_bt2_unique.sorted.bam > PH207_PH207_bt2_unique.sorted.coverage.txt

#######################
##   PH207 vs PH207  ##
#######################

# Note: Don't really need to do this; b/c you can use older data from CH

# cd /scratch.global/broha006/projects/frac/bt2realn/PH207vPH207
#
# echo "sorting file"
# sambamba sort -t 15 -p --tmpdir /scratch.global/broha006/projects/frac/bt2realn/PH207vPH207/tmp -F "mapping_quality >= 20" pvp.bt2.bam -o pvp.bt2.sorted.bam
#echo "getting coverage"
# bedtools coverage -sorted -d -a /scratch.global/broha006/projects/frac/bt2realn/PH207_exon_sorted.bed -b pvp.bt2.sorted.bam > pvp_exon_coverage.txt
#cd /panfs/roc/scratch/broha006/projects/frac/reads/reseq/PH207
#bedtools coverage -sorted -d -a /home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.gene_exons.exons.sorted.bed -b pvp.bt2.no6.unique.sorted2.bam > pvp_exon_coverage.txt
#cd /panfs/roc/scratch/broha006/projects/frac/reads/reseq/B73
#bedtools coverage -sorted -d -a /home/maize/shared/databases/genomes/Zea_mays/PH207/ZmaysPH207_443_v1.1.gene_exons.exons.sorted.bed -b bvp.bt2.unique.sorted.bam > bvp_exon_coverage.txt2
# echo "parsing coverage"
# perl /home/hirschc1/shared/projects/fractionation/munge/get_exon_coverage.pl -i /scratch.global/broha006/projects/frac/reads/reseq/PH207/PH207_PH207_bt2_unique.sorted.coverage.txt -g PH207 > /scratch.global/broha006/projects/frac/reads/reseq/PH207/PH207vsPH207_exon_coverage.txt
#
# ##################
# ## PH207 vs B73 ##
# ##################
#
# cd /scratch.global/broha006/projects/frac/bt2realn/PH207vB73
#
# echo "starting sort..."
# sambamba sort -t 10 -p --tmpdir /scratch.global/broha006/projects/frac/bt2realn/PH207vB73/tmp -F "mapping_quality >= 20" pvb.bt2.bam
# sambamba sort -t 10 -p --tmpdir /scratch.global/broha006/projects/frac/bt2realn/PH207vB73/tmp -F "mapping_quality >= 20" pvb.bt2.2.bam
# echo "finished sorting, merging..."
# sambamba merge pvb.bt2.merged.sorted.bam pvb.bt2.sorted.bam pvb.bt2.2.sorted.bam
# echo "finished merging, getting coverage..."
# bedtools coverage -sorted -d -a /scratch.global/broha006/projects/frac/bt2realn/B73_exon_sorted.gff -b pvb.bt2.merged.sorted.bam > pvb_exon_coverage.txt
# perl /home/hirschc1/shared/projects/fractionation/munge/get_exon_coverage.pl -i /scratch.global/broha006/projects/frac/reads/reseq/PH207/PH207_B73_bt2_unique.sorted.coverage.txt -g B73 > /scratch.global/broha006/projects/frac/reads/reseq/PH207/PH207vsB73_exon_coverage.txt
# echo "job finished"
#
# ##################
# ## B73 vs PH207 ##
# ##################
#
# cd /scratch.global/broha006/projects/frac/bt2realn/B73vPH207
#
# echo "starting job now"
# sambamba sort -t 5 -p -F "mapping_quality >= 20" bvp.bt2.bam
# echo "getting coverage"
# bedtools coverage -sorted -d -a /scratch.global/broha006/projects/frac/bt2realn/PH207_exon_sorted.gff -b bvp.bt2.sorted.bam > bvp_exon_coverage.txt
# perl /home/hirschc1/shared/projects/fractionation/munge/get_exon_coverage.pl -i bvp_exon_coverage.txt -g PH207 > B73vsPH207_exon_coverage.txt
# echo "job finished"

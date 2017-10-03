# Tandem Cluster Analysis
This document describes how tandem duplicate clusters were analyzed, and how
"representative" genes from each cluster were identified. Clusters are
identified on the basis of two proximal genes in one genome matching the same
gene in another genome. There is no explicit search for global sequence
similarity among members of a tandem duplicate cluster in CoGe.

Initial output from CoGe had 2,311 tandem duplicate clusters in B73, containing
5,702 protein coding genes. There were 3,730 tandem duplicate clusters in PH207,
containing 9,343 protein coding genes.

## "True" and "False" Tandem Cluster Classification
Because clusters were not identified by local similarity within one genome, we
apply alignment filters to each cluster, and assign them into "True" or "False"
tandem clusters. "True" clusters have genes that are at least 75% identical at
the nucleotide sequence level, and cover each other in at least 50% in an
alignment. Clusters failiing these critera are labeled as "false" clusters.

Tandem duplicate clusters handled in the following way:

1. The CDS from each tandem duplicate maize gene were translated to amino acids.
2. The translated sequences were multiply aligned with the putative ancestral
   gene, according to the "master file" generated by ABB. Alignments were
   generated with `clustal-omega`.
3. Aligned amino acid sequences were back-translated into nucleotides.
4. Pairwise similarity and coverage were calculated using the `ThetaPi` and
   `Sites` columns from the output of the `compute` program from the `analysis`
   package, written by K. Thornton.
5. Alignments with great than 25% diversity (less than 75% similarity) and with
   gaps representing greater than half the aligned length were classified as 
   "false" clusters. The rest were classified as "true."

A summary of the filtering critera is shown below.

| Filter | B73 Syntenic | B73 Nonsyntenic | PH207 Syntenic | PH207 Nonsyntenic |
|--------|--------------|-----------------|----------------|-------------------|
| True   | 818 (2,403)  | 471 (1,270)     | 900 (2,939)    | 679 (2,102)       |
| False  | 873 (1,746)  | 149 (298)       | 1,677 (3,354)  | 474 (948)         |
| Total  | 1,691 (4,149)| 620 (1,568)     | 2,577 (6,293)  | 1,153 (3,050)     |

*: Numbers in parentheses are the total number of protein coding genes within
the clusters.

Scripts that perform the cluster classification are

- `Tandem_CDS_Align.py`: translation, alignment, backtranslation
- `Classify_Tandem.R`: Write files of "true" and "false" tandem clusters

## "False" Cluster Handling
### SynMap Representative
Visual inspection of "false" tandem clusters showed that they often contain
what appear to be spurious associations between a maize gene and a *Sorghum*
gene. These are caused by local similarity, despite global dissimilarity. There
are also a few cases where two small maize genes align to non-overlapping
regions of a large *Sorghum* gene. Inspection of the raw output from CoGe showed
that it is possible in some cases to recover single-gene relationships between
*Sorghum* and maize. This is important because putative tandem duplicates are
filtered before assignment of syntenic regions between two genomes in CoGe, and
these "false" tandems may be part of larger syntenic relationships.

To recover genes from "false" tandems, we performed the following procedure.
Note that this only works on genes that have an ancestral gene assignment.

1. Tandem duplicate CDS and ancestral gene CDS were translated to amino acids.
2. Translated sequences were aligned with `clustal-omega` 1.2.1.
3. Aligned amino acid sequences were back-translated to nucleotides.
4. dS was estimated with the `yn00` routine in the PAML package. `yn00`
   estimates divergence with a maximum likelihood approach, accounting for both
   codon frequency and transition-transversion ratio. The parameters in the
   control file were:

    ```
    seqfile = backtranslated_alignment.fasta
    outfile = backtranslated_alignment.fasta.yn00_out
    verbose = 1
    icode = 0
    weighting = 0
    commonf3x4 = 0
    ```

5. Maize genes with dS greater than 2.5 from the ancestral gene were removed.
   This is the maximum observed value for synonymous divergence between maize
   and *Sorghum* in syntenic regions, as reported by Schanble 2011 in PNAS.
6. The maize gene with the lowest dS from the ancestral gene were returned. If
   the interval of (est. - SE, est. + SE) overlapped with other genes, multiple
   gene IDs were returned.

A summary of the recovery is shown below.

| Outcome                     | B73   | PH207 |
|-----------------------------|-------|-------|
| One best gene               | 707   | 1,408 |
| Mulitple best genes         | 100   | 103   |
| Unalignable or no ancestral | 215   | 609   |
| Total                       | 1,022 | 2,151 |

Scripts that perform the "false" tandem handling are:

- `Tandem_CDS_Align.py`: tanslation, alignment, backtranslation.
- `Tandem_dNdS.py`: estimate divergence with `yn00`
- `Summarize_Tandem_Divergence.py`: parse the `rst` output from `yn00` for each
  cluster, and compile into a nice table.
- `Tandem_Lowest_Anc_Div.py`: compare dS from ancestral gene, determine best
  maize genes from each cluster.

### Additional Orthologues with Orthofinder
Orthofinder v. 1.1.5 was used to identify additional orthologues for genes in
"false" tandem duplicate clusters. FASTA files with representative amino acid
sequences from the following species were downloded. 

| Species                   | Assembly Ver. | Anno. Ver. | Source            |
|---------------------------|---------------|------------|-------------------|
| *Aegilops tauschii*       | ASM34733v1    | ASM34733v1 | Ensembl Plants 34 |
| *Brachypodium distachyon* | 3.1           | 3.1        | Phytozome V12     |
| *Hordeum vulgare*         | ASM32608v1    | ASM32608v1 | Ensembl Plants 34 |
| *Leersia perrieri*        | Lperr_V1.4    | Lperr_V1.4 | Ensembl Plants 34 |
| *Oropetium thomaeum*      | 1.0           | 1.0        | Phytozome V12     |
| *Oryza sativa*            | IRGSP-1.0     | IRGSP-1.0  | Ensembl Plants 34 |
| *Panicum hallii*          | 2.0           | 2.0        | Phytozome V12     |
| *Panicum virgatum*        | 1.1           | 1.1        | Phytozome V12     |
| *Phyllostachys edulis*    | 1.0           | 1.0        | BambooGDB         |
| *Setaria italica*         | 2.2           | 2.2        | Phytozome V12     |
| *Sorghum bicolor*         | 3.1           | 3.1        | Phytozome V12     |
| *Triticum aestivum*       | TGACv1        | TGACv1     | Ensembl Plants 34 |
| *Triticum urartu*         | ASM34745v1    | ASM34745v1 | Ensembl Plants 34 |

Orthofinder was run with the default 'dendroblast' algorithm to infer
orthologous relationships, and default parameters for MCL clustering. Protein
sequences from B73 and PH207 were kept in separate files.

Pairwise orthologus relationships between maize and the two ancestral species,
*Sorghum bicolor* and *Oryza sativa*, were parsed from the Orthofinder output.
Genes from each "false" tandem duplicate cluster were compared to the
representative relationship identifed by divergence from SynMap assignment, and
orthologues inferred from Orthofinder output. Relationships between
non-representative maize genes and *S. bicolor* or *O. sativa* were added to
ABB's "master file."

Scripts that perform Orthofinder clustering and identification of additional
orthologous relationships:

- `Orthofinder.job`: Run Orthofinder on MSI
- `Compare_SynMap_to_Orthofinder.py`: Compare SynMap representative maize gene
and ancestral assignment to Orthofinder output, and generate a consensus
maize-ancestral relationship table.

## Step-by-Step Workflow
Download and install the necessary software. Names and version numbers are given
in the next section. Orthofinder and Analysis were run on the compute cluster
(MSI), and everything else was run on a Mac laptop.

All scripts are available in this GitHub repository.

### True/False Tandem Classification
1. Generate tandem duplicate gene lists for B73 and PH207, using *S. bicolor*
   as an outgroup. Use a C-value of 0.1, and a window size of 15 genes, in
   CoGe. The format is one cluster per line, with genes in each cluster
   separated by a comma.
2. Download the B73 v4, PH207 v1, *S. bicolor* v3.1, and *O. sativa* IRGSP-1.0
   protein sequences and CDS sequences. Store these separately from the
   Orthofinder input sequences.
3. Generate alignments between ancestral loci and maize tandem duplicates for
   each tandem duplicate cluster with `clustal-omega`. The `master_file.txt`
   is generated by ABB, and contains information on paralogue and ancestral
   orthologues.

   ```
   python Tandem_CDS_Align.py B73_CDS.fa Sorghum_CDS.fa Rice_CDS.fa master_file.txt B73_Tandems.csv
   python Tandem_CDS_Align.py PH207_CDS.fa Sorghum_CDS.fa Rice_CDS.fa master_file.txt PH207_Tandems.csv
   ```

   This creates directories, and fills them with multiple sequence alignmetns
   between the maize genes and the SynMap ancestral gene, if applicable.
4. Transfer the alignments with the outgroup sequences to MSI, or a system with
   `compute` from the `analysis` packge available. Summarize the similarity
   with `compute`:

   ```bash
   cd /path/to/B73/alignments
   compute -i 'Anc_*.fasta' -O 1 > B73_Outgroup_Summary.txt
   cd /path/to/PH207/alignments
   compute -i 'Anc_*.fasta' -O 1 > PH207_Outgroup_Summary.txt
   ```

   Note that the single quotes are necessary in the `compute` command, else the
   shell will glob the filenames, and it will throw an error.

5. Parse the outgroup summary table with an R script to determine which tandem
   clusters are "true" and "false." Be sure to edit the proper directories into
   the script.

   ```
   Rscript Classify_Tandem.R
   ```

   This will write `B73_Cluster_Status.txt` and `PH207_Cluster_Status.txt`,
   which contain filtering information for each tandem cluster. If a cluster
   has `Pass` in both the `Similarity` and `Coverage` columns, then it is
   labeled as a "true" tandem duplicate cluster. Otherwise, it is labeled as a
   "false" tandem duplicate cluster.
5. Generate "true" and "false" tandem clusters for each accession.

  ```bash
  awk '$2 == "Pass" && $3 == "Pass" {print $1 "\t" $7}' B73_Cluster_Status.txt > B73_True_Tandems_WithID.txt
  awk '$2 == "Pass" && $3 == "Pass" {print $1 "\t" $7}' PH207_Cluster_Status.txt > PH207_True_Tandems_WithID.txt
  awk '$2 == "Fail" || $3 == "Fail" {print $1 "\t" $7}' B73_Cluster_Status.txt > B73_False_Tandems_WithID.txt
  awk '$2 == "Fail" || $3 == "Fail" {print $1 "\t" $7}' PH207_Cluster_Status.txt > PH207_False_Tandems_WithID.txt
  ```

### SynMap Representative Gene From "False" Tandem Duplicate Clusters
1. Calculate dN/dS between maize genes and ancestral genes with `yn00`, from
   PAML. This script will generate anaylsis directories for each cluster, and
   automatically generate a PAML control file for the analysis.

   ```bash
   python Tandem_dNdS.py /path/to/B73/alignments /path/to/B73/output
   python Tandem_dNdS.py /path/to/PH207/alignments /path/to/PH207/output
   ```

2. Parse the `yn00` output files from each cluster, and print out which
   maize-ancestral pair has the lowest divergence.

   ```bash
   python Tandem_Lowest_Anc_Div.py \
       B73_False_Tandems_withID.txt \
       /path/to/B73/alignments \
       > B73_False_Tandem_Ancestral_Assignment.txt

   python Tandem_Lowest_Anc_Div.py \
       PH207_False_Tandems_withID.txt \
       /path/to/PH207/alignments \
       > PH207_False_Tandem_Ancestral_Assignment.txt
   ```

### Orthofinder Clustering
1. Download the representative protein sequences from the speceis and sources
   given in the table above. Rename the files in the format `Genus_species.fa`
   for easy parsing, and put them into a single directory. These should be on
   the compute cluster (MSI).
2. Because the SynMap assignments were originally done on an earlier version of
   *S. bicolor* genome assembly, the gene IDs do not match between the SynMap
   table and the Orthofinder files. Download the translation table from
   Phytozome. The filename is `Sbicolor_255_v2.1.locus_transcript_name_map.txt`
3. Run Orthofinder on the cluster. This may take several days.

   ```bash
   qusb Orthofinder.job
   ```

4. Copy the `Orthologues/` directory from the Orthofinder output to your local
   machine. The files you need in particulare are
   - `B73__v__Sorghum_bicolor.csv`
   - `B73__v__Oryza_sativa.csv`
   - `PH207__v__Sorghum_bicolor.csv`
   - `PH207__v__Oryza_sativa.csv`
5. Use the following files to determine "best" orthologous relationships between
   maize and the ancestral species. Be sure to edit the paths to the following
   files in the `Tandem_Rep.sh` script. It's just a wrapper for the Python
   script, since so many pieces of data are required.
   - B73 cluster numbers file
   - PH207 cluster numbers file
   - `B73__v__Sorghum_bicolor.csv`
   - `B73__v__Oryza_sativa.csv`
   - `PH207__v__Sorghum_bicolor.csv`
   - `PH207__v__Oryza_sativa.csv`
   - B73 amino acid sequences
   - B73 CDS sequences
   - PH207 amino acid sequences
   - PH207 CDS sequences
   - *S. bicolor* amino acid sequences
   - *S. bicolor* CDS sequences
   - *O. sativa* amino acid sequences
   - *O. sativa* CDS sequences
   - `Rep_Tandem_From_Orthogroup.py` script

  ```bash
  bash Tandem_Rep.sh
  ```

  This will generate `B73_All_Tandem_Fractionation_Representative.txt` and
  `PH207_All_Tandem_Fractionation_Representative.txt`. Note that these are not
  the representative assignments for each gene, despite what the name says. This
  is confusing.

6. Trim down the `B73_All_Tandem_Fractionation_Representative.txt` and
   `PH207_All_Tandem_Fractionation_Representative.txt` files to only have the
   information from "false" tandem duplicate clusters.

   ```bash
   grep -f B73_False_ClusterID.txt B73_All_Tandem_Fractionation_Representative.txt > B73_Orthofinder_False.txt
   grep -f PH207_False_ClusterID.txt PH207_All_Tandem_Fractionation_Representative.txt > PH207_Orthofinder_False.txt
   ```

7. Generate a consensus assignment for false tandem duplicates between SynMap
   and Orthofinder. Use the `Sbicolor_255_v2.1.locus_transcript_name_map.txt`
   file that was downloaded from Phytozome earlier.

   ```bash
   python Compare_SynMap_to_Orthofinder.py \
       Sbicolor_255_v2.1.locus_transcript_name_map.txt \
       B73_Cluster_Numbers.txt \
       B73_False_Tandem_Ancestral_Assignment.txt \
       B73_Orthofinder_False.txt \
       > B73_False_Tandem_Consensus.txt

   python Compare_SynMap_to_Orthofinder.py \
       Sbicolor_255_v2.1.locus_transcript_name_map.txt \
       PH207_Cluster_Numbers.txt \
       PH207_False_Tandem_Ancestral_Assignment.txt \
       PH207_Orthofinder_False.txt \
       > PH207_False_Tandem_Consensus.txt
   ```

   **These files are the consensus assignments for "false" tandem duplicate
   clusters.**


## Software Versions
The versions of all software used are given below.

- Linux (MSI): `Linux 2.6.32-642.15.1.el6.x86_64 #1 SMP Fri Feb 24 14:31:22 UTC 2017 x86_64 x86_64 x86_64 GNU/Linux`
- MacOS: `Darwin 16.5.0 Darwin Kernel Version 16.5.0: Fri Mar  3 16:52:33 PST 2017; root:xnu-3789.51.2~3/RELEASE_X86_64 x86_64`
- Python: `2.7.13`
- Biopython: `1.68`
- R: `R version 3.3.3 (2017-03-06) -- "Another Canoe"` 64 bit
- Clustal-omega: `1.2.1`
- PAML: `4.9d (February 2017)`
- Orthofinder: `1.1.5` (commit `203a861ff24b0f63aebba6c513829dff3bb3f15e`)
- Analysis (`module load molpopgen` in MSI): From April 3, 2002
#!/usr/bin/perl -w
use strict;
use DBI;
my $dbh = DBI->connect( "dbi:SQLite:dbname=/home/hirschc1/shared/projects/fractionation/data/sqlite/all_info.db", "", "",
        { RaiseError => 1 });

# .import /home/hirschc1/shared/projects/fractionation/data/b73-ph207-regions-mappedBack-sorted-merged-overlap-hits-wcov20-uniq.txt nah
$dbh->do("CREATE TABLE IF NOT EXISTS nah
  (	gene TEXT PRIMARY_KEY NOT NULL,
    homolog TEXT,
    cov INTEGER
)");

# .import /home/hirschc1/shared/projects/fractionation/data/sqlite/inputs/exon-coverage.txt reseq_cov
$dbh->do("CREATE TABLE IF NOT EXISTS reseq_cov
  (	gene TEXT PRIMARY_KEY NOT NULL,
    aa_cov INTEGER,
    ab_cov INTEGER
)");

# .import /home/hirschc1/shared/projects/fractionation/data/sqlite/inputs/ph207-gff-meta.txt gff
# .import /home/hirschc1/shared/projects/fractionation/data/sqlite/inputs/b73-gff-meta.txt gff
$dbh->do("CREATE TABLE IF NOT EXISTS gff
(	chr INTEGER,
	start INTEGER,
	stop INTEGER,
  gene TEXT PRIMARY_KEY NOT NULL,
	chrorder INTEGER,
	len INTEGER,
  exon_count INTEGER
)");

# .import /home/hirschc1/shared/projects/fractionation/data/sqlite/inputs/curated-gene-lists.txt curated
$dbh->do("CREATE TABLE IF NOT EXISTS curated
(	gene TEXT PRIMARY_KEY NOT NULL
)");

# .import /home/hirschc1/shared/projects/fractionation/data/sqlite/inputs/masterFile_summary.txt assignments
$dbh->do("CREATE TABLE IF NOT EXISTS assignments
( gene TEXT PRIMARY_KEY NOT NULL,
  subgenome TEXT,
  duplicate TEXT,
  cognate TEXT,
  reciprocal TEXT,
  genotype TEXT
)");

# .import /home/hirschc1/shared/projects/fractionation/data/sqlite/inputs/classical_gene_list_v4.txt classical_genes
$dbh->do("CREATE TABLE IF NOT EXISTS classical_genes
(	gene TEXT PRIMARY_KEY NOT NULL
)");

# .import /home/hirschc1/shared/projects/fractionation/data/sqlite/inputs/maizeGDB-curated-genelist_v4.txt maizegdb_curated
$dbh->do("CREATE TABLE IF NOT EXISTS maizegdb_curated
(	gene TEXT PRIMARY_KEY NOT NULL
)");

# .import /home/hirschc1/shared/projects/fractionation/data/sqlite/inputs/wallace-gwashits-v4.txt gwas_hits
$dbh->do("CREATE TABLE IF NOT EXISTS gwas_hits
(	gene TEXT PRIMARY_KEY NOT NULL,
	trait TEXT,
	gwas_hit TEXT
)");

# .import /home/hirschc1/shared/projects/fractionation/data/sqlite/inputs/walley_hub_genes_v4.txt hub_genes
$dbh->do("CREATE TABLE IF NOT EXISTS hub_genes
(	gene TEXT PRIMARY_KEY NOT NULL,
	edge_count INTEGER
)");

#.import /home/hirschc1/shared/projects/fractionation/cache/rnaseq/b73_ph207_CDS_htseq-count-matrix-len-corrected_htseq-normalized-avr.txt normCounts
$dbh->do("CREATE TABLE IF NOT EXISTS normCounts
(	gene TEXT PRIMARY_KEY NOT NULL,
  b73Bd TEXT,
  b73Cp TEXT,
  b73Gk TEXT,
  b73Rt TEXT,
  b73Sd TEXT,
  b73St TEXT,
	ph207Bd TEXT,
	ph207Cp TEXT,
	ph207Gk TEXT,
	ph207Rt TEXT,
	ph207Sd TEXT,
	ph207St TEXT
)");

$dbh->disconnect;
exit;
# .import /home/hirschc1/shared/projects/fractionation/cache/rnaseq/B73_PH207_de-gene-matrix-for-R.txt2 deCounts
# $dbh->do("CREATE TABLE IF NOT EXISTS deCounts
# (	gene TEXT PRIMARY_KEY NOT NULL,
#   bd TEXT,
#   cp TEXT,
#   gk TEXT,
#   rt TEXT,
#   sd TEXT,
#   st TEXT,
# 	id TEXT,
#   syntelog TEXT
# )");
# .import /home/hirschc1/shared/projects/fractionation/cache/rnaseq/B73_cds_RPKM_matrix_avr.txt b73RPKM
# $dbh->do("CREATE TABLE IF NOT EXISTS b73RPKM
# (	gene TEXT PRIMARY_KEY NOT NULL,
#   bd TEXT,
#   cp TEXT,
#   gk TEXT,
#   rt TEXT,
#   sd TEXT,
#   st TEXT
# )");
# .import /home/hirschc1/shared/projects/fractionation/cache/rnaseq/PH207_cds_RPKM_matrix_avr.txt ph207RPKM
# $dbh->do("CREATE TABLE IF NOT EXISTS ph207RPKM
# (	gene TEXT PRIMARY_KEY NOT NULL,
#   bd TEXT,
#   cp TEXT,
#   gk TEXT,
#   rt TEXT,
#   sd TEXT,
#   st TEXT
# )");

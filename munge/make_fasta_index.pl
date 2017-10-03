#!/usr/bin/perl
use strict;
use Getopt::Std;
use Data::Printer;
use Bio::Index::Fasta;
##ABOUT: This script will make a fasta index for multifasta files for bioperl blast runs
##RUNLIKE: perl make_fasta_index.pl -i /home/hirschc1/shared/projects/fractionation/data/assests/Zea_mays.AGPv4.cds.all.fa -o B73index

our ($opt_i, $opt_o);
getopts("i:o:");

#Bio::Seq object
my $index_name = $opt_o;

my $inx = Bio::Index::Fasta->new( -filename => $index_name,
                                  -write_flag => 1);
$inx->make_index($opt_i);
exit 0;

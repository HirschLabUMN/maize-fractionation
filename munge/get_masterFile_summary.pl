#!/usr/bin/perl -w
##ABOUT: This script is for the SQLite database. It produces a gene-wise summary of the master file for querying
##RUNLIKE: perl /home/hirschc1/shared/projects/fractionation/munge/get_masterFile_summary.pl -i /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC3.txt -o /home/hirschc1/shared/projects/fractionation/data/sqlite/inputs/masterFile_summary.txt
use strict;
use Getopt::Std;

our($opt_i,$opt_o);
getopts("i:o:");

open ( my $infile, '<', $opt_i ) or die;
open ( my $outfile, '>', $opt_o ) or die;

while ( my $line = <$infile> ) {
  chomp $line;

  my ( $anchr, $m1, $m2, $m3, $m4 ) = split '\t', $line;
  my @genes = ($m1,$m2,$m3,$m4);
  my ($subgenome, $duplicate, $cognate, $reciprocal);
  if ( $m1 !~ /NA/ ) {
    $subgenome = 'maize1';
    $duplicate = $m2;
    $cognate = $m3;
    $reciprocal = $m4;
    print $outfile "$m1\t$subgenome\t$duplicate\t$cognate\t$reciprocal\tB73\n";
  }

  if ( $m2 !~ /NA/ ) {
    $subgenome = 'maize2';
    $duplicate = $m1;
    $cognate = $m4;
    $reciprocal = $m3;
    print $outfile "$m2\t$subgenome\t$duplicate\t$cognate\t$reciprocal\tB73\n";
  }

  if ( $m3 !~ /NA/ ) {
    $subgenome = 'maize1';
    $duplicate = $m4;
    $cognate = $m1;
    $reciprocal = $m2;
    print $outfile "$m3\t$subgenome\t$duplicate\t$cognate\t$reciprocal\tPH207\n";
  }

  if ( $m4 !~ /NA/ ) {
    $subgenome = 'maize2';
    $duplicate = $m3;
    $cognate = $m2;
    $reciprocal = $m1;
    print $outfile "$m4\t$subgenome\t$duplicate\t$cognate\t$reciprocal\tPH207\n";
  }

}
close $infile;
close $outfile;

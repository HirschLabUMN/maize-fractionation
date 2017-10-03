#!/usr/bin/perl -w
use strict;
use feature qw/say/;
use Getopt::Std;
use Storable;

##ABOUT: This will add OS genes onto the Zm.vs.Sb_designated file. I.e. sorghum gene was deleted but retained in rice and there are maize orthologs to the rice genes.
##RUNLIKE: perl munge/add_OS_genes_to_combined_syntenlog_file.pl -i Zm.vs.Sb_designated.txt -p Zm.vs.Os_designated.txt -o Os_genes_to_add.txt


our ( $opt_i, $opt_p, $opt_o);
getopts('i:p:o:') or die "Couldn't get all options, $!\n";

open ( my $infile, '<', $opt_i ) or die "error opening $opt_i\n";
open ( my $infile2, '<', $opt_p ) or die "error opening $opt_p\n";
open ( my $outfile, '>', $opt_o ) or die "error opening $opt_o\n";

my %hash;
while ( my $line = <$infile> ) {
  chomp $line;

  my (undef, $m1, $m2, $m3, $m4) = split '\t', $line;

  # Capture genes in an array
  my @m1genes = split ',', $m1; # Need these b/c of the weird duplicates
  my @m2genes = split ',', $m2;
  my @m3genes = split ',', $m3;
  my @m4genes = split ',', $m4;

  # combined arrays into a single array
  my @allgenes = ( @m1genes, @m2genes, @m3genes, @m4genes );

  # Store contents of concatenated array into a hash
  foreach my $gene ( @allgenes ) {
    $hash{$gene} = '' unless $gene eq 'NA';
  }
}
close $infile;

## OPEN OS DESIGNATED FILE AND PRINT LINES WHERE NONE OF THE ZM GENES ARE IN THE LOOKUP HASH ##

while ( my $line = <$infile2> ) {
  chomp $line;

  my ($anchr, $m1, $m2, $m3, $m4) = split '\t', $line;

  my @m1genes = split ',', $m1;
  my @m2genes = split ',', $m2;
  my @m3genes = split ',', $m3;
  my @m4genes = split ',', $m4;

  # Combine all genes into single array
  my @genes = ( @m1genes, @m2genes, @m3genes, @m4genes );

  my $print; #flag that tells whether to print line or not

  #Loop through array, if none of genes are seen flag is set to 1 and line will be printed.
  foreach my $gene ( @genes ) {
    if ( $gene eq 'NA' ) {
      next;
    }
    elsif ( exists($hash{$gene}) ) {
      $print = 0;
      last;
    }
    else {
      $print = 1;
    }
  }

  #if the flag is set to 1, none of the genes were found in lookup hash so print
  if ( $print == 1 ) {
    say $outfile $line;
  }
}

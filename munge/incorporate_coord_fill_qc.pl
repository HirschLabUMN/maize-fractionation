#!/usr/bin/perl -w
use strict;
use Storable;
use Library;
use Getopt::Std;
use Data::Printer;
##ABOUT: This script will filter the output of bt2 realignments, which were used to explore false positive
# DFGs. Upon running bedtools intersect, I found that many reads map uniquely to a B73 gene model which was
# missed by SynMap. I want to incorporate these 'missed' genes that are on the right chromosome into the master File.
##RUNLIKE: perl /home/hirschc1/shared/projects/fractionation/munge/incorporate_coord_fill_qc.pl -i /home/hirschc1/shared/projects/fractionation/data/coord_fill -m /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC.txt -o /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC2.txt

my $usage = "\n\n$0 -i <INPUT> -o <OUTPUT>\n\n";

our ($opt_i, $opt_o, $opt_m, $opt_h);
getopts("i:o:m:h") or die "couldn't get options";

open (my $mastrfile, '<', $opt_m) or die "Cannot open $opt_m $!\n\n";
open (my $outfile, '>', $opt_o) or die "Cannot open $opt_o $!\n\n";

my $href = &getAssignedGenes($opt_i);

while ( my $line = <$mastrfile> ) {
  chomp $line;

  my ( $sb, $bm1, $bm2, $pm1, $pm2 ) = split '\t', $line;

  if ( $bm1 =~ /:/ ) {
    if ( exists($href->{$bm1}) ) {
      $bm1 = $href->{$bm1};
    }
  }
  if ( $bm2 =~ /:/ ) {
    if ( exists($href->{$bm2}) ) {
      $bm2 = $href->{$bm2};
    }
  }
  if ( $pm1 =~ /:/ ) {
    $pm1 =~ s/chr0//g; $pm1 =~ s/chr//g;
    if ( exists($href->{$pm1}) ) {
      $pm1 = $href->{$pm1};
    }
  }
  if ( $pm2 =~ /:/ ) {
    $pm2 =~ s/chr0//g; $pm2 =~ s/chr//g;
    if ( exists($href->{$pm2}) ) {
      $pm2 = $href->{$pm2};
    }
  }
  print $outfile "$sb\t$bm1\t$bm2\t$pm1\t$pm2\n";
}
close $outfile;
close $mastrfile;

sub getAssignedGenes {
  my ( $handle ) = @_;
  my %geneHash;
  my @files = <$handle/*coordinate-overlap.txt>;
	p @files; 
  foreach my $file ( @files ) {
    open (my $fh, '<', $file) or die "error opening $handle\n";
    while ( my $line = <$fh> ) {
      chomp $line;
      my ( $gene, $loc ) = split '\t', $line;
      $loc =~ s/chr0//; $loc =~ s/chr//;
      $geneHash{$loc} = $gene;
    }
    close $fh;
  }
  return(\%geneHash);
}

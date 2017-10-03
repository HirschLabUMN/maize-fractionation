#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Printer;
##ABOUT: This script is used to insert genes that hit to scaffolds into the master file.
##RUNLIKE: perl munge/incorporate_scaffold_hits.pl \
#-i  /home/hirschc1/shared/projects/fractionation/data/blast_results/scaffold_gene_search/scaffold_hits.txt \
#-m /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC4Gaps.txt
#-o  /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC4GapsScaf.txt


our($opt_i, $opt_m);
getopts("i:m:o:");
my $tmp = 'foo.tmp' . $$;
open ( my $oh, ">$tmp");

my $href = GetHits($opt_i);


open (my $mf, '<', $opt_m) or die;
while ( my $line = <$mf> ) {
  chomp $line;

  my ( $sb, $bm1, $bm2, $pm1, $pm2 ) = split '\t', $line;

  if ( exists($href->{$bm1}) ) {
    $pm1 = $href->{$bm1};
  }
  if ( exists($href->{$bm2}) ) {
    $pm2 = $href->{$bm2};
  }
  if ( exists($href->{$pm1}) ) {
    $bm1 = $href->{$pm1};
  }
  if ( exists($href->{$pm2}) ) {
    $bm2 = $href->{$pm2};
  }

  print $oh "$sb\t$bm1\t$bm2\t$pm1\t$pm2\n";
}
close $mf;
close $oh;
rename( $tmp, $opt_m );
exit 0;


sub GetHits {
  my ($handle) = @_;
  my %hash;
  open ( my $infile, '<', $handle ) or die;
  while ( my $line = <$infile> ) {
    chomp $line;
    my @fields = split '\t', $line;
    my ($query, $hitMeta) = @fields[0,3];
    $query =~ s/_.*//;
    my @hitFields = split ',', $hitMeta;
    my $hit = $hitFields[0];
    $hit =~ s/_.*//;
    $hash{$query} = $hit;
  }
  close $handle;
  return(\%hash);
}

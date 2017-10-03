#!/usr/bin/perl -w
use strict;
use Storable;
use Getopt::Std;
use YAML::XS qw/LoadFile/;
use Data::Printer;
##ABOUT: This script is used as part of effort to determine how many of the DFGs
# could not be identified because of assembly gaps. After pulling out list of genes from
# R sqlite queries, I need to get 5kb upstream and downstream of gene. This script is used to get These
# coordinates which will then be pulled from the chromosomal multifasta using bedtools getfasta.
##RUNLIKE: perl get_coordinates_for_gap_analysis.pl -i /home/hirschc1/shared/projects/fractionation/data/sqlite/B73nonscaffoldHits_DFGs.txt -o /home/hirschc1/shared/projects/fractionation/data/gap_analysis/B73_DFG_coordinates.txt -g B73

#perl get_coordinates_for_gap_analysis.pl -i /home/hirschc1/shared/projects/fractionation/data/sqlite/PH207nonscaffoldHits_DFGs.txt -o /home/hirschc1/shared/projects/fractionation/data/gap_analysis/PH207_DFG_coordinates.txt -g PH207

our ($opt_i, $opt_c, $opt_o);
getopts("i:o:c:") or die;

my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";
my ( $query, $subject ) = split 'v', $opt_c;
$query =~ tr/A-Z/a-z/; $subject =~ tr/A-Z/a-z/;

# Get storable database
my $queryHashRef = retrieve($settings->{'db'}->{"$query"}->{'storable'});
my $subjHashRef = retrieve($settings->{'db'}->{"$subject"}->{'storable'});
my $queryOrderRef = &getOrderHash($queryHashRef);
my $subjOrderRef = &getOrderHash($subjHashRef);


open ( my $infile, '<', $opt_i ) or die;
open ( my $outfile, '>', $opt_o ) or die;
while ( my $line = <$infile> ) {
  chomp $line;

  my $coordRef = getCoordinates($line, $queryHashRef, $subjHashRef, $queryOrderRef, $subjOrderRef);
  if ( @{$coordRef} == 0 ) {
    next;
  }
 print $outfile "$coordRef->[0]\t$coordRef->[1]\t$coordRef->[2]\t$line $coordRef->[0]:$coordRef->[1]-$coordRef->[2]\n";
}
close $infile;
close $outfile;
exit 0;

sub getOrderHash {
  my ( $href ) = @_;
  my %ordHash;
  foreach my $key ( keys %$href ) {
    if ( $href->{$key}->{'synteny'} eq '1' ) {
      my $ord = $href->{$key}->{'order'};
      my $chr = $href->{$key}->{'chr'};
      my $st = $href->{$key}->{'start'};
      my $sp = $href->{$key}->{'stop'};
      $ordHash{$chr}{$ord}{'st'} = $st;
      $ordHash{$chr}{$ord}{'sp'} = $sp;
      $ordHash{$chr}{$ord}{'id'} = $key;
    }
  }
  return(\%ordHash);
}

sub getCoordinates {
  my ( $gene, $queryHRef, $subjHRef, $queryORef, $subjORef ) = @_;
  $gene =~ s/_(T|P)\d+//; # get rid of any transcript part of name
  my $order = $queryHRef->{$gene}->{'order'}; # get the gene order
  my $chr = $queryHRef->{$gene}->{'chr'}; # get the original chromosome
  # check to make sure hit is between the nearest up and downsteam syntenic genes.
  ## GET THE NEAREST UPSTREAM SYNTENIC GENE
  my $matchedEnd = 'foo';
  my ( $matchedChrEnd, $matchedChrStart );
  for my $i ( 1..50 ) { # check 50 genes upstream to find nearest match
    my $check = $order - $i; # if gene is order 30 then check 29,28,27,etc.
    if ( exists($queryORef->{$chr}->{$check})  ) { # only syntenic genes exist here
      my $nearestGene = $queryORef->{$chr}->{$check}->{'id'}; # this is the nearest upstream gene model
      if ( $queryHRef->{$nearestGene}->{'cognate'} ne 'NA' ) {
        my $matchedHomeolog = $queryHRef->{$nearestGene}->{'cognate'};
        $matchedChrEnd = $queryHRef->{$nearestGene}->{'chr'};
        if ( $matchedHomeolog !~ /Zm/ ) { next; }
        if ( $subjHRef->{$matchedHomeolog}->{'chr'} ne $chr ) { warn "$matchedHomeolog,$nearestGene, $gene do not match ", $subjHRef->{$matchedHomeolog}->{'chr'},":$chr\n"; next; }
        $matchedEnd = $subjHRef->{$matchedHomeolog}->{'stop'};
        last;
      }
      else {
        next;
      }
    }
    else {
      next;
    }
  }
  my $matchedStart = 'foo';
  for my $i ( 1..50 ) { # check 50 genes upstream to find nearest match
    my $check = $order + $i; # if gene is order 50 then check 49,48,47,etc.
    if ( exists($queryORef->{$chr}->{$check})  ) { # only syntenic genes exist here
      my $nearestGene = $queryORef->{$chr}->{$check}->{'id'}; # this is the nearest upstream gene model
      if ( $queryHRef->{$nearestGene}->{'cognate'} ne 'NA' ) {
        my $matchedHomeolog = $queryHRef->{$nearestGene}->{'cognate'};
        $matchedChrStart = $queryHRef->{$nearestGene}->{'chr'};
        if ( $matchedHomeolog !~ /Zm/ ) { next; }
        $matchedStart = $subjHRef->{$matchedHomeolog}->{'start'};
        last;
      }
      else {
        next;
      }
    }
    else {
      next;
    }
  }
  if ( $matchedEnd ne 'foo' && $matchedStart ne 'foo' && $matchedEnd < $matchedStart && ( $matchedChrEnd eq $matchedChrStart ) ) {

    if ( $gene =~ /Zm00001d/ ) { $chr = "chr0" . $chr; $chr =~ s/chr010/chr10/; }
    my @coord = ( $chr, $matchedEnd, $matchedStart );
    return(\@coord);
  }
  else {
    my @coord = ();
    return(\@coord);
  }
}

#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Printer;

##About: Filter the output of Quota Align to be more readible.
# Note: if importing file from CoGE you will need to convert the '||' delimiters into ','.

#Run like: perl ~/families/scripts/frac/coge_filtering/filter_Qalign_out.pl -i quota_aln_results.txt -o quota_aln_resultsLite.txt

my $usage = "\n\n$0 -i <quota_aln_results.txt> -o <quota_aln_resultsLite.txt>\n\n";

our ($opt_i, $opt_o, $opt_h);
getopts("i:o:h") or die "";

if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (defined $opt_h)) {die "$usage";}

open (my $infile, '<', $opt_i) or die "\nCannot open $opt_i $!\n\n";
open (my $outfile, '>', $opt_o) or die "\nCannot open $opt_o $!\n\n";

my $block_count = 0;
while (my $line = <$infile> ) {
  ++$block_count;
	chomp $line;

  #Skip the final output of quota align that gives coverage stats
  if ($line =~ /coverage/) {
    next;
  }

  my (%hash, $gene_count, $block_id, $firstaStart, $firstqStart, $lastaEnd, $lastqEnd);

  if ( $line !~ /^#/ ) {
    # This must be the first entry in a new block
    my (undef, $anchorMeta, undef, undef, undef, $queryMeta, undef, undef, $eVal, undef) = split('\t', $line);
    my ($aChr, $aStart, $aEnd, $aGene, $aStrand, $aType, $aDbID, undef, $aPID) = split('\|\|', $anchorMeta);
    my ($qChr, $qStart, $qEnd, $qGene, $qStrand, $qType, $qDbID, undef, $qPID) = split('\|\|', $queryMeta);
    # I don't want all info associated with block, just this relevant information
    my $info_line = join("\t", $aChr,$aStart,$aEnd,$aGene,$qChr,$qStart,$qEnd,$qGene);
    # push this info onto an array to print once you come to the end of the block
    push @{$hash{'info'}}, $info_line;
    # for the unique block identifier I want some info about the starting coordinates of the block
    $block_id = join('-', $block_count,$aChr,$qChr);
    $firstaStart = $aStart;
    $firstqStart = $qStart;
    $gene_count++;
  }

  # read though the lines of a specific block
  while ( my $block_line = <$infile> ) {
    chomp $block_line;
    # keep going unless you reach the start of the next block
    if ( $block_line !~ /^#/ ) {
      # keep track of the number of genes in the block
      ++$gene_count;
      # capture some relevant information
      my (undef, $anchorMeta, undef, undef, undef, $queryMeta, undef, undef, $eVal, undef) = split('\t', $block_line);
      my ($aChr, $aStart, $aEnd, $aGene, $aStrand, $aType, $aDbID, undef, $aPID) = split('\|\|', $anchorMeta);
      my ($qChr, $qStart, $qEnd, $qGene, $qStrand, $qType, $qDbID, undef, $qPID) = split('\|\|', $queryMeta);

      # Do this here for the very first block in the file
      if ( $block_count == 1 && $gene_count == 1) {
        $block_id = join('-', $block_count,$aChr,$qChr);
        $firstaStart = $aStart;
        $firstqStart = $qStart;
      }
      # I also want the end coordinates so will assign this each time and overwrite
      $lastaEnd = $aEnd;
      $lastqEnd = $qEnd;
      # I don't want all info associated with block, just this relevant information
      my $info_line = join("\t", $aChr,$aStart,$aEnd,$aGene,$qChr,$qStart,$qEnd,$qGene);
      # push this info onto an array to print once you come to the end of the block
      push @{$hash{'info'}}, $info_line;
    }
    # if you reached the end of the block, then exit loop and print info to output
    else {
      print $outfile "#$block_id\t$gene_count\t$firstaStart-$lastaEnd\t";
      if ( $firstqStart > $lastqEnd ) {
        print $outfile "$lastqEnd-$firstqStart\n";
      }
      else {
        print $outfile "$firstqStart-$lastqEnd\n";
      }
      for ( my $i = 0; $i<$gene_count; $i++ ) {
        print $outfile "$hash{info}[$i]\n";
      }
      last;
    }
  }
}
print "Num Blocks Seen: ",$block_count-1,"\n";
close $infile;
close $outfile;

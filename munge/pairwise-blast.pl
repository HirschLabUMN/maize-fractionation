#!/usr/bin/perl
use strict;
use Getopt::Std;
use Data::Printer;
use Bio::DB::Fasta;
use Bio::Index::Fasta;
use Bio::Tools::Run::StandAloneBlastPlus;


##ABOUT: This script takes clusters of duplicated genes after running through coordinate fill/bt2 realignment QC
# of DFGs. It will extract the clusters and put them into a hash of array data structure.
# Then for each cluster it will extract fasta files and make a mini blast database. This requires that you have
# a large fasta with all of the sequences you will need pre-made ahead of time.
# Once extracted, it will run tblastx on each cluster and print a summary of the top hits.
##RUNLIKE: perl munge/run_tblastx.pl -c duplicate_clusters.txt -s anchor_duplicates.fasta -b b73_duplicates.fasta -p ph207_duplicates.fasta -o tblastx_results.txt

our ($opt_i, $opt_q, $opt_s, $opt_r, $opt_o);
getopts("i:o:q:s:r:");

#Bio::Seq object
my $queryIndexName = $opt_q;
my $subjIndexName = $opt_s;

open ( my $infile, '<', $opt_i);
open ( my $outfile, '>', $opt_o);

print $outfile "header\n";

my @chars = ("A".."Z", "a".."z");
my $string;
$string .= $chars[rand @chars] for 1..8;

my $queryIdx = Bio::Index::Fasta->new($queryIndexName);
my $subjIdx = Bio::Index::Fasta->new($subjIndexName);

while ( my $line = <$infile> ) {
  chomp $line;
  my ($id1, $id2) = split '\t', $line;
  my $seq1 = $queryIdx->fetch($id1);
  my $seq2 = $subjIdx->fetch($id2);
  my $factory = Bio::Tools::Run::StandAloneBlastPlus->new();

  $factory->bl2seq(
                -method=>'blastn',
                -outfile=> $string,
                -query=> $seq1,
                -subject=> $seq2,
                -outformat=> 5,
                -method_args=> [ '-task' => 'blastn' ]);
  $factory->cleanup;

  my ($result,$hit);
  # Read-in BLAST output
  my $in = new Bio::SearchIO(-format => 'blastxml', -file => $string );

  # Loop through results (should only be one, since this is pairwise)
  while( $result = $in->next_result ) {
    my $id1 = $seq1->display_id;
    my $id2 = $seq2->display_id;
    # Get hits as long as they are defined
    if (!defined $result->next_hit) {
      print $outfile $id1,",",$id2,"\tno_hits_found";
    }
    else {
      #Go back to top
      $result->rewind;
      # Get defined hits
      while( $hit = $result->next_hit ) {
        print $outfile $id1,",",$id2,"\t";
        print $outfile $result->query_length,",",$hit->length,",";
        print $outfile $hit->frac_identical(),",",$hit->frac_aligned_query(),",",$hit->frac_aligned_hit(),"\t";
        my $hsp = $hit->next_hsp;
        print $outfile $hsp->length('total'),",",$hsp->percent_identity,",",$hit->significance . "\t" . $opt_r;
      }
    }
  }
  print $outfile "\n";
}
close $infile;
close $outfile;

#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;
use List::Util 'any';
use Data::Printer;

##ABOUT: This script was used to find how many DFGs hit to gene on scaffolds
##RUNLIKE: [/home/hirschc1/shared/projects/fractionation/data/blast_results/scaffold_gene_search] % perl /home/hirschc1/shared/projects/fractionation/munge/run_blastp.pl -p B73_dfGenes.fasta -d PH207 -t tmp.txt -o B73toPH207_blasp_results.txt

our ($opt_p, $opt_d, $opt_t, $opt_o, $opt_h);
getopts("i:p:d:t:o:c:h") or die "";


## READ-IN FASTA FILE OF SEQUENCES ##
my $seqio_obj = Bio::SeqIO->new(-file => $opt_p, -format => "fasta" );
open (my $outfile, '>', $opt_o) or die;

# Set-up BLAST factory
my $factory = Bio::Tools::Run::StandAloneBlastPlus->new( -db_name => $opt_d );

print $outfile "SeqID\tHitRank,QueryLen,HitPerID,FracAlnQuery,FracAlnHit\t";
print $outfile "AlnLen,HSPPerID,Eval\t";
print $outfile "Hit,HitLoc\n";
while ( my $seq1 = $seqio_obj->next_seq ) {
  my $seqid = $seq1->display_id;
  $seqid =~ s/lcl\|//;

  my $match = 0;
  # Call blastn and define parameters.
  $factory->blastp(
    -query => $seq1,
    -method_args => [ '-task' => 'blastp', '-max_target_seqs' => 1, '-num_threads' => 10, '-evalue' => 1e-30  ],
    -outformat => '5', #output is XML
    -outfile => $opt_t
  );
  #Get rid of temp files created through BLAST
  $factory->cleanup;

  ## PARSE THE OUTPUT OF THE ABOVE BLAST COMMAND ##
  # Define some global variables
  my ($result,$hit);

  # Read-in BLAST output
  my $in = new Bio::SearchIO(-format => 'blastxml', -file   => $opt_t );

  while( $result = $in->next_result ) {
    # Get hits as long as they are defined
    if (!defined $result->next_hit) {
      last;
    }
    else {
      #Go back to top
      $result->rewind;
      # Get defined hits
      while( $hit = $result->next_hit ) {
        my $hit_chr = $hit->name;
        my $hsp = $hit->next_hsp;
        if ( $hit->frac_identical >= 0.70 && $hit->frac_aligned_query > 0.40 ) {
        	print $outfile $seqid . "\t";

        	print $outfile $hit->rank . ",";
        	print $outfile $hit->query_length . ",";
        	print $outfile $hit->frac_identical . ",";
        	print $outfile $hit->frac_aligned_query . ",";
        	print $outfile $hit->frac_aligned_hit . "\t";

        	print $outfile $hit->length_aln  . ",";
        	print $outfile $hsp->percent_identity . ",";
        	print $outfile $hit->significance . "\t";

        	print $outfile $hit_chr . ",";
        	print $outfile $hsp->hit->start . ",";
        	print $outfile $hsp->hit->end . "\n";
      	}
      }
    }
  }
}
close $outfile;
exit 0;

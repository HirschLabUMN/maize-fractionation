#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlastPlus;
use Storable;
use YAML::XS qw/LoadFile/;
use Data::Printer;
##ABOUT: This script takes clusters of duplicated genes after running through coordinate fill/bt2 realignment QC
# of DFGs. It will extract the clusters and put them into a hash of array data structure.
# Then for each cluster it will extract fasta files and make a mini blast database. This requires that you have
# a large fasta with all of the sequences you will need pre-made ahead of time.
# Once extracted, it will run tblastx on each cluster and print a summary of the top hits.
##RUNLIKE: perl munge/run_tblastx.pl -c duplicate_clusters.txt -s anchor_duplicates.fasta -b b73_duplicates.fasta -p ph207_duplicates.fasta -o tblastx_results.txt

our ( $opt_c, $opt_s, $opt_b, $opt_p, $opt_o );
getopts("c:s:b:p:o:");
open ( my $outfile, '>', $opt_o );

my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";
my $in;

# Get storable database
my %b73Hash = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die "tried loading b73_gff storable but couldn't. If using a different genome you will need to update hardocoding in this script\n";
my %ph207Hash = %{ retrieve($settings->{'db'}->{'ph207'}->{'storable'}) } or die "tried loading ph207_gff storable but couldn't. If using a different genome you will need to update hardocoding in this script\n";
my %sbHash = %{ retrieve($settings->{'db'}->{'sb'}->{'storable'}) } or die;
my %osHash = %{ retrieve($settings->{'db'}->{'os'}->{'storable'}) } or die;

#  Check to see if this file exists already, if so delete it so old results aren't appended
if (-e 'anchor_sequences.fa')
{
`rm anchor_sequences.fa`;
}
if (-e 'maize_sequences.fa')
{
`rm maize_sequences.fa`;
}

# DEFINE SOME VARIABLES
my $sbfile = $opt_s;
my $b73file = $opt_b;
my $ph207file = $opt_p;
my $SbDb = Bio::DB::Fasta->new($sbfile);
my $B73db = Bio::DB::Fasta->new($b73file);
my $PH207db = Bio::DB::Fasta->new($ph207file);
my $href = &retrieveClusters($opt_c);
my ( @Sb, @B73, @PH207 );
# LOOP THROUGH CLUSTERS
foreach my $key ( keys %$href ) {
  # Empty arrays each time through the loop
  @Sb = (); @B73 = (); @PH207 = ();
  for ( my $i = 0; $i <= $#{ $href->{$key}->{'a'} }; $i++ ) {
    push @Sb, $href->{$key}->{'a'}->[$i];
  }
  for ( my $i = 0; $i <= $#{ $href->{$key}->{'b'} }; $i++ ) {
    if ( $href->{$key}->{'b'}->[$i] =~ /:/ ) {
      next;
    }
    if ( $href->{$key}->{'b'}->[$i] =~ /,/ ) {
      my @Bs = split ',', $href->{$key}->{'b'}->[$i];
      push @B73, $_ foreach @Bs;
    }
    else {
      push @B73, $href->{$key}->{'b'}->[$i];
    }
  }
  for ( my $i = 0; $i <= $#{ $href->{$key}->{'p'} }; $i++ ) {
    if ( $href->{$key}->{'p'}->[$i] =~ /:/ ) {
      next;
    }
    if ( $href->{$key}->{'p'}->[$i] =~ /,/ ) {
      my @Ps = split ',', $href->{$key}->{'p'}->[$i];
      push @PH207, $_ foreach @Ps;
    }
    else {
      push @PH207, $href->{$key}->{'p'}->[$i];
    }
  }

  ## BLAST COMMANDS START HERE ##
  foreach my $sb ( @Sb ) {
    my $anchrRep;
    if ( $sb =~ /Sobic/ ) {
      $anchrRep = $sbHash{$sb}{'rep_trans'};
    }
    else {
      $anchrRep = $osHash{$sb}{'rep_trans'};
    }
    &makeMultiFa($anchrRep, $SbDb, 'anchor_sequences.fa'); # Retrive sequence from larger database
  }
  foreach my $bgene ( @B73 ) { # starting with B73...
    my $b73rep = $b73Hash{$bgene}{'rep_trans'};
    &makeMultiFa($b73rep, $B73db, 'maize_sequences.fa');
  }
  foreach my $pgene ( @PH207 ) {
    my $ph207rep = $ph207Hash{$pgene}{'rep_trans'};
    &makeMultiFa($ph207rep, $PH207db, 'maize_sequences.fa');
  }
  &makeBlsDb('anchor_sequences.fa'); # Make a blast database from these seqs
  &runBls('maize_sequences.fa'); # Now, read in maize multifasta and set-up blast
  `rm anchor_sequences.fa`;
  `rm tmp_blast_output.bls`;
  `rm maize_sequences.fa`;
  `rm tmpBlastDb.n*`;
}
close $outfile;

# THIS IS TO COLLECT ALL OF THE CLUSTERS
# RETURN HASH->Cluster1=>'b'=>[b73gene1, b73gene2, etc.]
sub retrieveClusters {
  my ( $handle ) = @_;
  open ( my $fh, '<', $handle) or die;
  my ( %hash, $cluster, @b_genes, @p_genes, @anchors );
  while (my $line = <$fh> ) {
	 chomp $line;
    if ( $line =~ /Cluster\s(\d+)/ ) {
      $cluster = "Cluster" . $1;
    }

    while ( $line = <$fh> ) {
    # if this is the end of the block or the end of the file then exit
      if ( $line =~ /#/ ) {
        # since this is the end you should know the full gene-span
        $hash{$cluster}{'b'} = [ @b_genes ];
        $hash{$cluster}{'p'} = [ @p_genes ];
        $hash{$cluster}{'a'} = [ @anchors ];
        if ( $line =~ /Cluster\s(\d+)/ ) {
          $cluster = "Cluster" . $1;
          @b_genes = (); @p_genes = (); @anchors = ();
        }
        last;
      }
    # if this is not the end of the block...
      chomp $line;
      my @genes = split '\t', $line;
      push @anchors, shift @genes;
      foreach my $gene ( @genes ) {
        if ( $gene =~ /Zm00001d/ ) {
          push @b_genes, $gene unless grep{$_ eq $gene} @b_genes;
       }
        else {
          push @p_genes, $gene unless grep{$_ eq $gene } @p_genes;
        }
      }
    }
    # To prevent missing the last line of file enter this loop
    if (eof) {
      $hash{$cluster}{'b'} = [ @b_genes ];
      $hash{$cluster}{'p'} = [ @p_genes ];
      $hash{$cluster}{'a'} = [ @anchors ];
      last;
    }
  redo if $line =~ /#/;
  }
  return \%hash;
}

sub makeMultiFa {
  my ( $in, $db, $out_fh ) = @_;
  my $dbObj = $db->get_Seq_by_id($in);
  my $seqObj = Bio::Seq->new(-seq=>$dbObj->seq,
                             -display_id=>$dbObj->display_id);
  my $seqioObj = Bio::SeqIO->new(-file => ">>$out_fh",
                                 -format=>'fasta');
  $seqioObj->write_seq($seqObj);
}

sub makeBlsDb {
  my ( $fa ) = @_;
  my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(
    -db_name => 'tmpBlastDb',
    -db_data => $fa,
    -create => 1
  );
  $fac -> make_db();
}

sub runBls {
  my ( $maizeFasta ) = @_;
  my $factory = Bio::Tools::Run::StandAloneBlastPlus->new( -db_name => 'tmpBlastDb' );
  my $seqio_obj = Bio::SeqIO->new(-file => $maizeFasta, -format => "fasta" );
  while ( my $seq = $seqio_obj->next_seq ) {
    $factory->tblastx(
      -query => $seq,
      -method_args => [ '-num_threads' => 1  ],
      -outformat => '5', #output is XML
      -outfile => 'tmp_blast_output.bls'
    );
    $factory->cleanup;
    # Now, we're going to parse the blast output
    my $in = new Bio::SearchIO(-format => 'blastxml', -file => 'tmp_blast_output.bls' );
    my ( $result, $hit );
    while( $result = $in->next_result ) {
      print $outfile $result->query_description;
      while( $hit = $result->next_hit ) {
        print $outfile "\t",$hit->name,"\t",$hit->raw_score;
      }
    }
    print $outfile "\n";
  }
}

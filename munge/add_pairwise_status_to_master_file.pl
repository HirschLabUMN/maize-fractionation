#!/usr/bin/perl -w
use strict;
use feature qw/say/;
use Getopt::Long qw(GetOptions);
use Storable;
use Library;
use Data::Printer;

##ABOUT: This script will add a confidence status to each line. It will also add a column with the fractionation class.
##RUNLIKE: perl add_pairwise_status_to_master_file.pl -i ../cache/coge/Zm.vs.Sb_designated_wOsgenes_duplicates_split.txt -b ../data/blast_results/concat_bl2seq_scored.txt -p ../data/blast_results/concat_blsnRealn_scored.txt -n ../data/blast_results/concat_blsnNA_scored.txt -o ../cache/coge/masterFile_wStatus.txt

my ( %pwiseHash, %genewiseHash, %chrwiseHash );
my ( $master, $bl2seq, $blastn, $blastna, $output );
GetOptions(
  'master=s' => \$master,
  'blast2seq:s' => \$bl2seq,
  'blastn:s' => \$blastn,
  'blastNA=s' => \$blastna,
	'output=s' => \$output,
) or die;

open ( my $infile, '<', $master ) or die "couldn't open $master\n";
open ( my $naFile, '<', $blastna ) or die "couldn't open $blastna\n";
open ( my $outfile, '>', $output ) or die "couldn't open $output\n";

my ( $bl2seqFile, $blsnFile );
if ( $bl2seq ) {
  open ( $bl2seqFile, '<', $bl2seq ) or die "couldn't open $bl2seq\n";
  open ( $blsnFile, '<', $blastn ) or die "couldn't open $blastn\n";
}

#################################
## STORE BLAST RESULTS IN HASH ##
#################################

## BL2SEQ ALIGNMENT
if ( $bl2seq ) {
  while ( my $line = <$bl2seqFile> ) {
    chomp $line;
    my ( $g1, $g2, $perID, $g1Cov, $g2Cov, $status, $run ) = split '\t', $line;
    $g1 =~ s/_.*//; $g2 =~ s/_.*//;
    if ( $status eq 'PASS' ) {
      $pwiseHash{$g1}{$g2}{'status'} = $status;
      $pwiseHash{$g2}{$g1}{'status'} = $status;
    }
    elsif ( $status eq 'FAIL' && $perID < 0.75 ) {
      next;
    }
    elsif ( $status eq 'FAIL' && $g1Cov > 0.50 ) {
      $pwiseHash{$g2}{$g1}{'status'} = 'PASS';
    }
    elsif ( $status eq 'FAIL' && $g2Cov > 0.50 ) {
      $pwiseHash{$g1}{$g2}{'status'} = 'PASS';
    }
    else {
      next;
    }
  }
  close $bl2seqFile;
}

## BLASTN REALIGNMENT
if ( $blastn ) {
  while ( my $line = <$blsnFile> ) {
    chomp $line;
    my ( $g1, undef, undef, undef, $alnStatus, $loc, $run, $sub ) = split '\t', $line;
    $g1 =~ s/_.*//;
    if ( $alnStatus eq 'PASS' ) {
      $genewiseHash{$g1}{$run}{$sub}{'status'} = 'PASS';
      $genewiseHash{$g1}{$run}{$sub}{'loc'} = $loc;
    }
  }
  close $blsnFile;
}

##  BLASTN NA ALIGNMENT
while ( my $line = <$naFile> ) {
  if ( $line =~ /no_hits_found/ || $line !~ /Zm/ ) {
    next;
  }
  chomp $line;
  my ( $g1, undef, undef, $locStatus, $alnStatus, $loc, $run, $sub ) = split '\t', $line;
  $g1 =~ s/_.*//;
  if ( $locStatus eq 'PASS' && $alnStatus eq 'PASS' ) {
    $chrwiseHash{$g1}{$run}{$sub}{'status'} = 'PASS';
    $chrwiseHash{$g1}{$run}{$sub}{'loc'} = $loc;
  }
}
close $naFile;

#####################################################
## READ-THRU MASTER FILE AND FIND MAPPING CATEGORY ##
#####################################################

<$infile>;
while ( my $line = <$infile> ) {
  chomp $line;
  my @objects; my @laps = qw/bm1 bm2 pm1 pm2/;
  my ( $anchr, $b1, $b2, $p1, $p2 ) = split '\t', $line;
  foreach my $lap (@laps) {
		my ( $gene, $horiz, $vert, $diag, $geno, $opposite, $subgenome, $oppSubgenome, $finalGene ); # initially set statuses as 0
		# Define some default variables
		if ($lap eq 'bm1') {
      $gene = $b1; $vert = $p1; $horiz = $b2; $diag = $p2; $geno = 'B73'; $opposite = 'PH207'; $subgenome = 'maize1'; $oppSubgenome = 'maize2';
    }
    elsif ( $lap eq 'bm2' ) {
      $gene = $b2; $vert = $p2; $horiz = $b1; $diag = $p1; $geno = 'B73'; $opposite = 'PH207'; $subgenome = 'maize2'; $oppSubgenome = 'maize1';
    }
    elsif ( $lap eq 'pm1' ) {
      $gene = $p1; $vert = $b1; $horiz = $p2; $diag = $b2; $geno = 'PH207'; $opposite = 'B73'; $subgenome = 'maize1'; $oppSubgenome = 'maize2';
    }
    elsif ( $lap eq 'pm2' ) {
      $gene = $p2; $vert = $b2; $horiz = $p1; $diag = $b1; $geno = 'PH207'; $opposite = 'B73'; $subgenome = 'maize2'; $oppSubgenome = 'maize1';
    }
    else {
      die "unknown lap\n";
    }

    if ( $bl2seq && $gene ne 'NA' ) {
      if ( $gene =~ /:/ ) {
        $finalGene = $gene;
      }
      elsif ( $gene ne 'NA' ) {
        $finalGene = &scorePairwise( $gene, $horiz, $vert, $geno, $opposite, $subgenome, $oppSubgenome );
      }
      else {
        $finalGene = &scoreNA( $gene, $horiz, $vert, $diag, $geno, $opposite, $subgenome, $oppSubgenome );
      }
    }
    else {
      if ( $gene eq 'NA' ) {
        $finalGene = &scoreNA( $gene, $horiz, $vert, $diag, $geno, $opposite, $subgenome, $oppSubgenome );
      }
      else {
        $finalGene=$gene;
      }
    }
    push @objects, $finalGene;
  }

  print $outfile "$anchr";
  foreach my $i (0 .. $#objects) {
    print $outfile "\t$objects[$i]";
  }
  print $outfile "\n";
}

#################
## SUBROUTINES ##
#################

sub scorePairwise {

  my ( $gene, $horizontal, $vertical, $genotype, $opp, $sub, $oppSub ) = @_;
  my ( $verticalStat, $verticalLoc, $horizontalStat, $horizontalLoc, $diagonalStat, $diagonalLoc );
  # If Pm1 and Bm2 are present get bl2seq results
  my ($genotypeInfo, $subgenomeInfo);
  my $pairwiseRef = \%pwiseHash;
  my $genewiseRef = \%genewiseHash;
  if ( $vertical eq 'NA' && $horizontal eq 'NA' ) {
    return($gene);
  }

  # 1.) Check vertical
  if ( exists($pairwiseRef->{$gene}->{$vertical}) && $pairwiseRef->{$gene}->{$vertical}->{'status'} eq 'PASS' ) {
    return($gene);
  }
  # 2.) Check horizontal
  if ( exists($pairwiseRef->{$gene}->{$horizontal}) && $pairwiseRef->{$gene}->{$horizontal}->{'status'} eq 'PASS' ) {
    return($gene);
  }
  # 3.) Check vertical to gene space  ( Did B2 align to P2? )
  $genotypeInfo = $opp . "v" . $genotype; # If B1: B73 vs. PH207
  $subgenomeInfo = $sub . "v" . $sub; # If B1: maize1 vs. maize1
  if ( exists($genewiseRef->{$vertical}->{$genotypeInfo}->{$subgenomeInfo}) && $genewiseRef->{$vertical}->{$genotypeInfo}->{$subgenomeInfo}->{'status'} eq 'PASS' ) {
    return($genewiseRef->{$vertical}->{$genotypeInfo}->{$subgenomeInfo}->{'loc'});
  }
  # 4.) Check horizontal to gene space ( Did B2 align to B1? )
  $genotypeInfo = $genotype . "v" . $genotype; # If B1: B73 vs. PH207
  $subgenomeInfo = $oppSub . "v" . $sub; # If B1: maize1 vs. maize1
  if ( exists($genewiseRef->{$horizontal}->{$genotypeInfo}->{$subgenomeInfo}) && $genewiseRef->{$horizontal}->{$genotypeInfo}->{$subgenomeInfo}->{'status'} eq 'PASS' ) {
    return($genewiseRef->{$horizontal}->{$genotypeInfo}->{$subgenomeInfo}->{'loc'}); #return($genewiseRef->{$g1}->{$genotypeInfo}->{$subgenomeInfo}->{'loc'});
  }
  return('NA');
}

sub scoreNA {
  my ( $gene, $horizontal, $vertical, $diagonal, $genotype, $opp, $sub, $oppSub ) = @_;
  my ( $vertNAloc, $horizNAloc, $diagNAloc );
  my $geneFinal = 'NA';

  $vertNAloc = naBls($vertical, $opp . "v" . $genotype, $sub . "v" . $sub); # did the P1 gene map ?
  if ( $vertNAloc ne 'NA' ) {
    return($vertNAloc);
  }

  $horizNAloc = naBls($horizontal, $genotype . "v" . $genotype, $oppSub . "v" . $sub); # did the P1 gene map ?
  if ( $horizNAloc ne 'NA' ) {
    return($horizNAloc);
  }
  $diagNAloc = naBls($diagonal, $opp . "v" . $genotype, $oppSub . "v" . $sub); # did the P1 gene map ?
  if ( $diagNAloc ne 'NA' ) {
    return($diagNAloc);
  }
  return($geneFinal);
}

sub naBls {
  my ( $g1, $genoInfo, $subInfo ) = @_;
  my $chrwiseRef = \%chrwiseHash;

  if ( exists($chrwiseRef->{$g1}->{$genoInfo}->{$subInfo}) && $chrwiseRef->{$g1}->{$genoInfo}->{$subInfo}->{'status'} eq 'PASS') {
      return($chrwiseRef->{$g1}->{$genoInfo}->{$subInfo}->{'loc'});
  }
  else {
    return('NA');
  }
}

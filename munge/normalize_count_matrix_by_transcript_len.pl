#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Data::Printer;
use List::Util qw(sum);
use Sort::Naturally;
use Storable;

##ABOUT: This script will convert raw htseq counts to a normalized length count so that I can compare homologous genes from B73 and PH207 when mapping B73 to B73 and PH207 to PH207. The idea is the next step, you would import this matrix into DESeq2 and then carry out the differential expression analysis.
## the format of the matrix is put every B73 and PH207 gene model as the first column and then the tissues for both across the top with 0's all the way across for the tissues corresponding to the opposite genotype.
## 						B73 tissue1 B73 tissue2 PH207 tissue1 PH207 tissue2
## B73_gene1			4						5						0							0
## PH_gene2				0						0						10						16
## THE OUTPUT FORMAT IS DIFFERENT. I concatenate the B73 gene model and PH207 gene models together here.
##RUNLIKE: perl ${CODE}/rnaseq/normalize_count_matrix_by_transcript_len.pl -c /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus.txt -i /home/hirschc1/shared/projects/fractionation/cache/rnaseq/B73andPH207_mRNA_htseq_count_matrices.txt -o /home/hirschc1/shared/projects/fractionation/cache/rnaseq/B73andPH207_mRNA_htseq_count_matrices_lengthCorrected.txt2

my $usage = "\n\n$0 -i <INPUT> -o <OUTPUT>\n\n";

our ($opt_i, $opt_o, $opt_c, $opt_h);
getopts("i:o:c:h") or die "couldn't get options";

open (my $coge_file, '<', $opt_c) or die "Cannot open $opt_c $!\n\n";
open (my $infile, '<', $opt_i) or die "Cannot open $opt_i $!\n\n";
open (my $outfile, '>', $opt_o) or die "Cannot open $opt_o $!\n\n";
print $outfile "Gene\tB73bd1\tB73bd2\tB73bd3\tB73cp1\tB73cp2\tB73cp3\tB73gk1\tB73gk2\tB73gk3\tB73rt1\tB73rt2\tB73rt3\tB73sd1\tB73sd2\tB73sd3\tB73st1\tB73st2\tB73st3\tPH207bd1\tPH207bd2\tPH207cp1\tPH207cp2\tPH207gk1\tPH207gk2\tPH207rt1\tPH207rt2\tPH207sd1\tPH207sd2\tPH207st1\tPH207st2\n";


my %Bhash = %{ retrieve('/home/hirschc1/shared/projects/fractionation/config/b73_gff.storable')};
my %Phash = %{ retrieve('/home/hirschc1/shared/projects/fractionation/config/ph207_gff.storable')};

## GET SYNTELOG ASSIGNMENTS ##
my $href = &getAssignments($coge_file);
## FETCH UNCORRECTED COUNTS ##
my $HoAref = &getRawCounts($infile);
close $coge_file;
close $infile;

#######################
## LOOP AND CORRECT ##
######################

# TRAVERSE ASSIGNMENTS
foreach my $gene ( keys %{ $href } ) {
	my ($b73_gene, $ph207_gene) = split '_', $gene;

	if ( !exists($HoAref->{$b73_gene}) ) {
		warn "$b73_gene, $Bhash{$b73_gene}{'chr'} was skipped.\n";
		next;
	}
	if ( !exists($HoAref->{$ph207_gene}) ) {
		warn "$ph207_gene, $Phash{$ph207_gene}{'chr'} was skipped.\n";
		next;
	}


	print $outfile "$gene";

  # Get rep transcript lengths from storable db
  my $b_len = $Bhash{$b73_gene}{'cds_length'};
  my $p_len = $Phash{$ph207_gene}{'cds_length'};
  # calculate a ratio to correct by
  # Take the average gene length divided by individual gene length
  my $ratioB = $href->{$gene} / $b_len;
  my $ratioP = $href->{$gene} / $p_len;

  # the first 18 samples are from b73 ( 6 tissues x 3 reps )
  for my $i ( 0 .. 17 ) {
    # no need to correct here
    my $corrected_count = $HoAref->{$b73_gene}->[$i] > 0 ? $HoAref->{$b73_gene}->[$i] * $ratioB : '0';
    #print $outfile "\t$HoA{$b73_gene}[$i]";
    printf $outfile ("\t%.0f", $corrected_count);
  }
  # the next 12 samples are from ph207 ( 6 tissues x 2 reps )
  for my $p ( 18 .. $#{ ${$HoAref}{$ph207_gene} } ) {
    my $corrected_count = $HoAref->{$ph207_gene}->[$p] > 0 ? $HoAref->{$ph207_gene}->[$p] * $ratioP : '0';
    #print $outfile "\t$HoA{$ph207_gene}[$p]";
    printf $outfile ("\t%.0f", $corrected_count);
  }
  print $outfile "\n"
}
close $outfile;
exit 0;

#################
## SUBROUTINES ##
#################

sub getAvrLen {
  my ( $aref, $Bhashref, $Phashref ) = @_;
  my ( $totallen, $count ) = 0;

  foreach my $gene ( @{ $aref } ) {
    if ( $gene =~ /Zm00001d/ ) {
      my $len = ${$Bhashref}{$gene}{'cds_length'};
      $totallen += $len;
      ++$count;
    }
    elsif ( $gene =~ /Zm00008a/ ) {
      my $len = ${$Phashref}{$gene}{'cds_length'};
      $totallen += $len;
      ++$count;
    }
    else {
      next;
    }
  }
  my $avr_len = $totallen / $count;
  return($avr_len);
}

sub getAssignments {
  my ( $handle ) = @_;
  my %assigns;
  while ( my $line = <$handle> ) {
    chomp $line;
    my ( undef, $b1, $b2, $p1, $p2, undef ) = split '\t', $line;
    my @genes = ($b1, $b2, $p1, $p2);
    my $avrLen = &getAvrLen(\@genes, \%Bhash, \%Phash);
    if ( $b1 =~ /Zm/ && $p1 =~ /Zm/ ) {
      my $assign = join('_', $b1,$p1);
      $assigns{$assign} = $avrLen;
    }
    if ( $b2 =~ /Zm/ && $p2 =~ /Zm/ ) {
      my $assign = join('_', $b2,$p2);
      $assigns{$assign} = $avrLen;
    }
  }
  return(\%assigns);
}

sub getRawCounts {
  my ($input_fh) = @_;
  my %HoA;
  my $header = <$input_fh>;
  while ( my $line = <$input_fh> ) {
    chomp $line;

    my @fields = split '\t', $line;
    my $trans = shift @fields;
    my $gene = $trans;
    $gene =~ s/_T\d+//;
    if ( $gene =~ /^Zm00001d/ ) {
        $HoA{$gene} = [ @fields ];
      }
    if ( $gene =~ /^Zm00008a/ ) {
        $HoA{$gene} = [ @fields ];
    }
  }
  return(\%HoA);
}

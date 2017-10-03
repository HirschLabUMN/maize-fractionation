#!/usr/bin/perl -w
use strict;
use Library;
use Getopt::Std;
use Data::Printer;
use Storable;
##ABOUT: This script will parse the results of the tblastx runs but traversing the master file and identifying the duplicate lines. It will then look to see if the blast results support the assignment. If so it keeps the gene model and if not it replaces the gene model with 'NA'.
##RUNLIKE: perl ../../munge/parse_tblastx_results.pl -i /home/hirschc1/shared/projects/fractionation/cache/coge/tblastx_results.txt -m /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC2.txt -o /home/hirschc1/shared/projects/fractionation/cache/coge/masterFile_wStatus_filteredQC3.txt

our ( $opt_i, $opt_m, $opt_o );
getopts("i:m:o:");

open ( my $mfile, '<', $opt_m );
open ( my $outfile, '>', $opt_o );

my %bHash = %{ retrieve('/home/hirschc1/shared/projects/fractionation/config/b73_gff.storable')};
my %pHash = %{ retrieve('/home/hirschc1/shared/projects/fractionation/config/ph207_gff.storable')};
my %sHash = %{ retrieve('/home/hirschc1/shared/projects/fractionation/config/sb_gff.storable')};
my %oHash = %{ retrieve('/home/hirschc1/shared/projects/fractionation/config/os_gff.storable')};

my ( $href, $anchref ) = &storeResults( $opt_i );

while ( my $line = <$mfile> ) {
  chomp $line;
  my @fields = split '\t', $line;
  my $anchr = shift @fields;
  # was this anchor in the cluster file i.e. a duplicate? if so, then inspect
  if ( exists($anchref->{$anchr}) ) {
    # Loop through B73m1 B73m2 PH207m1 and PH207m2
    for ( my $i = 0; $i < $#fields+1; $i++ ) {
      if ( $i == 0 ) {
        print $outfile "$anchr";
      }
      # Did you store this gene as a hit ?
      if ( exists($href->{$fields[$i]})) {
        # If gene is equal to anchor that it was a best hit and belongs on this line
        if ( $href->{$fields[$i]}->{'hit'} eq $anchr ) {
          print $outfile "\t$fields[$i]";
        }
        else {
          # this gene mapped to another sorghum gene better. replace the gene with a coordinate
          my $coord = &getCoordinate($fields[$i]);
          #print $outfile "\t$$coord";
          print $outfile "\tNA"; # on second thought...replace with NA
        }
      }
      else {
        # this must be the case of tandem duplicates. Zm00008a013568,Zm00008a013567 both map to SbA
        if ( $fields[$i] =~ /Zm/ ) {
          # read both genes into an array
          my @subfields = split ',', $fields[$i];
          # Now we gonna see which is the better match and keep it
          my $match = 'NA';  my $score = 0;
          foreach my $subfield ( @subfields ) {
            # should be match here...
            if ( $href->{$subfield}->{'hit'} eq $anchr) {
              # if score is better than previous or default then this is the gene I like
              if ( $href->{$subfield}->{'score'} > $score) {
                $score = $score;
                $match = $subfield;
              }
            }
          }
          print $outfile "\t$match";
        }
        else {
          # If this gene was NA check to see if you had a match from a different gene...
          my $subgenome;
          # Get the subgenome you're looking for
          if ( $i == 0 ) {
            $subgenome = 'qmaize1';
          }
          elsif ( $i == 1 ) {
            $subgenome = 'qmaize2';
          }
          elsif ( $i == 2 ) {
            $subgenome = 'smaize1';
          }
          elsif ( $i == 3 ) {
            $subgenome = 'smaize2';
          }
          else {
            die "unknown counter\n";
          }
          if ( $anchr =~ /Sobic/ && $anchref->{$anchr}->{$subgenome}->{'hit'} && $anchref->{$anchr}->{$subgenome}->{'hit'} ne 'NA'  ) {
            print $outfile "\t",$anchref->{$anchr}->{$subgenome}->{'hit'};
          }
          else {
            print $outfile "\t$fields[$i]";
          }
        }
      }
    }
    print $outfile "\n";
  }
  else {
    # If this wasn't a problesome line, then simply print out the line.
    print $outfile "$line\n";
  }
}
close $outfile;

sub storeResults {
  my ( $handle ) = @_;
  open ( my $fh, '<', $handle) or die;
  my ( %hash, %anchors );
  while ( my $line = <$fh> ) {
    chomp $line;
    my @fields = split '\t', $line;
    my $query = shift @fields;
    $query =~ s/_T.*//;
    $hash{$query}{'hit'} = 'NA';
    $hash{$query}{'score'} = 0;
    for ( my $i = 0; $i < $#fields+1; $i = $i+2 ) {
      my $anchor;
      if ( $fields[$i] =~ /(.*\d+)\.\d/ ) {
        $anchor = $1;
      }
      my $subgenome = &getAssignment($anchor, $query);
      # Automatically store every anchor so I can tell it was a duplicate in next step
      $anchors{$anchor}{$subgenome}{'score'} = 'NA';
      $anchors{$anchor}{$subgenome}{'hit'} = 'NA';
      my $score = $fields[$i+1]; # the bit score should be the next array element
      # if this anchor has the highest score mark it
      # if its equal to the previous then store both
      if ( $score > $hash{$query}{'score'} ) {
        $hash{$query}{'score'} = $score;
        $hash{$query}{'hit'} = $anchor;
        if ( $subgenome ne 'foo' ) {
          if ( $anchors{$anchor}{$subgenome}{'hit'} ne 'NA' && $score > $anchors{$anchor}{$subgenome}{'score'} ) {
            $anchors{$anchor}{$subgenome}{'score'} = $score;
            $anchors{$anchor}{$subgenome}{'hit'} = $query;
          }
          else {
            $anchors{$anchor}{$subgenome}{'score'} = $score;
            $anchors{$anchor}{$subgenome}{'hit'} = $query;
          }
        }
      }
      elsif ( $score == $hash{$query}{'score'} ) {
        $hash{$query}{'hit'} = $fields[$i] . "-" . $hash{$query}{'hit'};
      }
      else {
        next;
      }
    }
  }
  return( \%hash, \%anchors);
}

sub getCoordinate {
  my ( $gene ) = @_;
  my $coord;
  if ( $gene =~ /Zm00001d/ ) {
    $coord = $bHash{$gene}{'chr'} . ":" . $bHash{$gene}{'start'} . "-" . $bHash{$gene}{'stop'};
  }
  elsif ( $gene =~ /Zm00008a/ ) {
    $coord = $pHash{$gene}{'chr'} . ":" . $pHash{$gene}{'start'} . "-" . $pHash{$gene}{'stop'};
  }
  else {
    die "Error getting coordinate for $gene\n";
  }
  return(\$coord);
}

sub getAssignment {
  my ( $anchrGene, $maizeGene ) = @_;

  my ( $sp, $maizeSp, $anchrChr, $maizeChr, $parsedAssignment );
  if ( $anchrGene =~ /Sobic/ ) {
    $sp = 'Sb';
    $anchrChr = $sHash{$anchrGene}{'chr'};
  }
  else {
    $sp = 'Os';
    $anchrChr = $oHash{$anchrGene}{'chr'};
  }
  if ( $maizeGene =~ /Zm00001d/ ) {
    $maizeChr = $bHash{$maizeGene}{'chr'};
    $maizeSp = 'q';
  }
  else {
    $maizeChr = $pHash{$maizeGene}{'chr'};
    $maizeSp = 's';
  }
  my $assignment = &Library::syntenicBlockAssignment($anchrChr, $maizeChr, $sp);

  if ( $assignment eq 'undef' ) {
    return('foo');
  }
  else {
    $parsedAssignment = $maizeSp . $assignment;
    return($parsedAssignment);
  }
}

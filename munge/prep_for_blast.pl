#!/usr/bin/perl -w
use strict;
use feature qw/say/;
use Getopt::Std;
use Storable;
use YAML::XS qw/LoadFile/;

##ABOUT: This will grab gene that can be aligned pairwise. Must combine the 'designated' files using `combine_condensed_syntelog_files.pl` before running.
##RUNLIKE: perl ../../munge/prep_for_blast_ii.pl -i Zm.vs.Sb_designated.txt -o /scratch.global/brohamm1/projects/frac/B73v4vPH207_B73ids.txt -p /scratch.global/brohamm1/projects/frac/B73v4vPH207_PH207ids.txt

our ($opt_i, $opt_o, $opt_r);
getopts("i:o:r:") or die;

my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";

# Get storable database
my (%bhash, %phash);
if ( $opt_r eq 'B73vB73' ) {
  %bhash = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die;
  %phash = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die;
}
elsif ( $opt_r eq 'B73vPH207'  ) {
  %bhash = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die;
  %phash = %{ retrieve($settings->{'db'}->{'ph207'}->{'storable'}) } or die;
}
elsif ( $opt_r eq 'PH207vPH207' ) {
  %bhash = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die;
  %phash = %{ retrieve($settings->{'db'}->{'ph207'}->{'storable'}) } or die;
}
else {
  exit "unknown run supplied to option r should be: B73vB73, B73vPH207, or PH207vPH207\n";
}

open (my $infile, '<', $opt_i) or die;
open (my $outfile, '>', $opt_o) or die;

my $run = $opt_r;

<$infile>;
while ( my $line = <$infile> ) {
  chomp $line;

  if ( $line =~ /^#/ ) {
    next;
  }

  my ($anchr, $m1, $m2, $m3, $m4) = split '\t', $line;

  my @m1genes = split ',', $m1;
  my @m2genes = split ',', $m2;
  my @m3genes = split ',', $m3;
  my @m4genes = split ',', $m4;

  if ( $run eq 'B73vPH207' ) {
    foreach my $m1gene ( @m1genes ) {
      foreach my $m3gene ( @m3genes ) {
        unless ( ( $m1gene eq 'NA' || $m1gene =~ /:/ ) || ( $m3gene eq 'NA' || $m3gene =~ /:/ ) ) {
          say $outfile "$bhash{$m1gene}{'rep_cds'}","\t","$phash{$m3gene}{'rep_cds'}";
        }
      }
    }
    foreach my $m4gene ( @m4genes ) {
      foreach my $m2gene ( @m2genes ) {
        unless ( ( $m2gene eq 'NA' || $m2gene =~ /:/ ) || ( $m4gene eq 'NA' || $m4gene =~ /:/ ) ) {
          say $outfile "$bhash{$m2gene}{'rep_cds'}","\t","$phash{$m4gene}{'rep_cds'}";
        }
      }
    }
  }
  elsif ( $run eq 'B73vB73' ) {
    foreach my $m1gene ( @m1genes ) {
      foreach my $m2gene ( @m2genes ) {
        unless ( ( $m1gene eq 'NA' || $m1gene =~ /:/ ) || ( $m2gene eq 'NA' || $m2gene =~ /:/ ) ) {
          say $outfile "$bhash{$m1gene}{'rep_cds'}","\t","$bhash{$m2gene}{'rep_cds'}";
        }
      }
    }
  }
  elsif ( $run eq 'PH207vPH207' ) {
    foreach my $m4gene ( @m4genes ) {
      foreach my $m3gene ( @m3genes ) {
        unless ( ( $m3gene eq 'NA' || $m3gene =~ /:/ ) || ( $m4gene eq 'NA' || $m4gene =~ /:/ ) ) {
          say $outfile "$phash{$m3gene}{'rep_cds'}","\t","$phash{$m4gene}{'rep_cds'}";
        }
      }
    }
  }
  else {
    die "\nrun format must be 1.) B73vPH207 or 2.) B73vB73, or 3.) PH207vPH207\n\b";
  }
}
close $outfile;
exit 0;

#! /usr/bin/perl -w
use strict;
use Getopt::Long();
use Data::Printer;
use Sort::Naturally;

##ABOUT: This script will make a matrix from expression data
##RUNLIKE:

my $usage = "\nUsage: $0 --matrix_out <output file for matrix> --workdir <WORKDIR->sample1->counts.htseq> --help <help on usage of script>\n\n";

my ($matrix_out, $workdir, @samples, $help);

Getopt::Long::GetOptions(
  'workdir=s' => \$workdir,
  'matrix_out=s' => \$matrix_out,
  'h|help' => \$help
);

if (defined($help)) {
    die $usage;
}

# Open output file
open (my $matrix_fh, '>', $matrix_out) or die "\nCan't open file $matrix_out\n\n";
my %hoh;

my @files = <$workdir/*/*CDS.extended.htseq>;
die "no *htseq files found\n" if (@files < 0);
foreach my $file (nsort @files) {
  open (my $fh, '<', $file) or die "error opening $file\n";
    my @samp = split '/', $file;

    # Get a list of samples for printing headers
    push @samples, $samp[-1];
    while (my $line = <$fh>) {
      chomp $line;
      my ($model, $count) = split '\t', $line;
      if ( $model !~ /(Zm|:)/ ) {
        next;
      }
      $hoh{$model}{$samp[-1]} = $count;
    }
  close $fh;
}

print $matrix_fh "Gene";
p @samples;
foreach my $s (nsort @samples) {
  print $matrix_fh "\t$s";
}

foreach my $key (nsort keys %hoh) {
  print $matrix_fh "\n$key";
  foreach my $key2 (nsort keys %{ $hoh{$key} } ) {
    print $matrix_fh "\t$hoh{$key}{$key2}"
  }
}

close $matrix_fh;
exit;

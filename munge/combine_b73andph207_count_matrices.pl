#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Data::Printer;
use List::Util qw(sum);
use Sort::Naturally;

##ABOUT: This script will convert normalized counts from HTseq into RPKMs. The denominator is the total number of unique reads mapped. This is regardless of whether they mapped to a feature (the representative transcript in this case) or not. This requires providing a separate file with these values to this script.
##RUNLIKE: perl ../../munge/rnaseq/normCounts2rpkm.pl -i b73v4ref_six_tissue_htseqcount_matrix_normalized.txt -o rpkm_normalized_matrix.txt -c uniq_counts_for_rpkm_conversion.txt

my $usage = "\n\n$0 -i <INPUT> -o <OUTPUT>\n\n";

our ($opt_i, $opt_o, $opt_p, $opt_h);
getopts("i:o:p:h") or die "";

open (my $infile, '<', $opt_i) or die "Cannot open $opt_i $!\n\n";
open (my $infile2, '<', $opt_p) or die "Cannot open $opt_p $!\n\n";
open (my $outfile, '>', $opt_o) or die "Cannot open $opt_o $!\n\n";

my $head1 = <$infile>;
my $head2 = <$infile2>;
chomp $head1; chomp $head2;
$head2 =~ s/Gene\t//;
print $outfile "$head1\t$head2\n";

while ( my $line = <$infile> ) {
	chomp $line;
	my @fields = split '\t', $line;
	print $outfile "$line\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
}
while ( my $line2 = <$infile2> ) {
	chomp $line2;
	my @fields = split '\t', $line2;
	my $gene = shift @fields;
	$gene =~ s/_T01.v1.1//;
	#print $outfile "$gene\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t"; #for p v p and p v b
	print $outfile "$gene\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";
	my $i = 0;
	map { print $outfile $_ . ($i == $#fields ? "\n" : "\t") ; $i++ } @fields;
}

close $infile;
close $infile2;
close $outfile;

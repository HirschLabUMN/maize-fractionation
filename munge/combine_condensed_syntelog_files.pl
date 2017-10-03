#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Data::Printer;
use Library;

##About: This files combines syntelog information for Sb-B73 and Sb-PH207 into a single file
# File1: Sb B73_1 B73_2 B73_3   File2: Sb PH207_1 PH207_2 PH207_3
# Outcome: Sb B73_1 B73_2 B73_3 PH207_1 PH207_2 PH207_3

# Run like: perl ~/families/scripts/frac/coge_filtering/combine_condensed_syntelog_files.pl -i ~/families/analysis/frac/coge/B73v2vSb/condensed_syntelogs_wGEvo_parsed.txt -p ~/families/analysis/frac/coge/PH207vSb/condensed_syntelogs_wGEvo_parsed_coded.txt -o ~/families/analysis/frac/syntelogs/combined_syntelogs.txt`

my $usage = "\n\n$0 -i <INPUT> -o <OUTPUT>\n\n";

our ($opt_i, $opt_p, $opt_o, $opt_h);
getopts("i:p:o:h") or die "";

if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (defined $opt_h)) {die "$usage";}

if (-e $opt_o) {
  die "\nThe output file, $opt_o, exists.\n\n";
}

open (my $infile, '<', $opt_i) or die "Cannot open $opt_i $!\n\n"; # B73 condensed syntelogs
open (my $infile2, '<', $opt_p) or die "Cannot open $opt_p $!\n\n"; # PH207 condensed syntelogs
open (my $outfile, '>', $opt_o) or die "Cannot open $opt_o $!\n\n";

## Store Sorghum--B73 Assignments from first file ##

<$infile>; # skip header
my %hash;
while (my $line = <$infile> ) {
	chomp $line;

  if ( $line =~ /^#/ ) {
    next;
  }

	# Define columns
	my ($sb, $m1, $m2, undef) = split("\t", $line);
  my $fm1 = Library::b73format($m1);
  my $fm2 = Library::b73format($m2);
  my $fsb = Library::b73format($sb);

	# Store Sb-Maize assigments and add NA for PH207 genes
	$hash{$fsb}{'Bm1'} = $fm1;
	$hash{$fsb}{'Bm2'} = $fm2;
	$hash{$fsb}{'Pm1'} = 'NA';
	$hash{$fsb}{'Pm2'} = 'NA';

}
close $infile;

## Store Sorghum--PH207 assignments from second file ##

#skip header
<$infile2>;
while (my $line = <$infile2> ) {
	chomp $line;

  if ( $line =~ /^#/ ) {
    next;
  }

	my ($sb, $m1, $m2, undef) = split("\t", $line);
  my $fm1 = Library::b73format($m1);
  my $fm2 = Library::b73format($m2);
  my $fsb = Library::b73format($sb);

	# Did you see this sorhum gene in the first file?
	# If so, just update PH207 information
	if (exists($hash{$fsb})) {
		$hash{$fsb}{'Pm1'} = $fm1;
		$hash{$fsb}{'Pm2'} = $fm2;
	}
	# If not, make B73 entries be NA
	else {
		$hash{$fsb}{'Pm1'} = $fm1;
		$hash{$fsb}{'Pm2'} = $fm2;
		$hash{$fsb}{'Bm1'} = 'NA';
		$hash{$fsb}{'Bm2'} = 'NA';
	}
}
close $infile2;

print $outfile "Sorghum\tB73_Maize1\tB73_Maize2\tPH207_Maize1\tPH207_Maize2\n";
foreach my $key (keys %hash) {
	if ($hash{$key}{'Bm1'} eq "NA" && $hash{$key}{'Bm2'} eq "NA" && $hash{$key}{'Pm1'} eq "NA" && $hash{$key}{'Pm2'} eq "NA" ) {
		next;
	}
	else {
		print $outfile "$key\t$hash{$key}{'Bm1'}\t$hash{$key}{'Bm2'}\t$hash{$key}{'Pm1'}\t$hash{$key}{'Pm2'}\n";
	}
}

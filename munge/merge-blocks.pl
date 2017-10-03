#!/usr/bin/perl -w
use strict;
use List::Util qw/max/;
use Data::Printer;
##ABOUT: This script will merge genomic coordinates given chromsome start stop. Two contiguous regions are merged.
# This was run to determine the amount of coding sequence occupied by syntenic blocks in B73 and PH207.
##RUNLIKE: perl ../../munge/merge-blocks.pl B73vSb_block_list.txt > B73vSb_merged_blocks.txt
# perl ../../munge/merge-blocks.pl B73vOs_block_list.txt > B73vOs_merged_blocks.txt
# cat B73vSb_merged_blocks.txt B73vOs_merged_blocks.txt >> B73vSbandOs_merged_blocks.txt
# perl ../../munge/merge-blocks.pl B73vSbandOs_merged_blocks.txt > B73vSbandOs_merged_blocks_final.txt

open (my $infile, '<', $ARGV[0]);

my %chroms;
while ( <$infile> ) {
	chomp;
	my ($chrom, $start, $end, $sub) = split ' ';
	# Anonymous hash of arrays of arrays
	push @{$chroms{$chrom}{$sub}}, [$start, $end];
}
close $infile;

for my $chrom ( sort keys %chroms ) {
	for my $sub (sort keys %{$chroms{$chrom}}) {
		my $ranges = $chroms{$chrom}{$sub};
		#		p $ranges; #[0] [ [0] start [1] stop ]
		#		p @{$ranges}; # ranges is a reference to an array of arrays
		@{$ranges} = sort {$a->[0] <=> $b->[0]} @{$ranges};
		#print "${$ranges}[0][0]\n";
		for my $i ( 0 .. $#$ranges-1 ) {
			if ( $ranges->[$i][1] >= $ranges->[$i+1][0] - 1) {
				$ranges->[$i+1][0] = $ranges->[$i][0];
				$ranges->[$i+1][1] = max($ranges->[$i][1], $ranges->[$i+1][1]);
				$ranges->[$i] = undef; # assign element that used to hold range as undef so its removed
			}
		}
		#p @{$ranges};
		@{$ranges} = grep {$_} @$ranges; # gets rid of undefineds
		## Loop through the array of arrays for this chromsome and print contents ##
		for my $p ( 0 .. $#$ranges ) {
			print "$chrom";
			for my $q ( 0 .. $#{ $ranges->[$p] } ) {
				print "\t$ranges->[$p]->[$q]";
			}
			print "\t$sub\n";
		}
	}
}

#!/usr/bin/perl
# Print the lengths of sequences
# Copyright 2012 Shaun Jackman

use strict;

while (<>) {
	next if /^#/;
	die unless /^>/;
	chomp;
	my ($id, $info) = split '\t', $_;
	$id =~ s/>//; $id =~ s/_T\d+//;
	my $seq = '';
	while (<>) {
		next if /^#/;
		last if /^>/;
		chomp;
		$seq .= $_;
	}

	my $len = $seq =~ tr/ACGTNacgtn//;
	my $ncount = $seq =~ tr/Nn//;
	print "$id\t$len\n";
#	printf "$id\t$info\t%.2f\n", $ncount/$len; 
	redo if /^>/;
	last if eof;
}

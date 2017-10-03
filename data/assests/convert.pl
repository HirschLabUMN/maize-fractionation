#!/usr/bin/perl -w
use strict; 
open ( my $in, '<', $ARGV[0] ) or die; 
my %hash;
while ( my $line = <$in> ) {
	chomp $line;
	my @fields = split '\t', $line; 
	my ( $v3, $v4 ) = @fields[0,4];
	$hash{$v3} = $v4; 
}
close $in; 

open ( my $in2, '<', $ARGV[1] ); 
while ( my $line = <$in2> ) {
	chomp $line; 
	my @fields = split '\t', $line; 
	if ( exists($hash{$fields[0]}) ) {
		print "$hash{$fields[0]}\n" unless $hash{$fields[0]} =~ /:/;
	}
}
close $in2; 

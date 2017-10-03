#!/usr/bin/perl -w
use strict;
use Storable;
use Data::Printer;
use Getopt::Std;
use YAML::XS qw/LoadFile/;

our ($opt_i, $opt_o);
getopts("i:o:");

my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";

my $hashref;
if ( $opt_i eq 'b73' ) {
	$hashref = retrieve($settings->{'db'}->{'b73'}->{'storable'});
}
elsif ( $opt_i eq 'ph207' ) {
	$hashref = retrieve($settings->{'db'}->{'ph207'}->{'storable'});
}
else {
	die "$opt_i should be 'b73' or 'ph207'\n\n";
}

open (my $out, '>', $opt_o);
foreach my $key ( keys %$hashref ) {
	if ( $hashref->{$key}->{'chr'} =~ /(ctg|scaffold)/ ) {
		my $rep_cds = $hashref->{$key}->{'rep_cds'};
		if ( $rep_cds =~ /Zm00008a/ ) {
			$rep_cds =~ s/T/P/;
		}
		print $out "$rep_cds\n";
	}
}
close $out;

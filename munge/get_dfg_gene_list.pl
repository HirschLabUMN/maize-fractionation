#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Storable;
use YAML::XS qw/LoadFile/;
my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml');

our ($opt_i, $opt_m, $opt_o);
getopts("i:m:o:");

my %hash = %{ retrieve($settings->{'db'}->{$opt_i}->{'storable'}) } or die;

open ( my $master, '<', $opt_m ) or die;
open ( my $out, '>', $opt_o ) or die;

while ( my $line = <$master> ) {
  chomp $line;
  my @fields = split '\t', $line;
  my ( $rep_cds );
  if ( $opt_i eq 'b73' ) {
  	if ( $fields[1] ne 'NA' && $fields[3] eq 'NA' && $fields[1] =~ /Zm/) {
  	  $rep_cds = $hash{$fields[1]}{'rep_cds'};
	  #$rep_cds =~ s/T/P/;
   	  print $out "$rep_cds\n";
  	}
  	elsif ( $fields[2] ne 'NA' && $fields[4] eq 'NA' && $fields[2] =~ /Zm/) {
    	  $rep_cds = $hash{$fields[2]}{'rep_cds'};
    	  #$rep_cds =~ s/T/P/;
    	  print $out "$rep_cds\n";
  	}
  }
  else {
	if ( $fields[1] eq 'NA' && $fields[3] ne 'NA' && $fields[3] =~ /Zm/) {
    	  $rep_cds = $hash{$fields[3]}{'rep_cds'};
    	  $rep_cds =~ s/T/P/;
    	  print $out "$rep_cds\n";
  	}
  	elsif ( $fields[2] eq 'NA' && $fields[4] ne 'NA' && $fields[4] =~ /Zm/) {
    	  $rep_cds = $hash{$fields[4]}{'rep_cds'};
    	  $rep_cds =~ s/T/P/;
    	  print $out "$rep_cds\n";
  	}
  }
}
close $master;
close $out;
exit 0;

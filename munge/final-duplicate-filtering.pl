#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Printer;
use Storable;
use YAML::XS qw/LoadFile/;


my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";

my %b73Hash = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die "tried loading b73_gff storable but couldn't. If using a different genome you will need to update hardocoding in this script\n";
my %ph207Hash = %{ retrieve($settings->{'db'}->{'ph207'}->{'storable'}) } or die "tried loading ph207_gff storable but couldn't. If using a different genome you will need to update hardocoding in this script\n";
my %sbHash = %{ retrieve($settings->{'db'}->{'sb'}->{'storable'}) } or die "tried loading sb_gff storable but couldn't. If using a different genome you will need to update hardocoding in this script\n";
my %osHash = %{ retrieve($settings->{'db'}->{'os'}->{'storable'}) } or die "tried loading os_gff storable but couldn't. If using a different genome you will need to update hardocoding in this script\n";

our ($opt_i, $opt_f, $opt_o);
getopts("i:f:o:") or die;

my $dupRef = getDuplicates($opt_i);

open ( my $infile, '<', $opt_i ) or die;
open ( my $outfile, '>', $opt_o ) or die;
open ( my $fused, '>', $opt_f ) or die;
my %dups;
while ( my $line = <$infile> ) {
	chomp $line;
	my @fields = split '\t', $line;
	foreach my $field ( @fields ) {
		if (exists($dupRef->{$field}) ) {
			my @geneOrders;
			my $duplicate1 = $dupRef->{$field}->[0];
			my $duplicate2 = $dupRef->{$field}->[1];
			# If a sorghum and rice gene share a duplicate then just print the sorghum gene and get rid of rice gene
			if ( $duplicate1 !~ /Sobic/ ) {
				$dups{$duplicate1} = '';
				next;
			}
			if ( $duplicate1 !~ /Sobic/ ) {
				$dups{$duplicate2} = '';
				next;
			}

			# If two sorghum genes share a duplicate, then see if they are consecutive genes
			foreach my $duplicate ( @{$dupRef->{$field}}) {
				my $order = &getOrder($duplicate);
				push @geneOrders, $order;
			}

			# If consecutive genes, then probably a gene fusion issue so filter out to separate file
			my $distance = abs ( $geneOrders[1] - $geneOrders[0] );
			if ( $distance <= 3 ) {
				$dups{$duplicate1} = '';
				$dups{$duplicate2} = '';
				print $fused "$line\n";
				last;
			}
			# If not consecutive, then probably just a false assignment and something I can take care of by re-blasting
		}
	}
	print $outfile "$line\n" unless exists($dups{$fields[0]});
}
close $infile;
close $outfile;
close $fused;

sub getDuplicates {
  my ( $handle ) = @_;
	my ( %seen, %duplicates );
  open ( my $in, '<', $handle ) or die;
	my $i = 0;
	{
		$i++;
    while ( my $line = <$in> ) {
      chomp $line;
      my @fields = split '\t', $line;
      my @genes = @fields[1,2,3,4];
      foreach my $gene ( @genes ) {
				if ( $i == 1 ) {
					++$seen{$gene}{'count'} unless $gene eq 'NA' || $gene =~/:/;
					push @{$seen{$gene}{'ortholog'}}, $fields[0] unless $gene eq 'NA' || $gene =~/:/;
				}
				else {
					if ( (exists($seen{$gene}) && $seen{$gene}{'count'} > 1) || $line =~ /,Zm/ ) {
						$duplicates{$gene} = [ @{$seen{$gene}{'ortholog'}} ];
					}
				}
			}
    }
		seek $in, 0, 0;
		redo if $i < 2;
	}
	close $in;
	return(\%duplicates);
}

sub getOrder {
	my ( $ortholog ) = @_;
	my $ord;
	if ( $ortholog =~ /Sobic/ ) {
		$ord = $sbHash{$ortholog}{'order'}; # Check Sorghum storable
	}
	else {
		$ord = $osHash{$ortholog}{'order'}; # Check Oryza storable
	}
	return($ord);
}

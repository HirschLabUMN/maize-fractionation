#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Printer;
use Storable;
use YAML::XS qw/LoadFile/;
##ABOUT: This is small script to collect duplicate genes from the master file after incorporating the coordinate fills and bowtie2 realignments.
##RUNLIKE: perl filter_duplicates.pl tmp1.txt masterFile_wStatus_filteredQC2.txt > tmp2.txt
#NOTE: See QC_master_file.md for info on creating tmp files.

our ($opt_i,$opt_s,$opt_o,$opt_b,$opt_p,$opt_c);
getopts("i:c:s:o:b:p:") or die;

my $dupRef = getDuplicates($opt_i);

#Sobic.006G257001	Zm00001d001929	NA	Zm00008a006012	Zm00008a038796
#Sobic.006G257100	Zm00001d001929	Zm00001d026546	NA	NA

my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";
my $in;

# Get storable database
my %b73Hash = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die "tried loading b73_gff storable but couldn't. If using a different genome you will need to update hardocoding in this script\n";
my %ph207Hash = %{ retrieve($settings->{'db'}->{'ph207'}->{'storable'}) } or die "tried loading ph207_gff storable but couldn't. If using a different genome you will need to update hardocoding in this script\n";
my %sbHash = %{ retrieve($settings->{'db'}->{'sb'}->{'storable'}) } or die "tried loading sb_gff storable but couldn't. If using a different genome you will need to update hardocoding in this script\n";
my %osHash = %{ retrieve($settings->{'db'}->{'os'}->{'storable'}) } or die "tried loading os_gff storable but couldn't. If using a different genome you will need to update hardocoding in this script\n";

#  Check to see if this file exists already, if so delete it so old results aren't appended
if (-e 'anchor_duplicates.fa')
{
`rm anchor_duplicates.fa`;
}

my ( %lineHash, %geneHash, @anchor_duplicates );
foreach my $key  ( keys %$dupRef ) {
	my @fields = split '\t', $key;
	my $anchr = $fields[0];
	my @genes = @fields[ 1 .. 4 ];
	my @finalGenes;

	foreach my $gene ( @genes ) {
		my @geneSubset = split ',', $gene;
		foreach my $subGene ( @geneSubset ) {
			push @finalGenes, $subGene;
		}
	}

	my @filtGenes = grep { $_ ne 'NA' } @finalGenes;
	# SorghumA { [0] Bm1 [1] Bm1 [2] Pm1 [3] Pm2 }
	push @{$lineHash{$anchr}}, $_ foreach @filtGenes;
	foreach my $gene ( @filtGenes ) {
		# Bm1 { [0] SorghumA}; Bm2 { [0] SorghumA , [1] SorghumB }
		push @{$geneHash{$gene}}, $anchr;
	}
}


# Now set up a hash
# this gene appears on x lines
# {SorghumA}{Bm2} = [SorghumA, SorghumB, etc.]
my %hash;
foreach my $line2 ( keys %lineHash ) { #a, b, c
	foreach my $gene ( @{$lineHash{$line2}} ) {
		foreach my $gene2 ( @{$geneHash{$gene}} ) {
			$hash{$line2}{$gene2} = '';
		}
	}
}

# Now make an array of arrays from this hash
# This is all of the anchor genes that contain overlapping genes
# [ [ SorghumA, SorghumB], [SorghumC], [SorghumD, SbE, SbF] ]
my @mainArray;
foreach my $key ( keys %hash ) {
my @subArray;
foreach my $key2 ( keys %{ $hash{$key}} ) {
	push @subArray, $key2;
}
push @mainArray, [@subArray];
}

# This code chunk will iterate through the array of array and do the clustering
ORIG: for my $i (0..$#mainArray) {
	for my $j ($i+1..$#mainArray) {
		my %seen;
		my @unique = grep {! $seen{$_}++} @{$mainArray[$i]}, @{$mainArray[$j]};
		if (@unique < @{$mainArray[$i]} + @{$mainArray[$j]}) {
			$mainArray[$i] = \@unique;
			splice @mainArray, $j, 1;
			redo ORIG;
		}
	}
}

# CREATE THE OUTPUT THAT CONTAINS CLUSTER INFORMATION
open( my $cluster_out, '>', $opt_c) or die;

for my $i ( 0 .. $#mainArray ) {
	print $cluster_out " # Cluster $i\n";
	for my $p ( 0 .. $#{$mainArray[$i]}) {
		my $anchor = $mainArray[$i][$p];
		print $cluster_out "$anchor";
		for my $q ( 0 .. $#{$lineHash{$anchor}}) {
			print $cluster_out "\t$lineHash{$anchor}[$q]";
		}
		print $cluster_out "\n";
	}
}

close $cluster_out;

# WRITE OUT CDS FILES TO MAKE DB FOR EACH GENE SET
open( my $sorghum_out, '>', $opt_s ) or die;
open( my $oryza_out, '>', $opt_o ) or die;
open( my $b73_out, '>', $opt_b ) or die;
open( my $ph207_out, '>', $opt_p ) or die;

for my $i ( 0 .. $#mainArray ) {
	for my $p ( 0 .. $#{$mainArray[$i]}) {
		my $anchor = $mainArray[$i][$p];

		if ( $anchor =~ /Sb/ || $anchor =~ /Sobic/ ) {
			print $sorghum_out "$sbHash{$anchor}{'rep_trans'}\n";
		}
		elsif ( $anchor =~ /Os/ ) {
			print $oryza_out "$osHash{$anchor}{'rep_trans'}\n";
		}
		else {
			warn("Saw unrecognized ortholog, should be Sorghum or Rice: $anchor\n");
			next;
		}
		for my $q ( 0 .. $#{$lineHash{$anchor}}) {
			if ( $lineHash{$anchor}[$q] =~ /Zm00001d/ ) {
				my @genes = split ',', $lineHash{$anchor}[$q];
				foreach my $gene ( @genes ) {
					print $b73_out "$b73Hash{$gene}{'rep_trans'}\n";
				}
			}
			elsif ( $lineHash{$anchor}[$q] =~ /Zm00008a/ ) {
				my @genes = split ',', $lineHash{$anchor}[$q];
				foreach my $gene ( @genes ) {
					print $ph207_out "$ph207Hash{$gene}{'rep_trans'}\n";
				}
			}
			else {
				warn("Saw unrecognized maize gene should be Zm00001d (B73) or Zm00008a (PH207): $lineHash{$anchor}[$q]\n");
				next;
			}
		}
	}
}
close $sorghum_out;
close $oryza_out;
close $b73_out;
close $ph207_out;

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
				my @subGenes = split ',', $gene;
				foreach my $subGene (@subGenes) {
					if ( $i == 1 ) {
						++$seen{$subGene}{'count'} unless $subGene eq 'NA' || $subGene =~/:/;
					}
					else {
						if ( (exists($seen{$subGene}) && $seen{$subGene}{'count'} > 1) || $line =~ /,Zm/ ) {
							$duplicates{$line} = '';
						}
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

#!/usr/bin/perl -w
use strict;
use Storable;
use Getopt::Std;
use YAML::XS qw/LoadFile/;
use Data::Printer;
##ABOUT: This script is used as part of effort to determine how many of the DFGs
# could not be identified because of assembly gaps. After pulling out list of genes from
# R sqlite queries, I need to get 5kb upstream and downstream of gene. This script is used to get These
# coordinates which will then be pulled from the chromosomal multifasta using bedtools getfasta.
##RUNLIKE: perl get_coordinates_for_gap_analysis.pl -i /home/hirschc1/shared/projects/fractionation/data/sqlite/B73nonscaffoldHits_DFGs.txt -o /home/hirschc1/shared/projects/fractionation/data/gap_analysis/B73_DFG_coordinates.txt -g B73

#perl get_coordinates_for_gap_analysis.pl -i /home/hirschc1/shared/projects/fractionation/data/sqlite/PH207nonscaffoldHits_DFGs.txt -o /home/hirschc1/shared/projects/fractionation/data/gap_analysis/PH207_DFG_coordinates.txt -g PH207

our ($opt_i, $opt_m, $opt_o);
getopts("i:m:o:") or die;

my $href = &populateHash($opt_i);
my $masterRef = &getMasterHash($opt_m);

open ( my $master, '<', $opt_m ) or die;

my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";


my %sbHash = %{ retrieve($settings->{'db'}->{'sb'}->{'storable'}) } or die;
my %b73Hash = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die;
my %ph207Hash = %{ retrieve($settings->{'db'}->{'ph207'}->{'storable'}) } or die;

my $orderHashRef = &getOrderHash(\%sbHash);
open ( my $outfile, '>', $opt_o );

while ( my $line = <$master> ) {
	chomp $line;
	my @fields = split '\t', $line;

	if ( $fields[0] !~ /Sobic/ ) {
		print $outfile "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\n";
		next;
	}

	if ( $fields[1] eq 'NA' && $fields[3] =~ /Zm/ ) {
		if ( exists($href->{$fields[3]}) ) {
			my $status = &checkSynteny($fields[0],$fields[3],'maize1');
			if ( $status ne 'NA' ) {
				$fields[1] = $status;
				print "$fields[3]\t$status\n";
			}
		}
	}
	if ( $fields[3] eq 'NA' && $fields[1] =~ /Zm/ ) {
		if ( exists($href->{$fields[1]}) ) {
			my $status = &checkSynteny($fields[0],$fields[1],'maize1');
			if ( $status ne 'NA' ) {
				$fields[3] = $status;
				print "$fields[1]\t$status\n";
			}
		}
	}

	if ( $fields[2] eq 'NA' && $fields[4] =~ /Zm/ ) {
		if ( exists($href->{$fields[4]}) ) {
			my $status = &checkSynteny($fields[0],$fields[4],'maize2');
			if ( $status ne 'NA' ) {
				$fields[2] = $status;
				print "$fields[4]\t$status\n";
			}
		}
	}
	if ( $fields[4] eq 'NA' && $fields[2] =~ /Zm/ ) {
		if ( exists($href->{$fields[2]}) ) {
			my $status = &checkSynteny($fields[0],$fields[2],'maize2');
			if ( $status ne 'NA' ) {
				$fields[4] = $status;
				print "$fields[2]\t$status\n";
			}
		}
	}

	print $outfile "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t$fields[4]\t$fields[5]\n";
}

sub populateHash {
	my ( $handle ) = @_;
	open ( my $infile, '<', $handle);
	my %hash;
	while ( my $line = <$infile> ) {
		chomp $line;
		my ( $query, $hit ) = split '\t', $line;
		push @{$hash{$query}}, $hit;
	}
	close $infile;
	return(\%hash);
}

sub getOrderHash {
  my ( $ordRef ) = @_;
  my %ordHash;
  foreach my $key ( keys %$ordRef ) {
    if ( $ordRef->{$key}->{'synteny'} eq '1' ) {
      my $chr = $ordRef->{$key}->{'chr'};
      my $ord = $ordRef->{$key}->{'order'};
      $ordHash{$chr}{$ord} = $key;
    }
  }
  return(\%ordHash);
}

sub getMasterHash {
	my ( $handle ) = @_;
	open ( my $infile, '<', $handle );
	my %hash;
	while ( my $line = <$infile> ) {
		chomp $line;
		my @fields = split '\t', $line;
		$hash{$fields[0]}{'bmaize1'} = $fields[1];
		$hash{$fields[0]}{'bmaize2'} = $fields[2];
		$hash{$fields[0]}{'pmaize1'} = $fields[3];
		$hash{$fields[0]}{'pmaize2'} = $fields[4];
	}
	close $infile;
	return(\%hash);
}


sub checkSynteny {
  # Sorghum_GENEA
  my ( $anchr, $queryGene, $sub ) = @_;

  my $anchrChr = $sbHash{$anchr}{'chr'};
  my $anchrOrd = $sbHash{$anchr}{'order'};
	my $subgenome; 	my $matchedDistance;
	my $status = 'NA';

	# Loop through all of the potential non-allelic homologs
	foreach my $newGene ( @{$href->{$queryGene}} ) {
		# Get the subgenome
		if ( $newGene =~ /Zm00001d/ ) {
			$subgenome = 'p' . $sub;
		}
		else {
			$subgenome = 'b' . $sub;
		}
  	my $leftBoundary = 'foo'; my $leftDistance; # default init
  	my $rightBoundary = 'foo'; my $rightDistance;


  	# Check upstream
  	for my $i ( 1..50 ) {
    	my $check = $anchrOrd - $i;
			my $cognate;
    	if ( exists($orderHashRef->{$anchrChr}->{$check}) ) {
      	my $matchedSbGene = $orderHashRef->{$anchrChr}->{$check}; # Nearest upstream anchor gene
      	# Is there a maize gene here to check with ?
				if ( !exists($masterRef->{$matchedSbGene}->{$subgenome}) ) { next; }
      	if ( $masterRef->{$matchedSbGene}->{$subgenome} ne 'NA' ) {
        	# Get the maize gene
        	my $matchedZmGene = $masterRef->{$matchedSbGene}->{$subgenome}; # Nearest upstream maize gene
					if ( $matchedZmGene !~ /Zm/ ) { next; }

					# get the order of this gene
					my ( $cognateOrd, $newGeneOrd );
					# If this is an actual B73 gene model...
					if ( $matchedZmGene =~ /Zm00001d/ ) {
						# Get it's cognate
						$cognate = $b73Hash{$matchedZmGene}{'cognate'};
						if ( $cognate !~ /Zm/ ) { # Check to see if the cognate is acutally a gene and not NA
							next;
						}
						# Ok, we have a cognate so get it's gene order as well as the matched genes' order
						$cognateOrd = $ph207Hash{$cognate}{'order'};
						$newGeneOrd = $ph207Hash{$newGene}{'order'};
					}
					else {
						# In this case must be dealing with a PH207 gene...
						my $cognate = $ph207Hash{$matchedZmGene}{'cognate'};
						$cognate = $ph207Hash{$matchedZmGene}{'cognate'};
						if ( $cognate !~ /Zm/ ) { # Make sure there was actually a gene here
							next;
						}
						$cognateOrd = $b73Hash{$cognate}{'order'};
						$newGeneOrd = $b73Hash{$newGene}{'order'};
					}

					# Now compare the newly matched genes' order to the nearest cognates order these should be within 15 genes of each other
					if ( abs($newGeneOrd - $cognateOrd) <= 15 ) {
						$leftBoundary = '1';
						$leftDistance = abs($newGeneOrd - $cognateOrd);
					}
					last;
      	}
      	else {
        	next;
      	}
			}
		}
		# Same thing as above but in opposite direction
		for my $i ( 1..50 ) {
    	my $check = $anchrOrd + $i;
			my $cognate;
    	if ( exists($orderHashRef->{$anchrChr}->{$check}) ) {
      	my $matchedSbGene = $orderHashRef->{$anchrChr}->{$check};
      	# Is there a maize gene here to check with ?
      	if ( $masterRef->{$matchedSbGene}->{$subgenome} ne 'NA' ) {
        	# Get the maize gene
        	my $matchedZmGene = $masterRef->{$matchedSbGene}->{$subgenome};
					if ( $matchedZmGene !~ /Zm/ ) { next; }
					# get the order of this gene
					my ( $cognateOrd, $newGeneOrd );
					if ( $matchedZmGene =~ /Zm00001d/ ) {
						$cognate = $b73Hash{$matchedZmGene}{'cognate'};
						if ( $cognate !~ /Zm/ ) {
							next;
						}
						$cognateOrd = $ph207Hash{$cognate}{'order'};
						$newGeneOrd = $ph207Hash{$newGene}{'order'};
					}
					else {
						my $cognate = $ph207Hash{$matchedZmGene}{'cognate'};
						$cognate = $ph207Hash{$matchedZmGene}{'cognate'};
						if ( $cognate !~ /Zm/ ) { # Make sure there was actually a gene here
							next;
						}
						$cognateOrd = $b73Hash{$cognate}{'order'};
						$newGeneOrd = $b73Hash{$newGene}{'order'};
					}
					if ( abs($newGeneOrd - $cognateOrd) <= 15 ) {
						$rightBoundary = '1';
						$rightDistance = abs($newGeneOrd - $cognateOrd);
					}
					last;
      	}
      	else {
        	next;
      	}
    	}
  	}

		# If you already had a hit...
		if ( $status ne 'NA' ) {
			# If you previously had a hit, but this time thru loop you didn't then skip
			if ( !$leftDistance || !$rightDistance ) {
				next;
			}
			# if you had two matches in the hash of arrays check to see if this new match is closer
			if ( $leftDistance < $matchedDistance || $rightDistance < $matchedDistance) {
				$status = $newGene;
			}
		}
		else  {
  		if ( $leftBoundary eq 'foo' && $rightBoundary eq 'foo') {
    		$status = 'NA';
  		}
			else {
    		$status = $newGene;
				$matchedDistance = $leftDistance;
  		}
		}

	}
	return($status);
}

#!/usr/bin/perl -w
use strict;
use feature qw/say/;
use Getopt::Std;
use Storable;
use Library;
use YAML::XS qw/LoadFile/;

our ($opt_i, $opt_o);
getopts("i:o:");

open my $outfile, '>', $opt_o or die "error opening $opt_o $!\n\n";

my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";
my $in;

# Get storable database
my $sHashRef = retrieve($settings->{'db'}->{'sb'}->{'storable'});
my $oHashRef = retrieve($settings->{'db'}->{'os'}->{'storable'});
my $bHashRef = retrieve($settings->{'db'}->{'b73'}->{'storable'});
my $pHashRef = retrieve($settings->{'db'}->{'ph207'}->{'storable'});
my $sOrderRef = &getOrderHash($sHashRef);
my $oOrderRef = &getOrderHash($oHashRef);


my $resultRef = &populateMaster($opt_i);


foreach my $key ( keys %$resultRef ) {
	if ( $key !~ /Sobic/ ) { next; }

	my ( @result, $upChr, $upStart, $upDown, $downChr, $downStart, $downDown, $varGene, $varCDS, $run );
	# B73 maize1 gone
	if ( $resultRef->{$key}->{'bmaize1'} eq 'NA' && $resultRef->{$key}->{'pmaize1'} =~ /Zm/ ) {
		@result = &getCoordinates($key, 'bmaize1', $sOrderRef);
		$varGene = $resultRef->{$key}->{'pmaize1'};
	}
	# B73 maize2 gone
	if ( $resultRef->{$key}->{'bmaize2'} eq 'NA' && $resultRef->{$key}->{'pmaize2'} =~ /Zm/ ) {
		@result = &getCoordinates($key, 'bmaize2', $sOrderRef);
		$varGene = $resultRef->{$key}->{'pmaize2'};
	}
	# PH207 maize1 gone
	if ( $resultRef->{$key}->{'pmaize1'} eq 'NA' && $resultRef->{$key}->{'bmaize1'} =~ /Zm/ ) {
		@result = &getCoordinates($key, 'pmaize1', $sOrderRef);
		$varGene = $resultRef->{$key}->{'bmaize1'};
	}
	# PH207 maize2 gone
	if ( $resultRef->{$key}->{'pmaize2'} eq 'NA' && $resultRef->{$key}->{'bmaize2'} =~ /Zm/ ) {
		@result = &getCoordinates($key, 'pmaize2', $sOrderRef);
		$varGene = $resultRef->{$key}->{'bmaize2'};
	}

	# Result should contain two elements: an upstream maize gene and a downstream maize gene. 
	if ( @result ) {
		if ( $result[0] eq 'foo' || $result[1] eq 'foo' ) {
			next;
		}
		if ( $result[0] =~ /Zm00001d/ ) {
			$upChr = $bHashRef->{$result[0]}->{'chr'};
			$upStart = $bHashRef->{$result[0]}->{'start'};
			$upDown = $bHashRef->{$result[0]}->{'stop'};
			$downChr = $bHashRef->{$result[1]}->{'chr'};
			$downStart = $bHashRef->{$result[1]}->{'start'};
			$downDown = $bHashRef->{$result[1]}->{'stop'};
			$run = "PH207vB73";
			$varCDS = $pHashRef->{$varGene}->{'rep_cds'};
		}
		else {
			$upChr = $pHashRef->{$result[0]}->{'chr'};
			$upStart = $pHashRef->{$result[0]}->{'start'};
			$upDown = $pHashRef->{$result[0]}->{'stop'};
			$downChr = $pHashRef->{$result[1]}->{'chr'};
			$downStart = $pHashRef->{$result[1]}->{'start'};
			$downDown = $pHashRef->{$result[1]}->{'stop'};
			$run = "B73vPH207";
			$varCDS = $bHashRef->{$varGene}->{'rep_cds'};
		}
		if ( $upChr ne $downChr ) {
			warn "WARNING: chromsomes do not match: $key,$result[0]-$upChr,$result[1]-$downChr\n";
			next;
		}

		if ( $downDown < $upStart ) {
			print $outfile "$varCDS\t$upChr:$downDown-$upStart\t$run\n";
		}
		else {
			print $outfile "$varCDS\t$upChr:$upStart-$downDown\t$run\n";
		}
	}
}
close $outfile;


	sub populateMaster {
	  my ( $handle ) = @_;
	  open ( my $master, '<', $handle ) or die;
	  my %hash;
	  while ( my $line = <$master> ) {
	    chomp $line;
	    my ( $anchr, $m1, $m2, $m3, $m4 )  = split '\t', $line;
			$hash{$anchr}{'bmaize1'} = $m1;
			$hash{$anchr}{'bmaize2'} = $m2;
			$hash{$anchr}{'pmaize1'} = $m3;
			$hash{$anchr}{'pmaize2'} = $m4;
		}
		close $master;
	  return(\%hash);
	}

	sub getOrderHash {
	  my ( $href ) = @_;
	  my %ordHash;
	  foreach my $key ( keys %$href ) {
	    if ( $href->{$key}->{'synteny'} eq '1' ) {
	      my $ord = $href->{$key}->{'order'};
	      my $chr = $href->{$key}->{'chr'};
	      my $st = $href->{$key}->{'start'};
	      my $sp = $href->{$key}->{'stop'};
	      $ordHash{$chr}{$ord}{'st'} = $st;
	      $ordHash{$chr}{$ord}{'sp'} = $sp;
	      $ordHash{$chr}{$ord}{'id'} = $key;
	    }
	  }
	  return(\%ordHash);
	}


	sub getCoordinates {
	  # Sorghum_GENEA
	  my ( $anchr, $subgenome ) = @_;

	  my $upstreamGene = 'foo';
	  my $downstreamGene = 'foo';
		my $anchrChr = $sHashRef->{$anchr}->{'chr'};
		my $anchrOrd = $sHashRef->{$anchr}->{'order'};

	  # Check upstream
	  for my $i ( 1..50 ) {
	    my $check = $anchrOrd - $i;
	    if ( exists($sOrderRef->{$anchrChr}->{$check}) ) {
	      my $matchedSbGene = $sOrderRef->{$anchrChr}->{$check}->{'id'};
	      # Is there a maize gene here to check with ?
				if ( !$resultRef->{$matchedSbGene} ) { next; } # if you filtered it out as a possible fusion
	      if ( $resultRef->{$matchedSbGene}->{$subgenome} =~ /Zm/ ) {
	        # Get the maize gene
	        $upstreamGene = $resultRef->{$matchedSbGene}->{$subgenome};
	        last;
	      }
	      else {
	        next;
	      }
	    }
	  }

		# Check downstream
	  for my $i ( 1..50 ) {
	    my $check = $anchrOrd + $i;
	    if ( exists($sOrderRef->{$anchrChr}->{$check}) ) {
	      my $matchedSbGene = $sOrderRef->{$anchrChr}->{$check}->{'id'};
	      # Is there a maize gene here to check with ?
	      if ( $resultRef->{$matchedSbGene}->{$subgenome} =~ /Zm/ ) {
	        # Get the maize gene
	        $downstreamGene = $resultRef->{$matchedSbGene}->{$subgenome};
	        last;
	      }
	      else {
	        next;
	      }
	    }
	  }
		return($upstreamGene,$downstreamGene);
	}

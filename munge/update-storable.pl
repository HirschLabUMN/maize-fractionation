#!/usr/bin/perl -w
use strict;
use Storable;
use Library;
use Getopt::Std;
use Data::Printer;
use YAML::XS qw/LoadFile/;
#Zm00001d047013
getopts("i:b:p:s:o:");
our ($opt_i,$opt_b,$opt_p,$opt_s,$opt_o);

my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";

# Get storable database
my %primaryHash = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die;
my %secHash = %{ retrieve($settings->{'db'}->{'ph207'}->{'storable'}) } or die;
my %primAncHash = %{ retrieve($settings->{'db'}->{'sb'}->{'storable'}) } or die;
my %secAncHash = %{ retrieve($settings->{'db'}->{'os'}->{'storable'}) } or die;

my ( $primAnchrHashRef, $secAnchrHashRef, $primHashRef, $secHashRef) = &getUpdate($opt_i);

foreach my $key ( keys %primAncHash ) {
	if ( exists($primAnchrHashRef->{$key}) ) {
		$primAncHash{$key}{'synteny'} = '1';
	}
	else {
		$primAncHash{$key}{'synteny'} = 'NA';
	}
}

foreach my $key ( keys %secAncHash ) {
	if ( exists($secAnchrHashRef->{$key}) ) {
		$secAncHash{$key}{'synteny'} = '1';
	}
	else {
		$secAncHash{$key}{'synteny'} = 'NA';
	}
}

foreach my $key ( keys %primaryHash ) {
	if ( exists($primHashRef->{$key}) ) {
		$primaryHash{$key}{'synteny'} = '1';
		$primaryHash{$key}{'homeolog'} = $primHashRef->{$key}->{'homeolog'};
		$primaryHash{$key}{'cognate'} = $primHashRef->{$key}->{'cognate'};
		$primaryHash{$key}{'reciprocal'} = $primHashRef->{$key}->{'reciprocal'};
		$primaryHash{$key}{'ortholog'} = $primHashRef->{$key}->{'ortholog'};
		$primaryHash{$key}{'subgenome'} = $primHashRef->{$key}->{'subgenome'};
	}
	else {
		$primaryHash{$key}{'synteny'} = 'NA';
		$primaryHash{$key}{'homeolog'} = 'NA';
		$primaryHash{$key}{'cognate'} = 'NA';
		$primaryHash{$key}{'reciprocal'} = 'NA';
		$primaryHash{$key}{'ortholog'} = 'NA';
		$primaryHash{$key}{'subgenome'} = 'NA';
	}
}

foreach my $key ( keys %secHash ) {
	if ( exists($secHashRef->{$key}) ) {
		$secHash{$key}{'synteny'} = '1';
		$secHash{$key}{'homeolog'} = $secHashRef->{$key}->{'homeolog'};
		$secHash{$key}{'cognate'} = $secHashRef->{$key}->{'cognate'};
		$secHash{$key}{'reciprocal'} = $secHashRef->{$key}->{'reciprocal'};
		$secHash{$key}{'ortholog'} = $secHashRef->{$key}->{'ortholog'};
		$secHash{$key}{'subgenome'} = $secHashRef->{$key}->{'subgenome'};
	}
	else {
		$secHash{$key}{'synteny'} = 'NA';
		$secHash{$key}{'homeolog'} = 'NA';
		$secHash{$key}{'cognate'} = 'NA';
		$secHash{$key}{'reciprocal'} = 'NA';
		$secHash{$key}{'ortholog'} = 'NA';
		$secHash{$key}{'subgenome'} = 'NA';
	}
}


store(\%primAncHash, $opt_s);
store(\%secAncHash, $opt_o);
store(\%primaryHash, $opt_b);
store(\%secHash, $opt_p);

sub getUpdate {
	my ( $handle ) = @_;
	open ( my $infile, '<', $handle ) or die;
	my ( %sbHash, %osHash, %bHash, %pHash );
	while ( my $line = <$infile> ) {
  	chomp $line;
  	my @objects; my @laps = qw/bm1 bm2 pm1 pm2/;
  	my ( $anchr, $b1, $b2, $p1, $p2 ) = split '\t', $line;
		if ( $anchr =~ /Sobic/ ) {
			$sbHash{$anchr} = '';
		}
		else {
			$osHash{$anchr} = '';
		}

  	foreach my $lap (@laps) {
			my ( $gene, $horiz, $vert, $diag, $subgenome);
    	# Define some default variables
    	if ($lap eq 'bm1') {
      	$gene = $b1; $vert = $p1; $horiz = $b2; $diag = $p2; $subgenome = 'maize1'; ;
    	}
    	elsif ( $lap eq 'bm2' ) {
      	$gene = $b2; $vert = $p2; $horiz = $b1; $diag = $p1; $subgenome = 'maize2';
    	}
    	elsif ( $lap eq 'pm1' ) {
      	$gene = $p1; $vert = $b1; $horiz = $p2; $diag = $b2; $subgenome = 'maize1';
    	}
    	elsif ( $lap eq 'pm2' ) {
      	$gene = $p2; $vert = $b2; $horiz = $p1; $diag = $b1; $subgenome = 'maize2';
    	}
    	else {
      	die "unknown lap\n";
    	}
			if ( $lap eq 'bm1' || $lap eq 'bm2' ) {
				if ( $gene =~ /Zm/ ) {
					$bHash{$gene}{'homeolog'} = $horiz;
					$bHash{$gene}{'cognate'} = $vert;
					$bHash{$gene}{'reciprocal'} = $diag;
					$bHash{$gene}{'ortholog'} = $anchr;
					$bHash{$gene}{'synteny'} = '1';
					$bHash{$gene}{'subgenome'} = $subgenome;
				}
			}
			elsif ( $lap eq 'pm1' || $lap eq 'pm2' ) {
				if ( $gene =~ /Zm/ ) {
					$pHash{$gene}{'homeolog'} = $horiz;
					$pHash{$gene}{'cognate'} = $vert;
					$pHash{$gene}{'reciprocal'} = $diag;
					$pHash{$gene}{'ortholog'} = $anchr;
					$pHash{$gene}{'synteny'} = '1';
					$pHash{$gene}{'subgenome'} = $subgenome;
				}
			}
		}
	}
	close $infile;
	return(\%sbHash, \%osHash, \%bHash, \%pHash);
}

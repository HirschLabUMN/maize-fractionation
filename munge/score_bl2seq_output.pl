#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Printer;
use feature qw/say/;
use Storable;
use Library;
use YAML::XS qw/LoadFile/;
##ABOUT: THIS READS IN THE PARTIALLY PARSED BL2SEQ OUTPUTS AND ASSIGNS EVERYTHING AS PASS OR FAIL
##RUNLIKE: perl score_bl2seq_output.pl -i /home/hirschc1/shared/projects/fractionation/data/blast_results/B73vSb_bl2seq.txt -o /home/hirschc1/shared/projects/fractionation/data/blast_results/B73vSb_bl2seq_scored.txt

my $usage = "\n\n$0 -i <INPUT> -o <OUTPUT>\n\n";

our ($opt_i, $opt_o, $opt_m, $opt_r, $opt_t, $opt_h);
getopts("i:o:m:r:t:h") or die "";

my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";

my $targetFile;
if ( defined($opt_t) ) {
  open ($targetFile, '>', $opt_t) or die "Cannot open $opt_t $!\n\n";
  print $targetFile "header\n";
}

# Get storable database
my (%hashA, %hashB);
my $run_mode = $opt_r;
if ( $run_mode eq 'B73vB73' ) {
  %hashA = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die;
  %hashB = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die;
}
elsif ( $run_mode eq 'B73vPH207'  ) {
  %hashA = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die;
  %hashB = %{ retrieve($settings->{'db'}->{'ph207'}->{'storable'}) } or die;
}
elsif ( $run_mode eq 'PH207vPH207' ) {
  %hashA = %{ retrieve($settings->{'db'}->{'ph207'}->{'storable'}) } or die;
  %hashB = %{ retrieve($settings->{'db'}->{'ph207'}->{'storable'}) } or die;
}
elsif ( $run_mode eq 'PH207vB73' ) {
  %hashA = %{ retrieve($settings->{'db'}->{'ph207'}->{'storable'}) } or die;
  %hashB = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die;
}
else {
  die "unknown run supplied to option r should be: B73vB73, B73vPH207, or PH207vPH207\n";
}

open (my $infile, '<', $opt_i) or die "Cannot open $opt_i $!\n\n";
open (my $outfile, '>', $opt_o) or die "Cannot open $opt_o $!\n\n";

my %hash;
if ( $opt_m eq 'bls2seq' ) {
	say $outfile "B73\tPH207\tPerID\tB73cov\tPH207cov\tAlnLen\tStatus";
}
elsif ( $opt_m eq 'blsn' ) {
	say $outfile "gene\tper_id\tcov\thsp-len\thit-len\tstatus\taln-loc\trun-mode";
}

#skip header
<$infile>;
while (my $line = <$infile> ) {
	chomp $line;

  # Are you processing blast2seq output ?
	if ($opt_m eq 'bls2seq' ) {
		# Get statement and pass/fail status from subroutine
		my ($stmnt, $arrayRef) = &processBl2seq(\$line);
		# print the statement to output
		say $outfile $stmnt;
		# was the status a fail ?
		if ( $arrayRef ne 'NULL' ) {
			# get the reason for failure from subroutine
			# loop through array reference and print both genes to target file if both need to be realigned or just one
			foreach my $string ( @$arrayRef ) {
				print $targetFile "$string\n";
			}
		}
	}
	# Are you processing blastn output ?
	elsif ($opt_m eq 'blsn' ) {
		my ($stmnt) = &processBlsn(\$line);
		say $outfile $stmnt;
	}
	else {
		die "unknown opt_m option\n";
	}
}
close $infile;
close $outfile;
if ( defined($opt_t) ) {
  close $targetFile;
}
exit;

# subroutine to process blast2seq output for summary file that has pass/fail listed next to alignment info
sub processBl2seq {
	my ( $line ) = @_;
	my ( $string, $aref );
	if (${$line} =~ /no_hits_found/) {
		if (${$line} =~ /(.*),(.*)\t(\w+)/) {
			$string = "$1\t$2\tNA\tNA\tNA\tNO_HIT\t$run_mode";
			$aref = &reasonForFail($1,$2,'NA','NA','NA',$run_mode);
		}
	}
	else {
		my ($genes, $meta, undef) = split("\t", ${$line});
		my ($gene1, $gene2) = split ',', $genes;
		my ( undef, undef, $fracIdentical, $fracAlnQuery, $fracAlnHit ) = split ',', $meta;

		if ($fracIdentical >= 0.75 && $fracAlnQuery >= 0.50 && $fracAlnHit >= 0.50)  {
			$string = "$gene1\t$gene2\t$fracIdentical\t$fracAlnQuery\t$fracAlnHit\tPASS\t$run_mode";
			$aref = 'NULL';
		}
		else {
			$string = "$gene1\t$gene2\t$fracIdentical\t$fracAlnQuery\t$fracAlnHit\tFAIL\t$run_mode";
			$aref = &reasonForFail($gene1, $gene2, $fracIdentical, $fracAlnQuery, $fracAlnHit, $run_mode);
		}
	}
	return($string, $aref);
}

# subroutine to process blast2seq output for summary file that has pass/fail listed next to alignment info
sub processBlsn {
	my ($line) = @_;
	my ( $string, $aref);

	if (${$line} =~ /no_matches/ || ${$line} =~ /no_hits_found/) {
    my ($gene, undef) = split '\t', ${$line};
    $gene =~ s/_.*//;
		$string = "$gene\tNA\tNA\tNA\tNA\tNO_HIT\tNA\t$run_mode\tNA\t";
	}
	else {
    my ( $g1, $meta1, $meta2, $loc, $status ) = split '\t', ${$line};
    my ( $rank, $queryLen, $fracIdentical, $fracAlnQuery ) = split ',', $meta1;
    my ($subjChr, undef) = split ':', $loc;
		$g1 =~ s/_.*//; $subjChr =~ s/chr0//; $subjChr =~ s/chr//;
		my $queryChr = $hashA{$g1}{'chr'};
    $queryChr =~ s/chr0//; $queryChr =~ s/chr//;
    my ( $run1, $run2, $run12 );
    if ( $queryChr eq $subjChr ) {
      $run1 = $hashA{$g1}{'subgenome'};
      $run12 = join('v', $run1,$run1);
		}
    else {
      $run1 = &Library::syntenicBlockAssignment($subjChr,$queryChr,'Zm');
      $run2 = &Library::syntenicBlockAssignment($queryChr, $subjChr,'Zm');
      $run12 = join('v', $run1,$run2);
    }
		if ( $fracIdentical >= 0.75 && $fracAlnQuery >= 0.50 )  {
			$string = "$g1\t$fracIdentical\t$fracAlnQuery\t$status\tPASS\t$loc\t$run_mode\t$run12";
		}
    elsif ( $fracIdentical >= 0.75 && $fracAlnQuery > 0.20 && ( $fracAlnQuery*$queryLen > 500 ) ) {
     $string = "$g1\t$fracIdentical\t$fracAlnQuery\t$status\tPASS\t$loc\t$run_mode\t$run12";
    }
		else {
			$string = "$g1\t$fracIdentical\t$fracAlnQuery\t$status\tFAIL\t$loc\t$run_mode\t$run12";
		}
	}
	return($string);
}

# GET REASON FOR FAILURE AND RETURN LINE FOR TARGET FILE
sub reasonForFail {
	my ($gene1, $gene2, $perID, $g1Cov, $g2Cov, $run) = @_;
	my ($status, @strings);
	my $g1 = Library::b73format($gene1); $g1 =~ s/_(T|C)\d+//;
  my $g2 = Library::b73format($gene2); $g2 =~ s/_(T|C)\d+//;
  my ( $genotype1, $genotype2 ) = split 'v', $run;

	if ( $perID =~ /NA/ ) {
  	$status = 'both';
	}
 	elsif ( $perID < 0.75 || ($g1Cov < 0.50 && $g2Cov < 0.50) ) {
  	$status = 'both';
 	}
 	elsif ( $g1Cov < 0.50 ) {
  	$status = 'g1';
 	}
 	elsif ( $g2Cov < 0.50 ) {
  	$status = 'g2';
 	}

	if ( $status eq 'both' ) {
		my $string1 = $hashA{$g1}{'rep_trans'} . "\t" . $hashB{$g2}{'chr'} . ":" . $hashB{$g2}{'start'} . "-" . $hashB{$g2}{'stop'} . "\t" . $genotype1 . "v" . $genotype2;
		my $string2 = $hashB{$g2}{'rep_trans'} . "\t" . $hashA{$g1}{'chr'} . ":" . $hashA{$g1}{'start'} . "-" . $hashA{$g1}{'stop'} . "\t" . $genotype2 . "v" . $genotype1;
		push @strings, $string1;
		push @strings, $string2;
	}
	elsif ( $status eq 'g1' ) {
		my $string1 = $hashA{$g1}{'rep_trans'} . "\t" . $hashB{$g2}{'chr'} . ":" . $hashB{$g2}{'start'} . "-" . $hashB{$g2}{'stop'} . "\t" . $genotype1 . "v" . $genotype2;
		push @strings, $string1;
	}
	elsif ( $status eq 'g2' ) {
		my $string1 =  $hashB{$g2}{'rep_trans'} . "\t" . $hashA{$g1}{'chr'} . ":" . $hashA{$g1}{'start'} . "-" . $hashA{$g1}{'stop'} . "\t" . $genotype2 . "v" . $genotype1;
		push @strings, $string1;
	}
	return(\@strings);
}

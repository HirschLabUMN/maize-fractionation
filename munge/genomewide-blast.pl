#!/usr/bin/perl -w
use strict;
use Library;
use Getopt::Std;
use Bio::Index::Fasta;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlastPlus;
use Storable;
use YAML::XS qw/LoadFile/;
##About: This script is used to blast genes that did not have a reciprocal match against the chromosome space
## broha006@ln0003 [/home/hirschc1/shared/projects/fractionation/munge] % perl genomewide-blast.pl -i /scratch.global/broha006/projects/frac/B73vB73_blsn-target-file.txt -o test.out -f /home/hirschc1/shared/projects/fractionation/data/assests/B73CDSFaIndex -b /home/maize/shared/databases/blast/Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa -c gene_region -r B73vB73

# perl genomewide-blast.pl -i /scratch.global/broha006/projects/frac/blastNA_targets.txt -o test.out -f /home/hirschc1/shared/projects/fractionation/data/assests/B73CDSFaIndex -b /home/maize/shared/databases/blast/Zea_mays/B73/Zea_mays.AGPv4.dna.toplevel.fa -c chr_region -r B73vB73

our ($opt_i, $opt_f, $opt_b, $opt_o, $opt_c, $opt_r);
getopts("i:o:f:b:c:r:");
my $config = '/home/hirschc1/shared/projects/fractionation/src/config.yml';
open ( my $infile, '<', $opt_i);
open ( my $outfile, '>', $opt_o) or die "counldn't open output\n";
print $outfile "header\n"; # print dummy line so snakemake doesn't choke looking for output

my $settings = LoadFile($config) or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";
my ( $query, $subject ) = split 'v', $opt_r;
$query =~ tr/A-Z/a-z/; $subject =~ tr/A-Z/a-z/;
my %queryHash = %{ retrieve($settings->{'db'}->{"$query"}->{'storable'}) } or die; # Get storable database
my %subjHash = %{ retrieve($settings->{'db'}->{"$subject"}->{'storable'}) } or die;
my $orderHashRef = &Library::fetchChrOrder($query, $config);

#Bio::Seq object fasta index
my $indexName = $opt_f;
my $fastaIndex = Bio::Index::Fasta->new($indexName);

<$infile>;
while ( my $line = <$infile> ) {
  chomp $line;
  my @line = split '\t', $line;
	my ( $gene, $loc ) = @line[0,1];
	next if ( ($line[2] ne $query || $line[4] ne $subject) && $opt_c eq 'chr_region' );
	my $results = &runBlast($fastaIndex, $gene, $opt_b);

	my @meta;
  if ( $opt_c eq 'gene_region' && $opt_r eq $line[2] ) {
    my ( $chr, $coord ) = split ':', $loc;
    my ( $start, $stop ) = split '-', $coord;
    $start = $start - 5000;  # we want to provide some additional buffer for the gene to align to
    $stop = $stop + 5000;    # deal with shifted annotation causing the pairwise alignment to fail.
		@meta = ( $chr, $start, $stop );
  }
  elsif ( $opt_c eq 'chr_region' ) {
		#	 $gene, $config, $query, $subject, $orderRef, $feature
		@meta = ( $gene, $loc, $config, $query, $subject, $orderHashRef, 'homeolog' );
  }
	else {
		print "unrecognized run mode, -c flag must be set to gene_region or chr_region\n\n";
		exit 1;
	}
	my $string = &parseBlast( \@meta, $gene, $$results );
	print $outfile "$string\n";
}
close $infile;
close $outfile;
exit 0;

sub runBlast {
	my ( $index, $gene, $db ) = @_;
	my $seq = $index->fetch($gene); # fetch gene sequence from fasta index

	my $factory = Bio::Tools::Run::StandAloneBlastPlus->new( -db_name => "$db" ) or die "couldn't get blast factory\n";
	# Generate a random string for the temporary output file name
	my $temp_file;
	my @chars = ("A".."Z", "a".."z");
	$temp_file .= $chars[rand @chars] for 1..8;

	# here we call blastn
	$factory->blastn(
							-outfile=> $temp_file,
							-query=> $seq,
							-outformat=> 5,
							-method_args=> [ '-task'=>'blastn', '-max_target_seqs'=>5, '-max_hsps'=>1, '-num_threads'=>16, '-evalue'=>1e-30 ]);
	$factory->cleanup;
	return \$temp_file;
}

sub parseBlast {
	my $metaInfo = $_[0];
	my $id = $_[1];
  my $in = new Bio::SearchIO(-format => 'blastxml', -file => $_[2] );

	my $resultString = "$id\tno_hits_found"; # the default is no hits found
  # Loop through results
  RESULTS: while( my $result = $in->next_result ) {
    # Get hits as long as they are defined
		if ( $result->num_hits < 1 ) {
			return $resultString;
		}

		while ( my $hit = $result->next_hit ) {
      my $hsp = $hit->next_hsp;
			my $status;

			if ( $opt_c eq 'chr_region' ) {
				next if $hit->name ne ${$metaInfo}[1];
				$status = &parseChrRegion( $metaInfo, $hit, $hsp );
			}
			else {
				$status = &parseGeneRegion( $metaInfo, $hit, $hsp )
			}

			# If we've found a PASS then exit loop
  		if ( $status eq 'PASS' ) {
				my $hitInfo = join(",", $hit->rank,$hit->query_length,$hit->frac_identical,$hit->frac_aligned_query);
				my $hspInfo = join(",", $hsp->length('total'),$hit->length_aln,$hsp->percent_identity,$hit->significance);
				my $locInfo = $hit->name . ":" . $hsp->hit->start . "-" . $hsp->hit->end;
				$resultString = join("\t", $id,$hitInfo,$hspInfo,$locInfo,$status);
				last RESULTS;
			}
		}
	}
	`rm $_[2]`; # get rid of temp alignment file.
	return $resultString;
}

sub parseGeneRegion {
	my ( $META, $HIT, $HSP ) = @_;
	# CHECK TO MAKRE SURE HSP IS IN THE RIGHT REGION
	my $hit_chr = $HIT->name;
	my ( $CHR, $START, $STOP ) = @$META[0,1,2];
	$hit_chr =~ s/chr0//; $hit_chr =~ s/chr10/10/; # parse to get rid of chr prefix
	if ( ($hit_chr eq $CHR) && ( ($HSP->hit->start > $START && $HSP->hit->start < $STOP) || ($HSP->hit->end > $START &&  $HSP->hit->end < $STOP) ) ) {
		return 'PASS';
	}
	else {
		return 'FAIL';
	}
}

sub parseChrRegion {
	my ( $META, $HIT, $HSP ) = @_;
	my ( $gene, $loc, $config, $query, $subject, $orderRef, $feature ) = @$META[0..$#$META];
	my $hit_chr = $HIT->name;
	$hit_chr =~ s/chr0//; $hit_chr =~ s/chr10/10/; # parse to get rid of chr prefix
	my $start = $HSP->hit->start;
	my $stop = $HSP->hit->end;

	# we need to check to see if the gene maps in a collinear regions. see subroutine
	# Zm00001d... ; config.yml ; b73 ; b73 ; b73orderRef ; duplicate
	my $maploc = &Library::getUpgetDown( $gene, $config, $query, $subject, $orderRef, $feature );
	if ( $maploc ) {
		my $leftBoundary = $maploc->{'left'}->{'start'} - 50000;
		my $rightBoundary = $maploc->{'right'}->{'stop'} + 50000;
		if ( ( $start > $leftBoundary && $start < $rightBoundary ) || ( $stop > $leftBoundary && $stop < $rightBoundary ) ) {
			return 'PASS';
		}
		else {
			return 'FAIL';
		}
	}
	else {
		return 'FAIL';
	}
}

# Ex. B73maize1 gene. Cognate- PH207 m1 ; Homoelog- B73 m2 ; Reciprocal- PH207m2
sub getFeature {
  my ( $ref ) = @_;
  my ( $queryGenotype, $querySubgenome, $subjGenotype, $subjSubgenome ) = @$ref[2..5];
  if ( $queryGenotype ne $subjGenotype ) {
    if ( $querySubgenome ne $subjSubgenome ) {
      return('reciprocal');
    }
    else {
      return('cognate');
    }
  }
  else {
    return('homeolog');
  }
}

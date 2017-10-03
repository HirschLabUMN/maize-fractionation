#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Printer;
use Library;

##About: Takes output from CoGE that has list of syntelogs between two genomes as input.
# Output is Sorghum \t Maize1_Gene \t Maize2_Gene \t Maize3_Gene
# NOTE: If two Zm genes are assigned syntelogs to the same Sorghum gene, they will be separated by comma.
# To separate later these I used grep "," foo.txt > duplicates.txt

# Run like: perl ~/families/scripts/frac/coge_filtering/filter_condensed_syntelog.pl -i condensed_syntelogs_wGEvo.txt -o ~/families/analysis/frac/coge/B73v2vSb/tmp.txt

my $usage = "\n\n$0 -i <final_syntenic_set_wGEvo_links.txt> -o <B73vSb_tmp.txt> -c <Os||Sb>\n\n";

our ($opt_i, $opt_o, $opt_c, $opt_h);
getopts("i:o:c:h") or die "";

if ( (!(defined $opt_i)) || (!(defined $opt_c)) || (!(defined $opt_o)) || (defined $opt_h)) {die "$usage";}

open (my $infile, '<', $opt_i) or die "Cannot open $opt_i $!\n\n";
open (my $outfile, '>', $opt_o) or die "Cannot open $opt_o $!\n\n";

my %hash;
#Create a header
print $outfile "Anchor\tMaize1\tMaize2\n";

my $ref;
if ( $opt_c =~ /Sb/ ) {
	$ref = 'Sb';
}
else {
	$ref = 'Os';
}

## READ-THRU CONDENSED SYNTELOG INPUT ##
while (my $line = <$infile> ) {
	chomp $line;

	#Skip the final output of quota align that gives coverage stats
	if ($line =~ /coverage/) {
		next;
	}

	#Capture clusters (these are ignored here)
	if ($line =~ /^#/) {
		next;
	}

	# Split lines to get columns
	my (undef, $query, undef, undef, undef, $subj, undef, undef, undef, undef, undef) = split("\t", $line);
  # Split meta from Anchor and Query to get further info
	my ( $anchrChr, $anchrStart, $anchrEnd, $anchrGene, $anchrStrand );
	my ( $queryChr, $queryStart, $queryEnd, $queryGene, $queryStrand );

	if ( $query =~ /Zm0000/ ) {
		($queryChr, $queryStart, $queryEnd, $queryGene, $queryStrand, undef, undef, undef, undef) = split('\|\|', $query);
		($anchrChr, $anchrStart, $anchrEnd, $anchrGene, $anchrStrand, undef, undef, undef, undef) = split('\|\|', $subj);
	}
	else {
		($queryChr, $queryStart, $queryEnd, $queryGene, $queryStrand, undef, undef, undef, undef) = split('\|\|', $subj);
		($anchrChr, $anchrStart, $anchrEnd, $anchrGene, $anchrStrand, undef, undef, undef, undef) = split('\|\|', $query);
	}

	my $block_id = $anchrChr . "_" . $queryChr; # Create an id to check against hash for chromsome asssignment
	# Fetch assignment
	my $assignment = Library::syntenicBlockAssignment($anchrChr, $queryChr, $ref);

	# print "\n$line\n$block_id\t$assignment\n\n";
	# die;

	# Have you seen this anchor already ?
  ## If already seen, update designations
  if ( exists($hash{$anchrGene}) ) {
    # Is the assignment maize1 ?
    if ( $assignment eq 'maize1' ) {
      # Is this the first time seeing a gene for the subgenome ?
      if ( $hash{$anchrGene}{'m1'} eq 'NA' ) {
        $hash{$anchrGene}{'m1'} = $queryGene;
      }
      # If a gene already seen add it as a duplicate...
      else {
        $hash{$anchrGene}{'m1'} = join(',', $hash{$anchrGene}{'m1'}, $queryGene);
      }
    }
    elsif ( $assignment eq 'maize2' ) {
      # Is this the first time seeing a gene for the subgenome ?
      if ( $hash{$anchrGene}{'m2'} eq 'NA' ) {
        $hash{$anchrGene}{'m2'} = $queryGene;
      }
      # If a gene already seen add it as a duplicate...
      else {
        $hash{$anchrGene}{'m2'} = join(',', $hash{$anchrGene}{'m2'}, $queryGene);
      }
    }
    else {
      next;
    }
  }
  ## This is the first time seeing the anchor, initialize hash values ...
  else {
    if ( $assignment eq 'maize1' ) {
      $hash{$anchrGene}{'m1'} = $queryGene;
      $hash{$anchrGene}{'m2'} = 'NA';
    }
    elsif ( $assignment eq 'maize2') {
      $hash{$anchrGene}{'m2'} = $queryGene;
      $hash{$anchrGene}{'m1'} = 'NA';
    }
    else {
      next;
    }
  }
}

## PRINT TO OUTPUT ##

my ( $count, $bothretained, $m1only, $m2only );
foreach my $key (keys %hash) {
	print $outfile "$key\t$hash{$key}{'m1'}\t$hash{$key}{'m2'}\n";
	++$count;
	if ( $hash{$key}{'m1'} ne 'NA' && $hash{$key}{'m2'} ne 'NA' ) {
		++$bothretained;
	}
	elsif ( $hash{$key}{'m1'} ne 'NA' && $hash{$key}{'m2'} eq 'NA' ) {
		++$m1only;
	}
	elsif ( $hash{$key}{'m1'} eq 'NA' && $hash{$key}{'m2'} ne 'NA' ) {
		++$m2only;
	}
	else {
		next;
	}
}

print "Fetched syntelogs for $opt_c\n#Total: $count Both_Retained: $bothretained Maize1_Only: $m1only Maize2_Only: $m2only\n";

close $infile;
close $outfile;
exit;

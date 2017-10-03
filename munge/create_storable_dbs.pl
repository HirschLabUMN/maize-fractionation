#!/usr/bin/perl -w
use strict;
use Storable;
use 5.010;
use Data::Printer;
use YAML::XS qw/LoadFile/;
use Getopt::Std;
use Sort::Naturally;
#ABB: 10.28.16
##ABOUT: Script will insert data into the storable database. The script is coded to parse Phytozome
# and Gramene GFF syntax (why can't it just be standard?).
###RUNLIKE: perl munge/create_storable_dbs.pl -i ph207 -m Phytozome

our ($opt_i,$opt_m,$opt_h);
getopts("i:m:h");

my $usage = "\n\n$0 -i <ph207,b73,sb,os> -m <Phytozome,Gramene>\n\n";

if (defined($opt_h) ) {
	die $usage;
}

my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";
my $in;
if ( $opt_i eq 'B73' || $opt_i eq 'b73' ) {
	open ( $in, '<', $settings->{'genome'}->{'b73'}->{'gff'} ) or die;
}
elsif ( $opt_i eq 'PH207' || $opt_i eq 'ph207' ) {
	open ( $in, '<', $settings->{'genome'}->{'ph207'}->{'gff'} ) or die;
}
elsif ( $opt_i eq 'Sb' || $opt_i eq 'sb' ) {
	open ( $in, '<', $settings->{'genome'}->{'sb'}->{'gff'} ) or die;
}
elsif ( $opt_i eq 'Os' ||  $opt_i eq 'os' ) {
	open ( $in, '<', $settings->{'genome'}->{'os'}->{'gff'} ) or die;
}
else {
	die "unknown option supplied: should be B73, PH207, Sb, or Os\n";
}

my ( %store, %cds_len_hash, %len_hash, %hier, %geneOrder);

while ( my $line = <$in>  ) {
	chomp $line;

  if ( $line =~ /^#/ ) {
    next;
  }

  my @fields = split '\t', $line;
	my @meta = split ';', $fields[8];
	my %metaHash;

	# $hash{Name} => Gene_X
	foreach my $entry ( @meta ) {
		my ( $id, $data ) = split '=', $entry;
		$metaHash{$id} = $data;
	}

	# clean-up chromosome name
	my $chr = $fields[0]; $chr =~ s/chr0//; $chr =~ s/Chr0//; $chr =~ s/chr//; $chr =~ s/Chr//;

  if ( $fields[2] eq 'gene' ) {
		# Gene_X = $hash{Name}
		my $name;
		if ( $opt_m eq 'Gramene' ) {
			$name = $metaHash{'ID'};
		}
		elsif ( $opt_m eq 'Phytozome' ) {
			$name = $metaHash{'Name'};
			my $id = $metaHash{'ID'};
			$hier{$id} = $name; #Zm00008a000001.v1.1 => Zm00008a000001
		}

		# what to store
		$store{$name}{'chr'} = $chr;
  	$store{$name}{'start'} = $fields[3];
  	$store{$name}{'stop'} = $fields[4];
  	$store{$name}{'rep_trans'} = 'foo';
  	$store{$name}{'length'} = '0';
		$store{$name}{'cds_length'} = '0';
		$store{$name}{'synteny'} = 'NA';
		$store{$name}{'order'} = 'NA';
		$store{$name}{'homeolog'} = 'NA';
		$store{$name}{'cognate'} = 'NA';
		$store{$name}{'reciprocal'} = 'NA';
		# store gene order, it two genes have the same start position then add one to the redundant gene
		if ( exists($geneOrder{$chr}{$fields[3]}) ) {
			$geneOrder{$chr}{$fields[3]+1} = $name;
		}
		else {
		$geneOrder{$chr}{$fields[3]} = $name;
		}
	}
	elsif ( $fields[2] eq 'mRNA' || $fields[2] eq 'transcript' || $fields[2] eq 'lincRNA' || $fields[2] eq 'miRNA' ) {
		if ( $opt_m eq 'Gramene') {
			my $id = $metaHash{'ID'};
			my $parent = $metaHash{'Parent'};
			$hier{$id} = $parent;
		}
		elsif ( $opt_m eq 'Phytozome' ) {
			my $id = $metaHash{'ID'}; # Zm00008a000001_T01.v1.1
			my $parent = $metaHash{'Parent'}; # Zm00008a000001.v1.1
			my $mRNA = $metaHash{'Name'}; # Zm00008a000001_T01
			$hier{$id} = $mRNA; # Zm00008a000001_T01.v1.1 => Zm00008a000001_T01
			$hier{$mRNA} = $parent; # Zm00008a000001_T01 => Zm00008a000001.v1.1
		}
	}
	elsif ( $fields[2] eq 'exon' ) {
		my ($gene, $grandparent);
		if ( $opt_m eq 'Gramene' ) {
			$grandparent = $metaHash{'Parent'};
			$gene = $hier{$grandparent};
		}
		elsif ( $opt_m eq 'Phytozome') {
			my $id = $metaHash{'ID'}; # Zm00008a000001_T01.v1.1.exon.1
	  	my $parent = $metaHash{'Parent'}; # Zm00008a000001_T01.v1.1
			$grandparent = $hier{$parent}; # Zm00008a000001_T01
			my $name = $hier{$grandparent}; # Zm00008a000001.v1.1
			$gene = $hier{$name}; # Zm00008a000001
		}
		my $exon_len = abs($fields[4] - $fields[3]) + 1;

		if ( exists($len_hash{$gene}{$grandparent}) ) {
			# add to the length
			$len_hash{$gene}{$grandparent} = $len_hash{$gene}{$grandparent} += $exon_len;
		}
		else {
			# otherwise, init
			$len_hash{$gene}{$grandparent} = $exon_len;
		}
	}
	elsif ( $fields[2] eq 'CDS' ) {
		my ($gene, $grandparent);
		if ( $opt_m eq 'Gramene' ) {
			$grandparent = $metaHash{'Parent'};
			$gene = $hier{$grandparent};
		}
		elsif ( $opt_m eq 'Phytozome') {
			my $id = $metaHash{'ID'}; # Zm00008a000001_T01.v1.1.exon.1
			my $parent = $metaHash{'Parent'}; # Zm00008a000001_T01.v1.1
			$grandparent = $hier{$parent}; # Zm00008a000001_T01
			my $name = $hier{$grandparent}; # Zm00008a000001.v1.1
			$gene = $hier{$name}; # Zm00008a000001
		}
		my $cds_len = abs($fields[4] - $fields[3]) + 1;

		if ( exists($cds_len_hash{$gene}{$grandparent}) ) {
			# add to the length
			$cds_len_hash{$gene}{$grandparent} = $cds_len_hash{$gene}{$grandparent} += $cds_len;
		}
		else {
			# otherwise, init
			$cds_len_hash{$gene}{$grandparent} = $cds_len;
		}
	}
	else {
		next;
	}
#p %metaHash;
}
close $in;

# Get the representative transcript
foreach my $key ( keys %len_hash ) {
  foreach my $key2 ( keys %{ $len_hash{$key} } ) {
		if (exists($store{$key})) {
    	if ($len_hash{$key}{$key2} > $store{$key}{'length'}) {
      	$store{$key}{'length'} = $len_hash{$key}{$key2};
      	$store{$key}{'rep_trans'} = $key2;
			}
    }
  }
}

# Get the representative transcript
foreach my $key ( keys %cds_len_hash ) {
  foreach my $key2 ( keys %{ $cds_len_hash{$key} } ) {
		if (exists($store{$key})) {
    	if ($cds_len_hash{$key}{$key2} > $store{$key}{'cds_length'}) {
      	$store{$key}{'cds_length'} = $cds_len_hash{$key}{$key2};
      	$store{$key}{'rep_cds'} = $key2;
			}
    }
  }
}

# Get gene order
foreach my $id (nsort keys %geneOrder ) {
  my $counter = 0;
  foreach my $id2 ( nsort keys %{ $geneOrder{$id} } ) {
    my $gene = $geneOrder{$id}{$id2};
		$store{$gene}{'order'} = $counter;
    ++$counter;
  }
}

my $file = '/home/hirschc1/shared/projects/fractionation/config/' . "$opt_i" . '_gff.storable';

store(\%store, $file);
exit 0;

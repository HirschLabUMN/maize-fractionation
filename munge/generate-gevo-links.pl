#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Printer;

# https://genomevolution.org/CoGe/GEvo.pl?prog=blastz;iw=1000;fh=20;padding=2;hsp_top=1;nt=1;cbc=0;spike_len=15;ca=1;skip_feat_overlap=1;skip_hsp_overlap=1;hs=0;bzW=8;bzK=3000;bzO=400;bzE=30;
#
# accn1=LOC_Os01g08760;fid1=303988874;dsid1=66093;dsgid1=16888;chr1=Chr1;dr1up=10000;dr1down=10000;ref1=1;
#
# accn2=Zm00001d034848;fid2=1088847885;dsid2=105450;dsgid2=34067;chr2=1;dr2up=10000;dr2down=10000;ref2=1;accn3=Zm00008a005802;fid3=1063955599;dsid3=105227;dsgid3=33916;chr3=1;dr3up=10000;dr3down=10000;ref3=1;num_seqs=3;hsp_overlap_limit=0;hsp_size_limit=0

our ($opt_i, $opt_o);
getopts("i:o:");
open ( my $infile, '<', $opt_i ) or die "$! couldn't open infile\n";
open ( my $outfile, '>', $opt_o ) or die "$! couldn't open outfile\n";

my $prefix = 'https://genomevolution.org/CoGe/GEvo.pl?prog=blastz;iw=1000;fh=20;padding=2;hsp_top=1;nt=1;cbc=0;spike_len=15;ca=1;skip_feat_overlap=1;skip_hsp_overlap=1;hs=0;bzW=8;bzK=3000;bzO=400;bzE=30;';
my $suffix = 'hsp_overlap_limit=0;hsp_size_limit=0;autogo=1';
my $sb_dsid = '99077';
my $sb_dsgid = '28853';
my $os_dsid = '66093;';
my $os_dsgid = '16888;';
my $b73_dsid = '105450;';
my $b73_dsgid = '34067;';
my $ph207_dsid = '105227;';
my $ph207_dsgid = '33916;';

<$infile>;
while ( my $line = <$infile> ) {
	chomp $line;
	my @fields = split '\t', $line;
#	my @genes = shift @fields;
	my @genes;
	push @genes, grep { $_ =~ /(Zm|Sobic|LOC)/ } @fields;
	my $stringseq = '';
	for ( my $i = 0; $i < @genes; $i++ ) {
		my $count = $i + 1;
		my $dbid;
		my $gene = $genes[$i];
		my $string = 'accn' . $count . "=$gene;";

		if ( $gene =~ /Sobic/ ) {
			$dbid = "dsid" . "$count=" . $sb_dsid . "dsgid" . "$count=" . $sb_dsgid;
		}
		elsif ( $gene =~ /LOC/ ) {
			$dbid = "dsid" . "$count=" . $os_dsid . "dsgid" . "$count=" . $os_dsgid;
		}
		elsif ( $gene =~ /Zm00001d/ ) {
			$dbid = "dsid" . "$count=" . $b73_dsid . "dsgid" . "$count=" . $b73_dsgid;
		}
		elsif ( $gene =~ /Zm00008a/ ) {
			$dbid = "dsid" . "$count=" . $ph207_dsid . "dsgid" . "$count=" . $ph207_dsgid;
		}
		else {
			die "unrecognized gene: $gene\n\n";
		}

		my $buffer = "dr$count" . "up=10000;dr$count" . "down=10000;ref$count" . "=1;";
		$stringseq = $stringseq . $string . $dbid . $buffer;
	}
	my $geneCount = @genes;
	my $link = $prefix . $stringseq . "num_seqs=$geneCount;" . $suffix;
	print $outfile "$line\t$link\n";
}
close $infile;
close $outfile;

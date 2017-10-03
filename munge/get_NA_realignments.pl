#!/usr/bin/perl -w
use strict;
use feature qw/say/;
use Getopt::Long;
use Storable;
use Library;
use YAML::XS qw/LoadFile/;

my $usage = "USAGE: $0 [--query --subj --out]\n\n";

my ($input, $query, $subj, $out, $infile);
GetOptions(
  'input=s' => \$input,
  'query=s' => \$query,
  'subj=s' => \$subj,
  'out=s' => \$out,
) or die "$usage\n";

open $infile, '<', $input or die "error openening $input\n\n";
open my $outfile, '>', $out or die "error opening $out $!\n\n";

my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";
my $in;

# Get storable database
my %queryHash = %{ retrieve($settings->{'db'}->{"$query"}->{'storable'}) } or die;
my %subjHash = %{ retrieve($settings->{'db'}->{"$subj"}->{'storable'}) } or die;
my %sbHash = %{ retrieve($settings->{'db'}->{'sb'}->{'storable'}) } or die;
my %osHash = %{ retrieve($settings->{'db'}->{'os'}->{'storable'}) } or die;

<$infile>;
while ( my $line = <$infile> ) {
  chomp $line;

  if ( $line =~ /,/ ) {
    next;
  }

  my ( $anchr, $b1, $b2, $p1, $p2, undef ) = split '\t', $line;

  ###############################
  ## GET THE ANCHOR CHROMOSOME ##
  ###############################

  my ( $sp, $chr, $target,);
  if ( $anchr =~ /Sb/ || $anchr =~ /Sobic/ ) {
    $sp = 'Sb';
    $chr = $sbHash{$anchr}{'chr'};
  }
  elsif ( $anchr =~ /Os/ ) {
    $sp = 'Os';
    $chr = $osHash{$anchr}{'chr'};
  }
  else {
    die "unrecognized anchor. this must be Sorghum or Rice\n";
  }

  ##############################
  ## GET TARGET CHR AND PRINT ##
  ##############################

  # CHECK FOR EVERY POSSIBLE SCERNARIO INVOLVING NA
  ## 1
  if ( $b1 ne 'NA' && $b2 eq 'NA' && $p1 eq 'NA' && $p2 eq 'NA' ) {
    # CALL SUBROUTINE TO GET TARGET CHR THAT GENE SHOULD ALIGN TO
    if ( $b1 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #b1 to p1
      say $outfile "$queryHash{$b1}{'rep_cds'}\t$target\tb73\tmaize1\tph207\tmaize1";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #b1 to p2
      say $outfile "$queryHash{$b1}{'rep_cds'}\t$target\tb73\tmaize1\tph207\tmaize2";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #b1 to b2
      say $outfile "$queryHash{$b1}{'rep_cds'}\t$target\tb73\tmaize1\tb73\tmaize2";
    }
  }
  ## 2
  elsif ( $b1 eq 'NA' && $b2 ne 'NA' && $p1 eq 'NA' && $p2 eq 'NA' ) {
    if ( $b2 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #b2 to p1
      say $outfile "$queryHash{$b2}{'rep_cds'}\t$target\tb73\tmaize2\tph207\tmaize1";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #b2 to p2
      say $outfile "$queryHash{$b2}{'rep_cds'}\t$target\tb73\tmaize2\tph207\tmaize2";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #b2 to b1
      say $outfile "$queryHash{$b2}{'rep_cds'}\t$target\tb73\tmaize2\tb73\tmaize1";
    }
  }
  ## 3
  elsif ( $b1 eq 'NA' && $b2 eq 'NA' && $p1 ne 'NA' && $p2 eq 'NA' ) {
    if ( $p1 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #p1 to b1
      say $outfile "$subjHash{$p1}{'rep_cds'}\t$target\tph207\tmaize1\tb73\tmaize1";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #p1 to b2
      say $outfile "$subjHash{$p1}{'rep_cds'}\t$target\tph207\tmaize1\tb73\tmaize2";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #p1 to p2
      say $outfile "$subjHash{$p1}{'rep_cds'}\t$target\tph207\tmaize1\tph207\tmaize2";
    }
  }
  ## 4
  elsif ( $b1 eq 'NA' && $b2 eq 'NA' && $p1 eq 'NA' && $p2 ne 'NA' ) {
    if ( $p2 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #p2 to b1
      say $outfile "$subjHash{$p2}{'rep_cds'}\t$target\tph207\tmaize2\tb73\tmaize1";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #p2 to b2
      say $outfile "$subjHash{$p2}{'rep_cds'}\t$target\tph207\tmaize2\tb73\tmaize2";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #p2 to p1
      say $outfile "$subjHash{$p2}{'rep_cds'}\t$target\tph207\tmaize2\tph207\tmaize1";
    }
  }
  ## 5
  elsif ( $b1 ne 'NA' && $b2 eq 'NA' && $p1 ne 'NA' && $p2 eq 'NA'  ) {
    if ( $b1 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #b1 to p2
      say $outfile "$queryHash{$b1}{'rep_cds'}\t$target\tb73\tmaize1\tph207\tmaize2";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #b1 to b2
      say $outfile "$queryHash{$b1}{'rep_cds'}\t$target\tb73\tmaize1\tb73\tmaize2";
    }
    if ( $p1 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #p1 to p2
      say $outfile "$subjHash{$p1}{'rep_cds'}\t$target\tph207\tmaize1\tph207\tmaize2";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #p1 to b2
      say $outfile "$subjHash{$p1}{'rep_cds'}\t$target\tph207\tmaize1\tb73\tmaize2";
    }
  }
  ## 6
  elsif ( $b1 eq 'NA' && $b2 ne 'NA' && $p1 eq 'NA' && $p2 ne 'NA'  ) {
    if ( $b2 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #b2 to b1
      say $outfile "$queryHash{$b2}{'rep_cds'}\t$target\tb73\tmaize2\tb73\tmaize1";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #b2 to p1
      say $outfile "$queryHash{$b2}{'rep_cds'}\t$target\tb73\tmaize2\tph207\tmaize1";
    }
    if ( $p2 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #p2 to p1
      say $outfile "$subjHash{$p2}{'rep_cds'}\t$target\tph207\tmaize2\tph207\tmaize1";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #p2 to b1
      say $outfile "$subjHash{$p2}{'rep_cds'}\t$target\tph207\tmaize2\tb73\tmaize1";
    }
  }
  ## 7
  elsif ( $b1 ne 'NA' && $b2 ne 'NA' && $p1 eq 'NA' && $p2 eq 'NA' ) {
    if ( $b1 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #b1 to p1
      say $outfile "$queryHash{$b1}{'rep_cds'}\t$target\tb73\tmaize1\tph207\tmaize1";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #b1 to p2
      say $outfile "$queryHash{$b1}{'rep_cds'}\t$target\tb73\tmaize1\tph207\tmaize2";
    }
    if ( $b2 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #b2 to p2
      say $outfile "$queryHash{$b2}{'rep_cds'}\t$target\tb73\tmaize2\tph207\tmaize2";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #b2 to p1
      say $outfile "$queryHash{$b2}{'rep_cds'}\t$target\tb73\tmaize2\tph207\tmaize1";
    }
  }
  ## 8
  elsif ( $b1 eq 'NA' && $b2 eq 'NA' && $p1 ne 'NA' && $p2 ne 'NA' )  {
    if ( $p1 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #p1 to b1
      say $outfile "$subjHash{$p1}{'rep_cds'}\t$target\tph207\tmaize1\tb73\tmaize1";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #p1 to b2
      say $outfile "$subjHash{$p1}{'rep_cds'}\t$target\tph207\tmaize1\tb73\tmaize2";
    }
    if ( $p2 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #p2 to b1
      say $outfile "$subjHash{$p2}{'rep_cds'}\t$target\tph207\tmaize2\tb73\tmaize1";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #p2 to b2
      say $outfile "$subjHash{$p2}{'rep_cds'}\t$target\tph207\tmaize2\tb73\tmaize2";
    }
  }
  ## 9
  elsif ( $b1 ne 'NA' && $b2 eq 'NA' && $p1 eq 'NA' && $p2 ne 'NA' ) {
    if ( $b1 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #b1 to p1
      say $outfile "$queryHash{$b1}{'rep_cds'}\t$target\tb73\tmaize1\tph207\tmaize1";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #b1 to b2
      say $outfile "$queryHash{$b1}{'rep_cds'}\t$target\tb73\tmaize1\tb73\tmaize2";
    }
    if ( $p2 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #p2 to b2
      say $outfile "$subjHash{$p2}{'rep_cds'}\t$target\tph207\tmaize2\tb73\tmaize2";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #p2 to p1
      say $outfile "$subjHash{$p2}{'rep_cds'}\t$target\tph207\tmaize2\tph207\tmaize1";
    }
  }
  ## 10
  elsif ( $b1 eq 'NA' && $b2 ne 'NA' && $p1 ne 'NA' && $p2 eq 'NA' ) {
    if ( $p1 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #p1 to b1
      say $outfile "$subjHash{$p1}{'rep_cds'}\t$target\tph207\tmaize1\tb73\tmaize1";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #p1 to p2
      say $outfile "$subjHash{$p1}{'rep_cds'}\t$target\tph207\tmaize1\tph207\tmaize2";
    }
    if ( $b2 !~ /:/ ) {
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #b2 to b1
      say $outfile "$queryHash{$b2}{'rep_cds'}\t$target\tb73\tmaize2\tb73\tmaize1";
      $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #b2 to p2
      say $outfile "$queryHash{$b2}{'rep_cds'}\t$target\tb73\tmaize2\tph207\tmaize2";
    }
  }
  ## 11
  elsif ( $b1 ne 'NA' && $b2 ne 'NA' && $p1 eq 'NA' && $p2 ne 'NA' ) {
    $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #b1 to p1
    say $outfile "$queryHash{$b1}{'rep_cds'}\t$target\tb73\tmaize1\tph207\tmaize1" unless $b1 =~ /:/;
    $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #p2 to p1
    say $outfile "$subjHash{$p2}{'rep_cds'}\t$target\tph207\tmaize2\tph207\tmaize1" unless $p2 =~ /:/;
  }
  ## 12
  elsif ( $b1 ne 'NA' && $b2 ne 'NA' && $p1 ne 'NA' && $p2 eq 'NA' ) {
    $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #b2 to p2
    say $outfile "$queryHash{$b2}{'rep_cds'}\t$target\tb73\tmaize2\tph207\tmaize2" unless $b2 =~ /:/;
    $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #p1 to p2
    say $outfile "$subjHash{$p1}{'rep_cds'}\t$target\tph207\tmaize1\tph207\tmaize2" unless $p1 =~ /:/;
  }
  ## 13
  elsif ( $b1 eq 'NA' && $b2 ne 'NA' && $p1 ne 'NA' && $p2 ne 'NA' ) {
    $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #p1 to b1
    say $outfile "$subjHash{$p1}{'rep_cds'}\t$target\tph207\tmaize1\tb73\tmaize1" unless $p1 =~ /:/;
    $target = Library::syntenicBlockTargets($sp, $chr, 'maize1'); #b2 to b1
    say $outfile "$queryHash{$b2}{'rep_cds'}\t$target\tb73\tmaize2\tb73\tmaize1" unless $b2 =~ /:/;
  }
  ## 14
  elsif ( $b1 ne 'NA' && $b2 eq 'NA' && $p1 ne 'NA' && $p2 ne 'NA' ) {
    $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #p2 to b2
    say $outfile "$subjHash{$p2}{'rep_cds'}\t$target\tph207\tmaize2\tb73\tmaize2" unless $p2 =~ /:/;
    $target = Library::syntenicBlockTargets($sp, $chr, 'maize2'); #b1 to b2
    say $outfile "$queryHash{$b1}{'rep_cds'}\t$target\tb73\tmaize1\tb73\tmaize2" unless $b1 =~ /:/;
  }
  else {
    next;
  }
}
close $infile;
close $outfile;

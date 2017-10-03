#!/usr/bin/perl  -w
use strict;
use Getopt::Std;
use Storable;
use Library;
use YAML::XS qw/LoadFile/;
use Data::Printer;

our ($opt_i,$opt_m,$opt_o);
getopts("i:m:o:");

my $settings = LoadFile('/home/hirschc1/shared/projects/fractionation/src/config.yml') or die "\nCouldn't load config file -- path should be hardcoded in script\n\n";

my %bStore = %{ retrieve($settings->{'db'}->{'b73'}->{'storable'}) } or die;
my %pStore = %{ retrieve($settings->{'db'}->{'ph207'}->{'storable'}) } or die;
my %sStore = %{ retrieve($settings->{'db'}->{'sb'}->{'storable'}) } or die;
my %oStore = %{ retrieve($settings->{'db'}->{'os'}->{'storable'}) } or die;

# Any gene in tandem duplicate cluster => list of rep maize genes => ortholog for each rep
my ($assignRef, $href) = &populateHash($opt_m);

open ( my $infile, '<', $opt_i ) or die;
open ( my $outfile, '>', $opt_o ) or die;

# READ THROUGH TANDEM DUPLICATE FILE
while ( my $line = <$infile> ) {
  chomp $line;
  my @fields = split '\t', $line;

  # Only care about the line if assignment made was actually syntenic
	# ( some of these assignments come from OrthoFinder which means not every assignment will be syntenic)
  if ( $fields[5] eq 'PASS' ) {
    my $ortholog = $fields[3];    # get some basic variables
    my $repGene = $fields[2];
    $ortholog =~ s/\.\d$//;

    my $subAssign = &getAssign($ortholog, $repGene);  # get subgenome assignment

    # Check to see if an assignment already existed for this ortholog/assignment combination
    if ( exists($assignRef->{$ortholog}->{$subAssign}) ) {
      my $previousAssignment = $assignRef->{$ortholog}->{$subAssign};   # get the previous maize gene
      # If there was already an assignment for this gene then append this as a duplicate to be sorted out with tblastx
      if ( $previousAssignment =~ /Zm/ ) {
        $assignRef->{$ortholog}->{$subAssign} = $previousAssignment . "," . $repGene unless $previousAssignment eq $repGene;
      }
      else {
        $assignRef->{$ortholog}->{$subAssign} = $repGene;
      }
    }
    else {          # If we haven't previously seen this ortholog:
      # Check to see if this is a sorghum ortholog replacing a rice ortholog
      # If so, then we need to convert hash entries associated with the rice ortholog to be with the sorghum ortholog
      if ( exists($href->{$repGene}) && ($href->{$repGene}->{'anchor'} !~ /Sobic/) ) {
        my $oldOrtholog = $href->{$repGene}->{'anchor'};
        $assignRef->{$ortholog}->{'bm1'} = delete $assignRef->{$oldOrtholog}->{'bm1'};
        $assignRef->{$ortholog}->{'bm2'} = delete $assignRef->{$oldOrtholog}->{'bm2'};
        $assignRef->{$ortholog}->{'pm1'} = delete $assignRef->{$oldOrtholog}->{'pm1'};
        $assignRef->{$ortholog}->{'pm2'} = delete $assignRef->{$oldOrtholog}->{'pm2'};
        delete $assignRef->{$oldOrtholog};
        $href->{ $assignRef->{$ortholog}->{'bm1'} }->{'anchor'} = $ortholog unless $assignRef->{$ortholog}->{'bm1'} eq 'NA';
        $href->{ $assignRef->{$ortholog}->{'bm2'} }->{'anchor'} = $ortholog unless $assignRef->{$ortholog}->{'bm2'} eq 'NA';
        $href->{ $assignRef->{$ortholog}->{'pm1'} }->{'anchor'} = $ortholog unless $assignRef->{$ortholog}->{'pm1'} eq 'NA';
        $href->{ $assignRef->{$ortholog}->{'pm2'} }->{'anchor'} = $ortholog unless $assignRef->{$ortholog}->{'pm2'} eq 'NA';
      }
      else {
        # If not replacing a rice ortholog then just define a new entry
        $assignRef->{$ortholog}->{'bm1'} = 'NA';
        $assignRef->{$ortholog}->{'bm2'} = 'NA';
        $assignRef->{$ortholog}->{'pm1'} = 'NA';
        $assignRef->{$ortholog}->{'pm2'} = 'NA';
        # Replace the NA at the current assignment
        $assignRef->{$ortholog}->{$subAssign} = $repGene;
      }
    }
  }
}

# Now we can ahead and print everything out
foreach my $key ( keys %$assignRef ) {
    print $outfile "$key";
  foreach my $sub ( qw/ bm1 bm2 pm1 pm2 / ) {
    if (exists($assignRef->{$key}->{$sub}) ) {
      print $outfile "\t$assignRef->{$key}->{$sub}";
    }
    else {
      print $outfile "\tNA";
    }
  }
  print $outfile "\n";
}
close $infile;
close $outfile;

# Put the working list of syntenic assignments into a hash
sub populateHash {
  my ( $handle ) = @_;
  open ( my $master, '<', $handle ) or die;
  my ( %hash, %assign );

  <$master>; # skip the header
  while ( my $line = <$master> ) {
    chomp $line;
    my @fields = split '\t', $line; # put each column into array
    my $anchr = shift @fields; # first column is always ancestral gene
    $assign{$anchr}{'bm1'} = $fields[0]; # second column is B73 maize1
    $hash{$fields[0]}{'anchor'} = $anchr unless $fields[0] eq 'NA';
    $hash{$fields[0]}{'subgenome'} = 'bm1' unless $fields[0] eq 'NA';
		# Repeat above for B73 maize2 and PH207 m1/m2...
    $assign{$anchr}{'bm2'} = $fields[1];
    $hash{$fields[1]}{'anchor'} = $anchr unless $fields[1] eq 'NA';
    $hash{$fields[1]}{'subgenome'} = 'bm2' unless $fields[1] eq 'NA';
    $assign{$anchr}{'pm1'} = $fields[2];
    $hash{$fields[2]}{'anchor'} = $anchr unless $fields[2] eq 'NA';
    $hash{$fields[2]}{'subgenome'} = 'pm1' unless $fields[2] eq 'NA';
    $assign{$anchr}{'pm2'}= $fields[3];
    $hash{$fields[3]}{'anchor'}  = $anchr unless $fields[3] eq 'NA';
    $hash{$fields[3]}{'subgenome'} = 'pm2' unless $fields[3] eq 'NA';
  }
	# exit and return hashes
  close $master;
  return \%assign, \%hash;
}

# To know which subgenome a gene belongs to, you need to know the sorghum and maize chromsomes. This routine will fetch these and return subgenome.
sub getAssign {
  my ( $olog, $rep ) = @_;
	# Get the sorghum chromosome
  my $anchrChr = $sStore{$olog}{'chr'};

  my ( $maizeChr, $maizeSp );

	# Is this a B73 gene or PH207 gene ?
  if ( $rep =~ /Zm00001d/ ) {
    $maizeChr = $bStore{$rep}{'chr'};
    $maizeSp = 'b';
  }
  else {
    $maizeChr = $pStore{$rep}{'chr'};
    $maizeSp = 'p';
  }

	# Call library to retrieve assignment for chromosome combination
  my $assignment = &Library::syntenicBlockAssignment($anchrChr, $maizeChr, 'Sb');
	# Reformat and return
  $assignment =~ s/maize1/m1/; $assignment =~ s/maize2/m2/;
  my $parsedAssignment = $maizeSp . $assignment; # bm1, pm1, etc.
  return $parsedAssignment;
}

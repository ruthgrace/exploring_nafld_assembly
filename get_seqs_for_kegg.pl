#!/usr/bin/env perl -w
#use strict;
my $help = "\tFirst argument is a path to the annotated counts table with zero features removed\n
Second argument is the assembled reference library\n
Third argument is the constructed reference library\n
Forth argument is the DNA output file (fasta format)\n
Fifth argument is the amino acid output file (.faa format)\n";
print $help if !$ARGV[0];exit if !$ARGV[0];
print $help if !$ARGV[1];exit if !$ARGV[1];
print $help if !$ARGV[2];exit if !$ARGV[2];
print $help if !$ARGV[3];exit if !$ARGV[3];
print $help if !$ARGV[4];exit if !$ARGV[4];

my $counts = $ARGV[0];
my $assembled = $ARGV[1];
my $constructed = $ARGV[2];
my $fa_out = $ARGV[3];
my $faa_out = $ARGV[4];

# read in all assembled refseqs
my %seqs;
my $id = "";
my $seq = "";
open (REFFILE, "< $assembled") or die "Could not open $assembled\n";
while(defined (my $l = <REFFILE>)) {
	chomp ($l);
  if ($l =~ /^>/) {
		$id = $l;
		# strip > character
		$id =~ s/^.//;
    if ($id ne "") {
      $seqs{$id} = $seq;
    }
    $seq = "";
	}
  else {
    $seq = $seq . $l;
  }
}
close REFFILE;

my %aa_seqs;
$id = "";
$seq = "";
open (REFFILE, "< $constructed") or die "Could not open $constructed\n";
while(defined (my $l = <REFFILE>)) {
	chomp ($l);
  if ($l =~ /^>/) {
		$id = $l;
		# strip > character
		$id =~ s/^.//;
    if ($id ne "") {
      $seqs{$id} = $seq;
    }
    $seq = "";
	}
  else {
    $seq = $seq . $l;
  }
}
close REFFILE;

open (COUNTFILE, "< $counts") or die "Could not open $counts\n";
open (OUTFILE, "> $fa_out") or die "Could not open $fa_out\n";
while(defined (my $l = <COUNTFILE>)) {
	chomp ($l);
  if ($l =~ /^([^\t]+)\t/) {
    $id = $1;
		print "id: " . $id . "\n";
    if (exists $seqs{$id}) {
      print OUTFILE ">" . $id . "\n";
      print OUTFILE $seqs{$id} . "\n";
			print "found id in seqs\n";
    }
    else {
      %seqs = %aa_seqs;
      print "Switched to amino acid reference at " . $l . "\n";
      close OUTFILE;
      open (OUTFILE, "> $faa_out") or die "Could not open $faa_out\n";
    }
  }
  else {
    print "Unable to parse " . $l . "\n";
  }
}
close OUTFILE;
close COUNTFILE;

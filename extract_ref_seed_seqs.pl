#!/usr/bin perl -w
use strict;

my $help = "\tFirst argument is a path to the BLAST output (outfmt 6).\n
Second argument is the refseqs fasta\n
Third argument is the output file for the sections of ref sequences matching the SEED database\n
Fourth argument is the output file for all remaining unmatched ref sequence segments over 500 nt long\n";
print $help if !$ARGV[0];exit if !$ARGV[0];
print $help if !$ARGV[1];exit if !$ARGV[1];
print $help if !$ARGV[2];exit if !$ARGV[2];
print $help if !$ARGV[3];exit if !$ARGV[3];

my $blastout = $ARGV[0];
my $refseq = $ARGV[1];
my $matchout = $ARGV[2];
my $unmatchout = $ARGV[3];

my %seqs;
my $id = "";
my $seq = "";
open (SEQFILE, "< $refseq") or die "Could not open $refseq\n";
while(defined (my $l = <SEQFILE>)) {
	chomp ($l);
	if ($l =~ /^>/) {
    if ($id ne "") {
      $seqs{$id} = $seq;
    }
		$id = $l;
		# strip > character
		$id =~ s/^.//;
    $seq = "";
	}
  else {
    $seq = $seq . $l;
  }
}
close SEQFILE;

print "Processed sequences in " . $blastout . "\n";

open (MATCH, "> $matchout") or die "Could not open $matchout\n";
open (UNMATCH, "> $unmatchout") or die "Could not open $unmatchout\n";
open (BLAST, "< $blastout") or die "Could not open $blastout\n";
my @lineitems;
my %segments;
my $previd = "";
my $start;
my $end;
my $overlap;
my @starts;
my $append;
my $segmentlength;
my $prevend;
my $key;
while(defined (my $l = <BLAST>)) {
	chomp ($l);
  @lineitems = split(/\t/, $l);
  if (@lineitems == 12) {
    if ($lineitems[6] > $lineitems[7]) {
      $start = $lineitems[7];
      $end = $lineitems[6];
    }
    else {
      $start = $lineitems[6];
      $end = $lineitems[7];
    }
    $id = $lineitems[0];
    if ($id ne $previd) {
      # append seqids of matched segments with numbers, and unmatched segments > 500 long with letters
      $append = 1;
      @starts = keys %segments;
      foreach $key (@starts) {
        print MATCH ">" . $previd . "_" . $append . "\n";
        $segmentlength = $segments{$key} - $key + 1;
        $seq = substr $seqs{$previd}, $key, $segmentlength;
        print MATCH $seq . "\n";
        ++$append;
      }
      @starts = sort @starts;
      $append = 'A';
      $prevend = -1;
      foreach $key (@starts) {
        if ($key - $prevend > 500) {
          print UNMATCH ">" . $previd . "_" . $append . "\n";
          $segmentlength = $key - $prevend - 1;
          $seq = substr $seqs{$previd}, $prevend + 1, $segmentlength;
          print UNMATCH $seq . "\n";
          ++$append;
        }
        $prevend = $segments{$key};
      }
      if (length($previd) - $prevend > 500) {
        print UNMATCH ">" . $previd . "_" . $append . "\n";
        $segmentlength = length($previd) - $prevend - 1;
        $seq = substr $seqs{$previd}, $prevend + 1, $segmentlength;
        print UNMATCH $seq . "\n";
        ++$append;
      }
      %segments = ();
      $segments{$start} = $end;
      $previd = $id;
    }
    else {
      $overlap = 0;
      @starts = keys %segments;
      foreach $key (@starts) {
        if (($start > $key and $start < $segments{$key}) or ($end > $key and $end < $segments{$key})) {
          $overlap = 1;
          last;
        }
      }
      if (!$overlap) {
        $segments{$start} = $end;
      }
    }    
  }
}
close MATCH;
close UNMATCH;
close BLAST;
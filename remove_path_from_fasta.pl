#!/usr/bin/env perl -w
#use strict;
my $help = "\tFirst argument is a path to the assembled sequences\n
Second argument is path to output file\n";

print $help if !$ARGV[0];exit if !$ARGV[0];
print $help if !$ARGV[1];exit if !$ARGV[1];

my $seqs = $ARGV[0];
my $output = $ARGV[1];
my $fullid;
open (SEQFILE, "< $seqs") or die "Could not open $seqs\n";
open (OUT, "> $output") or die "Could not open $output\n";
while(defined (my $l = <SEQFILE>)) {
	chomp ($l);
	if ($l =~ /^>/) {
		$fullid = $l;
    if ($fullid =~ /^(*)_path*/) {
      print OUT $1 . "\n";
    }
    else {
      print OUT $fullid .  "\n";
    }
	}
  else {
    print OUT $l . "\n";
  }
}
close SEQFILE;
close OUT;
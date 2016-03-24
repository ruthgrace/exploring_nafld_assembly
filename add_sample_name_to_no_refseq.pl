#!/usr/bin/env perl -w
#use strict;
my $help = "\tFirst argument is a path to the assembled sequences\n
Second argument is the sample name to be prepended to sequence ids\n";
print $help if !$ARGV[0];exit if !$ARGV[0];
print $help if !$ARGV[1];exit if !$ARGV[1];

my $seqs = $ARGV[0];
my $sample = $ARGV[1];
my $fullid;
open (SEQFILE, "< $seqs") or die "Could not open $seqs\n";
while(defined (my $l = <SEQFILE>)) {
	chomp ($l);
	if ($l =~ /^>/) {
		$fullid = $l;
		# strip > character
		$fullid =~ s/^.//;
    print ">" . $sample . "_" . $fullid .  s"\n";
	}
  else {
    print "$l\n";
  }
}
close SEQFILE;

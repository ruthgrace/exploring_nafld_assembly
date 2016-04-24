#!/usr/bin/env perl -w
#use strict;
my $help = "\tFirst argument is a path to the annotated counts table (first column is refseq ids)\n
Second argument is output file path\n";
print $help if !$ARGV[0];exit if !$ARGV[0];
print $help if !$ARGV[1];exit if !$ARGV[1];

my $counts = $ARGV[0];
my $out = $ARGV[1];

# read in all refseq lengths
my $seqlength = "";
my @countlineitems;
my $firstline = 1;
open (COUNTFILE, "< $counts") or die "Could not open $counts\n";
open (OUTFILE, "> $out") or die "Could not open $out\n";
while(defined (my $l = <COUNTFILE>)) {
	chomp ($l);
	if ($firstline) {
		$firstline = 0;
		my @lineitems = split(/\t/, $l, 2);
		print OUTFILE "$lineitems[0]\tlength\t$lineitems[1]\n";
	}
	else {
		@countlineitems = split(/\t/, $l, 2);
		if (@countlineitems == 2) {
      if ($countlineitems[0] =~ m/_len=([^_])+/g) {
        $seqlength = $1;
        print OUTFILE $countlineitems[0] . "\t" . $seqlength . "\t" . $countlineitems[1] . "\n";
      }
      else {
        print "was not able to parse seqlength for line $l\n";
      }
			
		} else {
			print "Couldn't parse line $l\n";
		}
	}
}
close COUNTFILE;
close OUTFILE;
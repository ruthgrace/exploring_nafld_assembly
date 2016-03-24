#!/usr/bin/env perl -w
#use strict;
my $help = "\tFirst argument is a path to the seed hierarchy file\n
Second argument is the sample name to be prepended to the lines in the seed hierachy file\n";
print $help if !$ARGV[0];exit if !$ARGV[0];
print $help if !$ARGV[1];exit if !$ARGV[1];

my $seed_hier = $ARGV[0];
my $sample = $ARGV[1];

my $firstline = 1;

open (SEEDHIER, "< $seed_hier") or die "Could not open $seed_hier\n";
while(defined (my $l = <SEEDHIER>)) {
	chomp ($l);
	if ($firstline) {
		$firstline = 0;
    print $l . "\n";
	}
	else {
    print $sample . "_" . $l . "\n";
	}
}
close SEEDHIER;

#!/usr/bin/perl -w
#use strict;
my $help = "\tFirst argument is a path to the sam file\n
Second argument is path to output file\n";

print $help if !$ARGV[0];exit if !$ARGV[0];
print $help if !$ARGV[1];exit if !$ARGV[1];

my $sam = $ARGV[0];
my $output = $ARGV[1];
my $fullid;
open (SAMFILE, "< $sam") or die "Could not open $sam\n";
open (OUT, "> $output") or die "Could not open $output\n";
while(defined (my $l = <SAMFILE>)) {
	chomp ($l);
  if ($l =~ /(TRINITY.+)_path[^\t]+(\t.*)$/) {
    print OUT $1 . $2 . "\n";
  }
  else {
    print OUT $l .  "\n";
  }
}
close SAMFILE;
close OUT;
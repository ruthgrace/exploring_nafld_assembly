#!/usr/bin/perl
#use strict;

# June 15, 2013 - JM
#--------------------------------------------------------------------------------------------------
# This will go through all _CDS_counts.txt tables and make list of unique refseqs
#	for BLASTing the SEED
#--------------------------------------------------------------------------------------------------

my %refseqs;

my @map = glob("*_map");
	foreach(@map){
	chomp($_);
		my @table = glob("$_/*_counts.txt");
			open (INPUT, "< $table[0]") or die "Could not open $table[0]: $!\n";
			while(my $l = <INPUT>) {
			chomp ($l);
				my @hold = split (/\t/, $l);
				if ($l !~ m/^refseq/ && !exists $refseqs{$hold[0]}){
					$refseqs{$hold[0]} = "";
				}
			}close INPUT;
	}

foreach(keys(%refseqs)){		
	print "$_\n";
}

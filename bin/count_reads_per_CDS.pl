#!/usr/bin/perl
#use strict;

#####
# Mar 2, 2011 - JM
#
# Use the best_hits file after mapping to add up number of reads per CDS (unique reference ID)
#
#####


my $best_hits_file = $ARGV[0];

my %count_fwd; my %count_rev; my %seen;


open(IN, "< $best_hits_file") or die;
while(defined(my $l = <IN>)){
	chomp $l;
		my @l = split(/\t/, $l);
		$seen{$l[0]} = "";
		$count_fwd{$l[0]} = 0 if !exists $count_fwd{$l[0]};
		$count_fwd{$l[0]} = $count_fwd{$l[0]} + $l[2];

		$count_rev{$l[0]} = 0 if !exists $count_rev{$l[0]};
		$count_rev{$l[0]} = $count_rev{$l[0]} + $l[3];
		
	}
close IN;

print "refseqID\tfwd_reads\trev_reads\n";
foreach (keys(%seen)){
	$count_fwd{$_} = 0 if !exists $count_fwd{$_};
	$count_rev{$_} = 0 if !exists $count_rev{$_};

	print "$_\t$count_fwd{$_}\t$count_rev{$_}\n";
}

#!/usr/bin/perl
use strict;

#####
# Mar 2, 2011 - JM
#
# 1. parse out names of FASTA files: grep '^@SQ' 4_map/map_best1_trim10_4.out > headers.txt
# 2. Use this to add gene length information to the read counts tables (assuming column order: ref_id;fwd_reads;rev_reads;total_reads)
#####

my $headers = $ARGV[0];
my $cds_counts = $ARGV[1];

my %length;

open(IN, "< $headers") or die;
while(defined(my $l = <IN>)){
	chomp $l;
		my @l = split(/\t/, $l);
		$l[1] =~ s/SN://;
		$l[-1] =~ s/LN://;
		
		$length{$l[1]} = $l[-1];

		
	}
close IN;

open(IN, "< $cds_counts") or die;
while(defined(my $l = <IN>)){
	chomp $l;
	print "refseqID\tlength\tfwd_reads\trev_reads\ttotal_reads\n" if $l =~ /^refseq/;
	if ($l !~ /^refseq/){
		my @l = split(/\t/, $l);

#L. crispatus fix: Apr 19, 2011		
#		my $blah = join("\t", @l[2..13]);
#		print "$l[0]\t$length{$l[0]}\t$blah\n";
		print "$l[0]\t$length{$l[0]}\t$l[1]\t$l[2]\t$l[3]\n";
		
	}
	}
close IN;

#!/usr/bin/env perl -w
use strict;

my $help = "\tFirst argument is a path to a table of refseq ids in the first column (up to the first underscore in the ID) and columns for subsys4, subsys1, subsys2, and subsys3.\n
Second argument is a path to directory with one folder with one count file inside per sample, where the count file has columns for refseqID, length, fwd_reads, rev_reads, and total_reads\n
Third argument is output file path\n";
print $help if !$ARGV[0];exit if !$ARGV[0];
print $help if !$ARGV[1];exit if !$ARGV[1];
print $help if !$ARGV[2];exit if !$ARGV[2];

my $annotationFile = $ARGV[0];
my $countsDir = $ARGV[1];
my $outfile = $ARGV[2];


# read in all annotations into hash, keyed by refseq ID
my %annotatedSeqs;

# make separate refseq full ID value hash (keyed just by the unique identifier number)
my %sampleSeqs;

# put refseqs from sample into hash
open (ANNOTATIONS, "< $annotationFile") or die "Could not open $annotationFile\n";
my $firstline = 1;
my $id;
my $annotatedSeqs;
my $annotation;
my @tempArray;
while(defined (my $l = <ANNOTATIONS>)) {
	chomp ($l);
	if ($firstline) {
		$firstline = 0;
	}
	else {
		@tempArray = split(/\t/, $l, 2);
		if (scalar(@tempArray)==2) {
			$id = $tempArray[0];
			$annotation = $tempArray[1];
			$annotatedSeqs{$id} = $annotation;
		}
		else {
			print "Unable to parse line $l\n";
		}
	}
}
close ANNOTATIONS;

my @header = ('subsys4','subsys1','subsys2','subsys3');

# process samples
$countsDir =~ s/\/$//;
my @mappedSamples = glob($countsDir . '/*map');
my $countFile;
my $sample;
my %fullIDs;
my %counts;
foreach my $mappedSample (@mappedSamples) {
	$sample = $mappedSample;
	print "sample dir: " . $sample . "\n";
	if (index($sample, ".txt") == -1) {
		$sample =~ s/\/$//; # remove trailing slash
		$sample = (split(/\//,$sample))[-1];
		$sample =~ s/_map$//;
		print "sample: " . $sample . "\n";
		$countFile = "${mappedSample}/${sample}_CDS_counts.txt";
		print "count file: " . $countFile . "\n";
		%fullIDs = ();
		%counts = ();
		open (COUNTS, "< $countFile") or die "Could not open $countFile\n";
		$firstline = 1;
		while(defined (my $l = <COUNTS>)) {
			chomp ($l);
			if ($firstline) {
				$firstline = 0;
			}
			elsif (length($l)>0) {
				$id = ( split('\t', $l ))[0];
				$counts{$id} = ( split('\t', $l ))[ -1 ];
				print "adding id " . $id . " with count " . $counts{$id} . "\n";
			}
		}
		close COUNTS;
		# prepend sample name onto header
		unshift @header, $sample;
		# for each annotation refseq, prepend 0 or count onto value, tab separated, add refseq full ID value hash if doesn't exist
		for my $key (keys %annotatedSeqs) {
			if (exists($counts{$key})) {
				$annotatedSeqs{$key} = "$counts{$key}\t$annotatedSeqs{$key}";
			}
			else {
				$annotatedSeqs{$key} = "0\t$annotatedSeqs{$key}";
			}
		}
	}
}

# print out full refseq with counts & annotations
open(OUT, '>', $outfile) or die "Could not open file '$outfile' $!";
unshift @header, "refseq";
print OUT join("\t", @header), "\n";
for my $key (keys %annotatedSeqs) {
	print OUT $key . "\t" . $annotatedSeqs{$key} . "\n";
}

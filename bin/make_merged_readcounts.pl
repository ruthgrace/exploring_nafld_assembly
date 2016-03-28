#!/usr/bin/perl
#use strict;

# June 18, 2013 - JM
#--------------------------------------------------------------------------------------------------
# MUST run this where you did the mapping (e.g. map_bac/). Assumes all mapping dirs are named *_map
#
# This will go through all _CDS_counts.txt tables and make list of unique refseqs. These will be
#	printed out in rows with the readcounts/sample in columns. If length column was present the
#	refseq length will be added as a column. Otherwise length will be 1.
#
# Run: /Volumes/rhamnosus/twntyfr/bin/make_merged_readcounts.pl > merged_readcounts.txt
#
# ALTERNATE: There are 3 lines of code to comment out if you want to keep the strand info. Otherwise
#	total read counts are put in the table. Look for "#total"
#JEAN: One day add flags to input so user can choose
#--------------------------------------------------------------------------------------------------

my %refseqs;				# Keep track of the seen refseqs
my @sample_list;			# Keep track of samples in order
my @header;				#for fwd/rev
my %counts;				# Read counts per sample/per refseq
my $length = "N";			# Check for length column


my @map = glob("*_map");						# Look for all directories named *_map
	foreach(@map){
	chomp($_);
	my @hold = split(/\_map/, $_);
	my $name = $hold[0];								# sample ID
	push (@sample_list, $name);							# keep list of sampleIDs in order
	push (@header, $name . "_fwd");						# keep list of sampleIDs in order for header. Add strand info.
	push (@header, $name . "_rev");						# ""
		my @table = glob("$_/*_counts.txt");			# Look for readcounts table
#refseqID	length	fwd_reads	rev_reads	total_reads

			open (INPUT, "< $table[0]") or die "Could not open $table[0]: $!\n";
			while(my $l = <INPUT>) {
			chomp ($l);
				my @hold = split (/\t/, $l);
				if ($l !~ m/^refseq/ && !exists $refseqs{$hold[0]}){			# hold[0] is the refseqID. Keep unique refseqIDs
					$refseqs{$hold[0]} = "1" if $length ne "Y";
					$refseqs{$hold[0]} = $hold[1] if $length eq "Y";
				}
				if ($l !~ m/^refseq/){
					my $lookup = join("\t", $name, $hold[0]);					# Make the lookup a combination of sampleId and refseqID: $001A	refseqidX
					$counts{$lookup} = $hold[-1];								# TOTAL counts for the sample+refseq
#1					$counts{$lookup} = "$hold[-3]\t$hold[-2]";						# separate fwd/rev counts					
				}elsif ($l =~ m/^refseq/ && $l =~ m/\tlength/){					# Check for length column
					$length = "Y";
				}
			}close INPUT;
	}

# List of ordered samples
#2	my $s_header = join("\t", @header);
	my $s_header = join("\t", @sample_list);
print "refseqID\tlength\t$s_header\n";

foreach my $ref (keys(%refseqs)){								# For each unique refseq
my @hold;
	foreach my $sample (@sample_list){							# For each sample ID (per refseqs)
		my $temp = join("\t", $sample, $ref);					# Put all read counts for that refseq in sample order
		push (@hold, $counts{$temp}) if (exists $counts{$temp});
		push (@hold, "0") if (!exists $counts{$temp});			# Add zero counts if that refseq wasn't mapped for a sample
#3		push (@hold, "0\t0") if (!exists $counts{$temp});			# Add zero counts (for fwd AND rev) if that refseq wasn't mapped for a sample
	}
	my $counts = join("\t", @hold);

	print "$ref\t$refseqs{$ref}\t$counts\n";					# Print tab-delimited table
}




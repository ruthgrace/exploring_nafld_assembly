#!/usr/bin/perl
##use strict;

# June 19, 2013 - JM

# Adds the .faa ID to the readcounts table (containing the refseqid)

my $table = $ARGV[0];
my $lookup = $ARGV[1];

#my $table = "/Volumes/rhamnosus/twntyfr/map_bac/merged_readcounts.txt";
#my $lookup = "/Volumes/rhamnosus/reference_genomes/may2013/wget/all_ffn_faa_lookup.txt";
#ffnID	faaID
#>


my %lookup;
open (INPUT, "< $lookup") or die "Could not open $lookup: $!\n";
while(defined (my $l = <INPUT>)) {
chomp ($l);
	if ($l !~ m/^ffnID/){
		my @hold = split(/\t/, $l);
		my $ffn = &fixheaders($hold[0]); my $faa = &fixheaders($hold[1]);
#print "$ffn\t$faa\n";close INPUT;exit;
		$lookup{$ffn} = $faa;				# faa lookup by ffn
	}
}close INPUT;

open (INPUT, "< $table") or die "Could not open $table: $!\n";
while(defined (my $l = <INPUT>)) {
chomp ($l);
	if ($l !~ m/^refseq/){
		my @hold = split(/\t/, $l);
		$hold[0] =~ s/>//;
		print "$lookup{$hold[0]}\t$l\n" if exists $lookup{$hold[0]};
		print "noID\t$l\n" if !exists $lookup{$hold[0]};
	}else{
		print "refseqIDfaa\t$l\n";
	}
}close INPUT;

sub fixheaders{
	my @hold = split(/\s+/, $_[0]);
	$hold[0] =~ s/>//;
	return $hold[0];
}

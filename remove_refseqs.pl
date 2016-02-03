#!/usr/bin/perl -w
use strict;
use Cwd;
use warnings;

my $infile = $ARGV[0] if $ARGV[0];
my $outfile = $ARGV[1] if $ARGV[1];
my $blastfile = $ARGV[2] if $ARGV[2];
my $badfile = $ARGV[3] if $ARGV[3];
print "This function takes in four arguments: the path of the input fasta file, and the path of the output fasta file, the path of the blast file with the sequences you want to remove, and the path of the file you want to put the removed sequences in\n" if !$ARGV[0];
print "This function takes in four arguments: the path of the input fasta file, and the path of the output fasta file, the path of the blast file with the sequences you want to remove, and the path of the file you want to put the removed sequences in\n" if !$ARGV[1];
print "This function takes in four arguments: the path of the input fasta file, and the path of the output fasta file, the path of the blast file with the sequences you want to remove, and the path of the file you want to put the removed sequences in\n" if !$ARGV[2];
print "This function takes in four arguments: the path of the input fasta file, and the path of the output fasta file, the path of the blast file with the sequences you want to remove, and the path of the file you want to put the removed sequences in\n" if !$ARGV[3];
exit if !$ARGV[0];
exit if !$ARGV[1];
exit if !$ARGV[2];
exit if !$ARGV[3];

my %badseqs;
my @items;
my $start;
my $end;
my $qlen;
my $trinityseqid;
# open the blast file first and extract the sequence IDs of all things with hit results
open (BLAST, "< $blastfile") or die "$!\n";
while(defined(my $l = <BLAST>)){
  chomp $l;
  @items=split(/\t/,$l);
  # blast tab format is
  # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
  # 0      1      2      3      4        5       6      7    8      9    10     11       12   13
  print "trinitiy seq id is $items[0]";
  $trinityseqid = $items[0];
  print "qstart is $items[6]";
  $start = $items[6];
  print "qend is $items[7]";
  $end = $items[7];
  print "qlen is $items[12]";
  $qlen = $items[12];
  if (abs($start - $end) / $qlen >= 0.9) {
    $badseqs{$trinityseqid} = "trinity seq id";
  }
}	
close BLAST;

my %seqs;
my $lastID = "";

# open the input file and output file and write only sequences that DON'T match a refseq
open (IN, "< $infile") or die "$!\n";
unless(open OUT, '>', $outfile) {
    print "\nUnable to open $outfile: $!\n";
}
unless(open BAD, '>', $badfile) {
    print "\nUnable to open $badfile: $!\n";
}
while(defined(my $l = <IN>)){
	chomp $l;
	if ($l =~ /^>/) {
		if ($lastID ne "") {
      if (exists $badseqs{$lastID}) {
				print BAD "$lastID\n";
				print BAD "$seqs{$lastID}";
      }
      else {
				print OUT "$lastID\n";
				print OUT "$seqs{$lastID}";
      }
		}
		$lastID = $l;
    $seqs{$lastID} = "";
	}
	else {
    $seqs{$lastID} .= "$l\n";
  }
}	
close IN;
close OUT;
close BAD;

#!/usr/bin/perl
#use strict;

#map ABI reads from Bowtie to chromosome positions
#output a table of reads per position for each unique accession number
#accession pos fr rr

#this output is fed into the feature table parser 

my %fcount; #holds a hash of arrays keyed by accession and holding counts
my %rcount; #holds a hash of arrays keyed by accession and holding counts

#this reads the SAM output and captures every forward and reverse mapping read

open (IN, "< $ARGV[0]") or die;
	while(defined(my $l = <IN>)){
		chomp $l;
		my @l = split/\t/, $l;
		my $acc = $l[2];
		my $orient = $l[1];
		my $pos = $l[3];
		my $seq= $l[9];
		
		#SAM output uses @ for comments. WHY!!!!!!!
		if ($l !~ /^\@/){
		
			#key is accession, postion in the array is the position in the chromosome
			#value is the sum of the counts of reads starting a given position in the accession
			if ( $orient == 0){
				${ $fcount{$acc} }[$pos]++;
			}	
			
			#BOWTIE RC's the sequence and reports mapping to the forward strand
			#from the manual
			#read starts at the end of the sequence for reverse reads
			if ( $orient == 16){
				my $n = length($seq);
				${ $rcount{$acc} }[$pos + $n]++;
			}		
			
		}
		
	}
close IN;

#get the accession numbers
my @fids = (keys(%fcount));
my @rids = keys(%rcount);
my %union;
foreach(@fids){ $union{$_} = ""; }
foreach(@rids){ $union{$_} = ""; }
my @union = sort(keys(%union));


foreach(@union){
	my $size = 0;
	$size = @{ $fcount{$_} } if exists $fcount{$_} ;
	$size = @{ $rcount{$_} } if (exists $rcount{$_} && @{ $rcount{$_} } > $size);
	
	for (my $i = 0; $i <$size; $i++){
		my $fr = my $rr = 0;
		if (exists ${ $fcount{$_} }[$i] or exists ${ $rcount{$_} }[$i]){
			$fr =  ${ $fcount{$_} }[$i] if exists  ${ $fcount{$_} }[$i];
		}
		if (exists ${ $rcount{$_} }[$i]){
			$rr = ${ $rcount{$_} }[$i];
		}
		print "$_\t$i\t$fr\t$rr\n";
	
		my %nt; $nt{A}=0;$nt{C}= 0; $nt{G}=0;$nt{T}=0;
	}
}

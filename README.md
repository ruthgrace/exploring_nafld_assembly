# exploring_nafld_assembly

## BLAST to annotations

Blasting all the assembled sequences to the SEED database for annotations would take a super long time because

1) there is a lot of data
2) the assembled sequences are nucleotide, and you'd have to blast all the reading frames

so the strategy we are going with is that we'll take all the annotated sequences, and exclude anything that matches our reference library (which has already been blasted with the SEED database) at 90% identity with 90% of the length of the sequence.

#### make blast db

I was super dumb and had spaces in my refseq IDs. I even had a space in between the '>' and the first thing in each ID. Since BLAST takes everything up to the first whitespace as the identifier, all my sequence IDs were blank the first time I ran this :(

To remove the space between the ">" and the unique number for each identifier, run:

```
sed -i.backup '/^>/ s/^> />/' all_genus_orfs_clustered_at_99_unique.fa 
```

To replace the rest of the spaces in the sequence identifiers with underscores, run:

```
sed -i.backup '/^>/ s/ /_/g' all_genus_orfs_clustered_at_99_unique.fa 
```

Trinity, the assembly program I used, also likes to make sequence identifiers with spaces. To remedy, run:

```
sed -i.backup '/^>/ s/ /_/g' ../scripts/trinity_CL_119_R1/Trinity.fasta 
```

To make the indexed blast database, run

```
makeblastdb -in all_genus_orfs_clustered_at_99_unique.fa -dbtype nucl
```

where `all_genus_orfs_clustered_at_99_unique.fa` is your reference sequences.

### Single sample

#### blast with refseq

I wrote a one off script just to test one sample. I'll just do one sample for now, and if the assembly is dredging up a lot of annotations that I've missed previously, I'll do all the other samples too in some batch script. The sample data file and blast database are hard coded in.

```
nohup ./assembly_test_blast.sh > CL_119_R1_blast_to_refseqs_nohup.out 2>&1&
```

#### remove the matches to refseq

This perl script takes in four parameters: the path of the input fasta file, and the path of the output fasta file, the path of the blast file with the sequences you want to remove, and the path of the file you want to put the removed sequences in.

```
nohup perl ./remove_refseqs.pl ../scripts/trinity_CL_119_R1/Trinity.fasta ./CL_119_R1_assembled_no_refseqs.fa CL_119_R1_blast_to_refseqs.out CL_119_R1_assembled_matched_refseqs.fa > remove_refseqs_CL_119_R1_nohup.out 2>&1&
```

#### translate blast remaining sequences to SEED

This is done with blastx, which translates all reading frames of the remaining assembled sequences to the SEED database for annotation.

```
nohup ./blastx_to_SEED.sh > blastx_to_SEED_CL_119_R1_nohup.out 2>&1&
```

### All samples

I wrote a batch script to do the steps outlined above for all the samples. In the script I am using 40 threads, so I've set the niceness to be 1 more than normal in case other people need to run things:

```
nohup nice -n 1 ./assembly_blast.sh > assembly_blast_nohup.out 2>&1&
```

## Create fasta file for bowtie index

Add sample names to all assembled refseq ids:

```
nohup perl add_all_samples_name_to_no_refseq.pl > add_all_samples_name_to_no_refseq_nohup.out 2>&1&
```

Concatenate all assembled refseqs:

```
cat */*_assembled_no_refseqs_with_sample_name.fa > all_assembled_no_refseqs_with_sample_name.fa
```

## Get seed hierarchy for all samples

```
nohup get_all_seed_hier.sh /Volumes/data/ruth/nafld_assembly/assembly_test_blast > get_all_seed_hier_nohup.out 2>&1&
```

## Map reads to assembled stuff

Run mapping:

```
nohup mapping.sh > mapping_nohup.out 2>&1&
```

Get counts:

```
nohup count_table_gen.sh > count_table_gen_nohup.out 2>&1&
```

## Generate count table for reads mapping to assembled annotated seqs

```
nohup perl ../../scripts/get_annotated_counts.pl $seedhier /Volumes/data/ruth/nafld_assembly/assembly_mapping get_annotated_counts_output.txt > get_annotated_counts_nohup.out 2>&1&
```

Add refseq lengths:

```
nohup perl add_assembled_length.pl get_annotated_counts_output.txt annotated_assembly_counts_with_seq_length.txt > add_assembled_length_nohup.out 2>&1&
```

## Amalgamate count tables to count tables for reads mapping to the library

This script concatenates the count tables, ensuring that the samples are in the same order.

```
Rscript amalgamate_counts.r annotated_assembly_counts_with_seq_length.txt annotated_counts_with_refseq_length.txt all_counts_with_seq_length.txt
```

## Remove zeros

## Aitchison transform

## Differential abundance analysis
# exploring_nafld_assembly

## BLAST to annotations

Blasting all the assembled sequences to the SEED database for annotations would take a super long time because

1) there is a lot of data
2) the assembled sequences are nucleotide, and you'd have to blast all the reading frames

so the strategy we are going with is that we'll take all the annotated sequences, and exclude anything that matches our reference library (which has already been blasted with the SEED database) at 90% identity with 90% of the length of the sequence.

### make blast db

I was super dumb and had spaces in my refseq IDs. I even had a space in between the '>' and the first thing in each ID. Since BLAST takes everything up to the first whitespace as the identifier, all my sequence IDs were blank the first time I ran this :(

To remove the space between the ">" and the unique number for each identifier, run:

```
sed -i.backup '/^>/ s/^> />/' all_genus_orfs_clustered_at_99_unique.fa 
```

To replace the rest of the spaces in the sequence identifiers with underscores, run:

```
sed -i.backup '/^>/ s/ /_/g' all_genus_orfs_clustered_at_99_unique.fa 
```

To make the indexed blast database, run

```
makeblastdb -in all_genus_orfs_clustered_at_99_unique.fa -dbtype nucl
```

where `all_genus_orfs_clustered_at_99_unique.fa` is your reference sequences.

### blast with refseq

I wrote a one off script just to test one sample. I'll just do one sample for now, and if the assembly is dredging up a lot of annotations that I've missed previously, I'll do all the other samples too in some batch script. The sample data file and blast database are hard coded in.

```
nohup ./assembly_test_blast.sh > CL_119_R1_blast_to_refseqs_nohup.out 2>&1&
```

### remove the matches to refseq

### translate blast remaining sequences to SEED


# exploring_nafld_assembly

## Shortening sequence id
I realized that I had some super long sequence IDs later in this process, so I removed everything `_path` and after in the sequence IDs in the reference library `all_assembled_no_refseqs_with_sample_name.fa` (using `remove_path_from_fasta.pl`) and the annotation `$sample"/"$sample"_SEED_blast.out"` (using `remove_path_from_all_blast_out.sh`). I also had to clean up the mapping data with `./remove_path_from_all_sam.sh` because I had already done the mapping when I decided to change the seq ids. I recommend that this be done on each sample for any Trinity output sequences, because sometimes the path value included in the sequence ID can be very long, and this seriously affects the size and processing speed for files later in the analysis process.

```
nohup perl remove_path_from_fasta.pl "../scripts/trinity_"$sample$"/Trinity.fasta" "../scripts/trinity_"$sample$"/Trinity_no_path.fasta" > remove_path_from_fasta_nohup.out 2>&1&
```

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

Extract non overlapping portions of reference sequence that match the SEED database. The argument is a directory containing a directory for each sample, which has the BLAST output.

```
nohup ./extract_all_ref_seed_seqs.sh /Volumes/data/ruth/nafld_assembly/assembly_test_blast > extract_all_ref_seed_seqs_nohup.out 2>&1&
```

Recursively BLAST non matching sequences > 500 long to SEED. Files output from each round of recursion will be output with the recursion count number.

```
nohup ./recursive_blast_to_seed.sh /Volumes/data/ruth/nafld_assembly/assembly_test_blast > recursive_blast_to_seed_nohup.out 2>&1&
```

Concatenate refseq matches (fasta file and seed blast out) per sample

```
nohup ./cat_all_refseq_matches.sh /Volumes/data/ruth/nafld_assembly/assembly_test_blast > cat_all_refseq_matches_nohup.out 2>&1&
```

Add sample names to refseq matches per sample

```
nohup ./add_all_sample_names_to_refseq_match_fasta.sh > add_all_sample_names_to_refseq_match_fasta_nohup.out 2>&1&
```

Concatenate all refseq SEED match sequences

```
cat */*all_refseq_matches_with_sample_name.fasta > all_recursive_refseq_matches.fasta
```

Build bowtie index inside `seed_matches` folder

```
mkdir seed_matches
cd seed_matches
bowtie2-build ../all_recursive_refseq_matches.fasta seed_matches
```

## Get seed hierarchy for all samples

```
nohup ./get_all_seed_hier.sh /Volumes/data/ruth/nafld_assembly/assembly_test_blast > get_all_seed_hier_nohup.out 2>&1&
```

## Map reads to assembled stuff

Run mapping:

```
nohup ./mapping.sh > mapping_nohup.out 2>&1&
```

Get counts:

```
nohup ./count_table_gen.sh > count_table_gen_nohup.out 2>&1&
```

If you mess up getting the counts and need to re run it, run this in your mapping folder first to clean up:

```
rm */best_hits*
rm */*_CDS_counts.txt
rm headers.txt
```

## Generate count table for reads mapping to assembled annotated seqs

Put all `$sample/$sample_get_seed_hier_output_with_sample_name.txt` into same file:

```
nohup cat */*_get_seed_hier_output_with_sample_name.txt > all_seed_hier_output_with_sample_name.txt 2>&1&
```

Manually remove the following lines throughout `all_seed_hier_output_with_sample_name.txt`:

```
nohup: ignoring input
refseq_id       subsys4 subsys1 subsys2 subsys3
```

```
nohup perl get_annotated_counts.pl all_seed_hier_output_with_sample_name.txt /Volumes/data/ruth/nafld_assembly/assembly_mapping get_annotated_counts_output.txt > get_annotated_counts_nohup.out 2>&1&
```

Add refseq lengths:

```
nohup perl add_assembled_length.pl get_annotated_counts_output.txt annotated_assembly_counts_with_seq_length.txt > add_assembled_length_nohup.out 2>&1&
```

## Amalgamate count tables to count tables for reads mapping to the library

This script concatenates the count tables, ensuring that the samples are in the same order.

```
nohup Rscript amalgamate_counts.r annotated_assembly_counts_with_seq_length.txt /Volumes/data/ruth/scripts/annotated_counts_with_refseq_length.txt all_counts_with_seq_length.txt > amalgamate_counts_nohup.out 2>&1&
```

## Remove zeros

Remove all features with zero counts for all samples

```
nohup Rscript remove_zeros.r all_counts_with_seq_length.txt all_counts_with_seq_length_zeros_removed.txt > remove_zeros_nohup.out 2>&1&
```

## Aitchison transform + differential abundance analysis

Run `run_aldex_for_stripcharts.r` line by line in R (for some reason it won't run all in one go but works fine line by line. You'll have to run the Aitchison stuff line by line too.)

The Aitchison transform script amalgamates counts to unique subsys hierarchies. The ALDEx (Anova-Like Differential Expression) analysis calculates the variance within groups for each feature and the variance between groups for each feature to arrive upon an effect size.

After running this script, you will get some files that are output in intermediate stages (so that you can read them in and start the script from the middle if necessary), plus `ALDEx_all_hierarchies_output.pdf`, which has a difference within vs. difference between plot for subsys4, and `ALDEx_output_for_stripcharts_ordered_by_effect.txt`, which you can examine to see which features had the greatest measured effect, and which you can also feed into the stripchart plots in the next section.

## Stripchart plots

Features can fall into multiple subsys1, subsys2, and subsys3 categories for each subsys4 category, and can be better visualized with stripcharts. Stripcharts can be generated using the `SEED_stripcharts_aldex2_update.R` script.

## Carb Lipid analysis

Code has been inserted in the `run_aldex_for_stripcharts.r` and `SEED_stripcharts_aldex2_update.R` scripts to run analysis on a carb only or lipid only subset (based on the subsys4 category). This script makes PCA biplots.

`Rscript carb_lipid_pca_plots.r`

## Annotate reference library sequences that matched at least one read with KEGG

This is really annoying because the assembled reference library is in nucleotides and the constructed reference library is in amino acid residues. This perl script extracts all the reference sequences that matched at least one read (in `all_counts_with_seq_length_zeros_removed.txt`), from the assembled library `all_assembled_no_refseqs_with_sample_name.fa` and the constructed library `/Volumes/data/ruth/refseqs/all_genus_orfs_clustered_at_99_unique.fa`, and then puts the DNA reference sequences into `assembled_library_for_kegg.fa` and the amino acid reference sequences into `constructed_library_for_kegg.faa`. This is a super hacky and non reuseable script, but you can use just half of it if you only have one reference library.

```
nohup perl get_seqs_for_kegg.pl all_counts_with_seq_length_zeros_removed.txt all_assembled_no_refseqs_with_sample_name_no_path.fa /Volumes/data/ruth/refseqs/all_genus_orfs_clustered_at_99_unique.fa assembled_library_for_kegg.fa constructed_library_for_kegg.faa > get_seqs_for_kegg_nohup.out 2>&1&
```

Use the online web tool to run: http://www.genome.jp/tools/kaas/

I had [fill this in] DNA sequences and [fill this in] amino acid sequences, and they took [fill this in] and [fill this in] to run on KAAS respectively.

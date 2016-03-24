my $d;
for $d in ../scripts/trinity_* ; do
  trinityFasta="$d/Trinity.fasta";
  sample=${d:19};
  noRefseqs=$sample"/"$sample"_assembled_no_refseqs.fa";
  perl add_sample_name_to_no_refseq.pl $sample"/"$sample"_assembled_no_refseqs.fa" $sample > $sample"/"$sample"_assembled_no_refseqs_with_sample_name.fa"
done

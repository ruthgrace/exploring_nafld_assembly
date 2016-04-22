dir=$1; # first parameter is the directory with one folder per sample
for D in `find $dir -type d`
do
  sample=$(basename $D);
  echo "processing $sample";
  seed_blast_out=$sample"/"$sample"_all_refseq_matches_SEED_BLAST.out";
  perl ../../scripts/get_subsys4.pl -i $seed_blast_out > $sample"/"$sample"_get_subsys4_output.txt";
  perl ../../scripts/get_seed_hier.pl -i $sample"/"$sample"_get_subsys4_output.txt" > $sample"/"$sample"_get_seed_hier_output.txt";
  perl add_sample_name_to_seed_hier_output.pl $sample"/"$sample"_get_seed_hier_output.txt" $sample > $sample"/"$sample"_get_seed_hier_output_with_sample_name.txt";
done

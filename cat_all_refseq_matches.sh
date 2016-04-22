dir=$1; # first parameter is the directory with one folder per sample
for D in `find $dir -type d`
do
  sample=$(basename $D);
  echo "processing $sample";
  cat $sample"/"$sample"_refseq_matches"*".fasta" > $sample"/"$sample"_all_refseq_matches.fasta";
  cat $sample"/"$sample"_refseq_matches"*"SEED"*".out" > $sample"/"$sample"_all_refseq_matches_SEED_BLAST.out";
done


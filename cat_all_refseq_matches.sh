dir=$1; # first parameter is the directory with one folder per sample
for D in `find $dir -type d`
do
  sample=$(basename $D);
  echo "processing $sample";
  cat $sample"/"$sample"_refseq_matches_*" > $sample"/"$sample"_all_refseq_matches.fasta";
done

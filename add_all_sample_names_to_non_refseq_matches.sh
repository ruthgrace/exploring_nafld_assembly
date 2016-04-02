for d in "$1"/*/ ; do
  sample=$(basename $d);
  echo "processing $sample";
  nohup add_sample_name_to_no_refseq.pl $d$sample"_refseq_non_matches.fasta" $d$sample"_refseq_non_matches_with_sample_name.fasta" &
done
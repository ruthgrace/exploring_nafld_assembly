for d in "$1"/*/ ; do
  sample=$(basename $d);
  echo "processing $sample";
  nohup perl add_sample_name_to_no_refseq.pl $d$sample"_all_refseq_matches.fasta" $sample > $d$sample"CL_119_R1_all_refseq_matches_with_sample_name.fasta" 2>&1&
done

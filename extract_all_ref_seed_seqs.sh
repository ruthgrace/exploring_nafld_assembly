for d in "$1"/*/ ; do
  sample=$(basename $d);
  echo "processing $sample";
  ./remove_path_from_fasta.pl "../scripts/trinity_"$sample"/Trinity.fasta" "../scripts/trinity_"$sample"/Trinity_no_path.fasta";
  nohup perl extract_ref_seed_seqs.pl $d$sample"_SEED_blast_no_path.out" "../scripts/trinity_"$sample"/Trinity_no_path.fasta" $d$sample"_refseq_matches.fasta" $d$sample"_refseq_non_matches.fasta" $d$sample"_refseq_matches_SEED_blast.out" > $d$sample"extract_ref_seed_seqs_nohup.out" 2>&1&
done

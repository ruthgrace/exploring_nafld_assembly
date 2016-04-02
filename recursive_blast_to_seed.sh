counter=1;
nonmatch="all_refseq_non_matches";
match="all_refseq_matches";
nonmatchfile=$nonmatch"_"$counter".fasta";
seqsize=$(wc -c <"$nonmatchfile")
while [ $seqsize -ge 0 ]; do
  echo "recursion round "$counter", seqfile size "$seqsize;
  counter=$counter+1;
  blastout=$nonmatch"_"$counter"_SEED_blast.out";
  #blastx
  blastx -db ../../../SEED_database/db_fastas.complex_plusold.faa -query $nonmatchfile -out $blastout -outfmt 6 -evalue 1e-3 -num_alignments 1000 -num_threads 44
  #extract matches and non matches
  counter=$counter+1;
  newnonmatchfile=$nonmatch"_"$counter".fasta";
  matchfile=$match"_"$counter".fasta";
  perl extract_ref_seed_seqs.pl $blastout $nonmatchfile $matchfile $newnonmatchfile;
  #update seqsize for non matches file
  nonmatchfile=$newnonmatchfile;
  seqsize=$(wc -c <"$nonmatchfile")
done

for d in "$1"/[CH]*/ ; do
  sample=$(basename $d);
  echo "processing $sample";
  counter=1;
  nonmatch="refseq_non_matches";
  match="refseq_matches";
  nonmatchfile=$d$sample"_"$nonmatch"_"$counter".fasta";
  cp $nonmatch".fasta" $nonmatchfile;
  seqsize=$(wc -c <"$nonmatchfile")
  while [ $seqsize -ge 0 ]; do
    echo "recursion round "$counter", seqfile size "$seqsize;
    counter=$((counter+1));
    blastout=$d$sample"_"$nonmatch"_"$counter"_SEED_blast.out";
    #blastx
    blastx -db ../../../SEED_database/db_fastas.complex_plusold.faa -query $nonmatchfile -out $blastout -outfmt 6 -evalue 1e-3 -num_alignments 100 -num_threads 44
    #extract matches and non matches
    newnonmatchfile=$d$sample"_"$nonmatch"_"$counter".fasta";
    matchfile=$d$sample"_"$match"_"$counter".fasta";
    matchblastout=$d$sample"_"$match"_"$counter"_SEED_BLAST.out";
    perl extract_ref_seed_seqs.pl $blastout $nonmatchfile $matchfile $newnonmatchfile $matchblastout;
    #update seqsize for non matches file
    nonmatchfile=$newnonmatchfile;
    seqsize=$(wc -c <"$nonmatchfile")
  done
done
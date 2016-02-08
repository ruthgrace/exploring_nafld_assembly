for d in ../scripts/trinity_* ; do
  if [ $d != "../scripts/trinity_CL_119_R1" ]; then
    trinityFasta="$d/Trinity.fasta";
    sample=${d:19};
    mkdir $sample;
    refseqBlastOut=$sample"/"$sample"_blast_to_refseqs.out";
    noRefseqs=$sample"/"$sample"_assembled_no_refseqs.fa";
    matchedRefseqs=$sample"/"$sample"_assembled_matched_refseqs.fa";
    seedBlast=$sample"/"$sample"_SEED_blast.out";
    # get rid of spaces in sequence identifiers
    sed -i.backup '/^>/ s/ /_/g' $trinityFasta;
    # blast with refseqs
    blastn -db all_genus_orfs_clustered_at_99_unique.fa -query $trinityFasta -out $refseqBlastOut -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -evalue 1e-3 -num_alignments 1 -num_threads 40 -perc_identity 90
    # remove refseqs
    perl ./remove_refseqs.pl $trinityFasta $noRefseqs $refseqBlastOut $matchedRefseqs;
    # blast with SEED
    blastx -db ../../../SEED_database/db_fastas.complex_plusold.faa -query $noRefseqs -out $seedBlast -outfmt 6 -evalue 1e-3 -num_alignments 10 -num_threads 40
  fi
done

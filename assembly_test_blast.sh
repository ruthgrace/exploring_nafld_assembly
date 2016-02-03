blastn -db all_genus_orfs_clustered_at_99_unique.fa -query ../scripts/trinity_CL_119_R1/Trinity.fasta -out CL_119_R1_blast_to_refseqs.out -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -evalue 1e-3 -num_alignments 1 -num_threads 10 -perc_identity 90


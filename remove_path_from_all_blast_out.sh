
for D in "$1"/*
do
  if [ -d "${D}" ]; then
    sample=${D##*/};
    echo $sample;
    blastfile=$D/$sample"_SEED_blast.out";
    blastfile_no_path=$D/$sample"_SEED_blast_no_path.out";
    ./remove_path_from_blast_out.pl $blastfile $blastfile_no_path;
  fi
done

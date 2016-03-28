dir="$1";

for d in $dir ; do
  echo $d;
  trinityFasta="$d/Trinity.fasta";
  sample=${d:19};
done

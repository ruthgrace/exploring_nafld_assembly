
for d in "$1"/*map
do
  echo $d;
  echo ${d##*/};
  sample=${d##*/};
  sample=${sample::-4};
  echo $sample;
done

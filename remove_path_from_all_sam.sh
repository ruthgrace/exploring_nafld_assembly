
for d in "$1"/*map
do
  echo $d;
  echo ${d##*/};
done

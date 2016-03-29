
for d in "$1"/*map
do
  sample=${d##*/};
  sample=${sample::-4};
  echo $sample;
  samfile=$d/$sample.sam;
  samfile_no_path=$d/$sample"_no_path.sam";
  ./remove_path_from_sam.pl $samfile $samfile_no_path;
  rm $samfile;
  mv $samfile_no_path $samfile;
done

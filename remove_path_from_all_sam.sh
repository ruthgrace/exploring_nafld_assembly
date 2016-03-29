
for d in "$1"/*map
do
  sample=${d##*/};
  sample=${sample::-4};
  echo $sample;
  samfile=$d/$sample.sam;
  samfile_no_path=$d/$sample"_no_path.sam";
  
  nohup remove_path_from_sam.sh $samfile $samfile_no_path > remove_path_from_sam_${sample}_nohup.out 2>&1&
done

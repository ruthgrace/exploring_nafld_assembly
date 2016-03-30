samfile="$1";
samfile_no_path="$2";
./remove_path_from_sam.pl $samfile $samfile_no_path;
rm $samfile;
mv $samfile_no_path $samfile;
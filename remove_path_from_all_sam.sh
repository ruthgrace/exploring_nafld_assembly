my $help = "\tFirst argument is a path to folder with mapped dirs\n";

print $help if !$ARGV[0];exit if !$ARGV[0];
my $dir = $ARGV[0];

for d in $dir ; do
  echo $d;
  trinityFasta="$d/Trinity.fasta";
  sample=${d:19};
done

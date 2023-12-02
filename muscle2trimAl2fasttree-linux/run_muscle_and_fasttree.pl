use strict;

### input region
my $infile = $ARGV[0]; ### fasta file name
my $outfile = $ARGV[1]; ### tree file name 
my $type = $ARGV[2]; ### file type in pep or cds 
####
system("muscle -in $infile -fastaout in.fasta -phyiout in.phyi");
#### 
if($type eq "cds")
{
	system("FastTree -gtr -nt in.fasta > $outfile");
}
if($type eq "pep")
{
	system("FastTree -out $outfile in.fasta");
}
####

#
#
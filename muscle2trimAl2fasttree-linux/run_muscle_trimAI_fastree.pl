use strict;

### input region
my $infile = $ARGV[0]; ### fasta file name
my $outfile = $ARGV[1]; ### tree file name 
my $type = $ARGV[2]; ### file type in pep or cds 
#### MSA
system("muscle -in $infile -fastaout in.fasta -phyiout in.phyi");
#### Triming
system("trimal -in in.fasta -out trim.fasta -htmlout trimout.html -gt 0.9 -cons 60");

# 移除列中间隙占比90%以上的列, 除非移除后剩余序列的长度小于60%
# trimal -in <inputfile> -out <outputfile> -gt 0.9 -cons 60
# 自动选择最佳的阈值来删除列
# trimal -in <inputfile> -out <outputfile> -strictplus

#### Tree
if($type eq "cds")
{
	system("FastTree -gtr -nt trim.fasta > $outfile");
}
if($type eq "pep")
{
	system("FastTree -out $outfile trim.fasta");
}
####

#
#
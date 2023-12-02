use strict;

### input region
my $infile = $ARGV[0]; ### fasta file name
my $outfile = $ARGV[1]; ### tree file name 
#my $type = $ARGV[2]; ### file type in pep or cds 
#### MSA
system("muscle -in $infile -fastaout in.fasta -phyiout in.phyi");

#### Triming
system("trimal -in in.phyi -out trim.phyi -htmlout trimout.html -gt 0.9 -cons 60");

# 移除列中间隙占比90%以上的列, 除非移除后剩余序列的长度小于60%
# trimal -in <inputfile> -out <outputfile> -gt 0.9 -cons 60
# 自动选择最佳的阈值来删除列
# trimal -in <inputfile> -out <outputfile> -strictplus

#### Make tree
#### 超快 bootstrap 参数 -bb ,推荐最小复制数 1000
#### 模型冲突的情况下，快速BS会高估BS值，推荐加上参数-bnni
#### 默认识别未完成任务，如需覆盖原任务，需加上参数-redo
#### -m MF 参数可不加，默认使用ModelFinder选模型，也可使用jModeltest或ProtTest，只需-m TEST即可
#### -nt对应的是CPU线程数，AUTO时软件自选，将AUTO改为4即设置4线程
system("iqtree -s trim.phyi -bb 1000 -bnni -redo -nt AUTO");

#### 正常的（标准非参数自举法-standard nonparametric bootstrap）BS(较慢)参数 -b,推荐最小复制数 1000
#### SH近似似然比检验(SH-like approximate likelihood ratio test) 参数 -alrt,推荐最小复制数 1000
#### 以上-b、-bb可与-alrt同时使用,会在树枝上产生两个检验值.
# system("iqtree -s $infile.phyi -b 1000 -redo -nt AUTO");
# system("iqtree -s $infile.phyi -b 1000 -alrt 1000 -redo -nt AUTO");
# system("iqtree -s $infile.phyi -bb 1000 -alrt 1000 -redo -nt AUTO");
# ####
system("cp trim.phyi.treefile $infile.phyi.treefile.nwk");
#####

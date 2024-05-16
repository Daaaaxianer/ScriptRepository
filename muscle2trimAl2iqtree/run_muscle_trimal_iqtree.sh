#! /bin/bash

### use Muscle Multiple Sequences Alignment、trimAL to Trimming and Iqtree make gene tree
### input_fasta : fasta in protein or nucleic acid
### output_tree_name : output tree file name in xxx.tree or xxx.nwk
### xian 2021-10-09 ###

if [ $# -lt 2 ];then
    echo "Two parameters needed, but only $# given!"
    echo "Usage: sh $0 <input_fasta> <output_tree_name>"
exit 1;
fi


infasta=$1
out_tree=$2 ### xxx.tree or xxx.nwk
#type=$3     ### cds or pep

### Multiple Sequences Alignment - MSA  
muscle -in ${infasta} -fastaout in.fasta -phyiout in.phyi
echo "MSA completed!"
#### Triming
trimal -in in.fasta -out trim.fasta -htmlout trimout.html -gt 0.9 -cons 60
echo "Triming completed!"
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
iqtree -s trim.phyi -bb 1000 -bnni -redo -nt AUTO

#### 正常的（标准非参数自举法-standard nonparametric bootstrap）BS(较慢)参数 -b,推荐最小复制数 1000
#### SH近似似然比检验(SH-like approximate likelihood ratio test) 参数 -alrt,推荐最小复制数 1000
#### 以上-b、-bb可与-alrt同时使用,会在树枝上产生两个检验值.
# system("iqtree -s $infile.phyi -b 1000 -redo -nt AUTO");
# system("iqtree -s $infile.phyi -b 1000 -alrt 1000 -redo -nt AUTO");
# system("iqtree -s $infile.phyi -bb 1000 -alrt 1000 -redo -nt AUTO");
# ####
mv trim.phyi.treefile ${out_tree}.treefile.nwk
echo "Outtree is : ${out_tree}"

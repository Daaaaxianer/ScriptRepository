#! /bin/bash

### use Muscle Multiple Sequences Alignment、trimAL to Trimming and Fastree make gene tree
### input_fasta : fasta in protein or nucleic acid
### output_tree_name : output tree file name in xxx.tree or xxx.nwk
### xian 2021-10-09 ###

if [ $# -lt 3 ];then
    echo "Three parameters needed, but only $# given!"
    echo "Usage: sh $0 <input_fasta> <output_tree_name> <file_type_cds_or_pep>"
exit 1;
fi


infasta=$1
out_tree=$2 ### xxx.tree or xxx.nwk
type=$3     ### cds or pep

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
if [ ${type} = "cds" ]
then
	FastTree -gtr -nt trim.fasta > ${out_tree}
fi

if [ ${type} = "pep" ]
then
	FastTree -out ${out_tree} trim.fasta
fi
#### 
echo "Outtree is : ${out_tree}"

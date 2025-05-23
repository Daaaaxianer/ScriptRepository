
#!/usr/bin/env python3

# -*- coding:UTF-8 -*-
# FileName  :del_seq.py
# Time      :2025/5/21 
# Author    :xian

from Bio import SeqIO
import argparse

def filter_sequences(cds_file, pep_file, min_length=200, exclude_genes=None, 
                    cds_output="filtered_cds.fasta", pep_output="filtered_pep.fasta"):
    """
    过滤CDS和PEP序列，根据条件移除特定基因
    
    参数:
        cds_file: 输入的CDS序列文件(FASTA格式)
        pep_file: 输入的PEP序列文件(FASTA格式)
        min_length: 最小长度阈值(默认200)，当exclude_genes=None时使用
        exclude_genes: 包含要排除基因列表的文件(单列基因名)
        cds_output: 输出的CDS序列文件名
        pep_output: 输出的PEP序列文件名
    """
    # 读取要排除的基因列表(如果提供了)
    exclude_set = set()
    if exclude_genes:
        with open(exclude_genes) as f:
            exclude_set = {line.strip() for line in f if line.strip()}
    
    # 读取CDS序列并筛选记录
    cds_records = []
    with open(cds_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # 如果有排除基因列表，只检查基因是否在排除列表中
            if exclude_genes:
                if record.id not in exclude_set:
                    cds_records.append(record.id)
            # 否则按长度过滤
            else:
                if len(record.seq) >= min_length:
                    cds_records.append(record.id)
    
    # 从PEP文件中提取对应的序列
    pep_records = []
    with open(pep_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in cds_records:
                pep_records.append(record)
    
    # 再次读取CDS文件，确保顺序一致
    filtered_cds = []
    with open(cds_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in cds_records:
                filtered_cds.append(record)
    
    # 写入输出文件
    SeqIO.write(filtered_cds, cds_output, "fasta")
    SeqIO.write(pep_records, pep_output, "fasta")
    
    print(f"筛选完成! 保留了{len(filtered_cds)}条序列。")
    if exclude_genes:
        print(f"已排除 {len(exclude_set)} 个指定基因")
    else:
        print(f"已过滤掉长度小于 {min_length} 的CDS序列")
    print(f"CDS序列已保存到: {cds_output}")
    print(f"蛋白质序列已保存到: {pep_output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='过滤CDS和PEP序列，可按长度过滤或排除指定基因')
    parser.add_argument('cds_file', help='输入的CDS序列文件(FASTA格式)')
    parser.add_argument('pep_file', help='输入的PEP序列文件(FASTA格式)')
    parser.add_argument('--min_length', type=int, default=200, 
                       help='最小长度阈值(默认200)，当不使用--exclude_genes时有效')
    parser.add_argument('--exclude_genes', 
                       help='包含要排除基因列表的文件(单列基因名)，如果使用此参数，将忽略--min_length')
    parser.add_argument('--cds_output', default="filtered_cds.fasta",
                       help='输出的CDS序列文件名(默认filtered_cds.fasta)')
    parser.add_argument('--pep_output', default="filtered_pep.fasta",
                       help='输出的PEP序列文件名(默认filtered_pep.fasta)')
    
    args = parser.parse_args()
    
    filter_sequences(args.cds_file, args.pep_file, 
                    args.min_length, args.exclude_genes,
                    args.cds_output, args.pep_output)
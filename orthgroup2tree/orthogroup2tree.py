# -*- coding: UTF-8 -*-
# FileName  : orthogroup2tree.py
# Time      : 2024/7/3
# Author    : xian

import os
import shutil
import pandas as pd
from Bio import SeqIO
import subprocess
from multiprocessing import Pool
import argparse
import numpy as np

# 设置命令行参数
parser = argparse.ArgumentParser(description='Process orthogroups and generate phylogenetic trees.')
parser.add_argument('orthogroups_file', metavar='orthogroups_file', type=str,
                    help='Path to the Orthogroups file')
parser.add_argument('-i','--indir', type=str, default='genome',
                    help='Directory containing genome files (default: genome)')
parser.add_argument('-o','--outdir', type=str, default='out',
                    help='Output directory name (default: out)')
parser.add_argument('-p', '--processes', type=int, default=5,
                    help='Number of processes to use (default: 5)')
parser.add_argument('-e', '--end', type=str, default='.pep',
                    help='Suffix for FASTA files (default: .pep)')
parser.add_argument('-n', '--min_count', type=int, default=5,
                    help='Minimum count for row comparison (default: 5)')
parser.add_argument('-x', '--max_count', type=int, default=20,
                    help='Maximum count for row comparison (default: 20)')
args = parser.parse_args()

# 解析命令行参数
genome_dir = args.indir
orthogroups_file = args.orthogroups_file
output_dir = args.outdir
num_processes = args.processes
file_suffix = args.end
min_count = args.min_count
max_count = args.max_count

# 如果输出目录存在，则删除其中所有内容
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)

# 创建输出目录
os.makedirs(output_dir)

# 步骤1：读取指定后缀的 FASTA 文件，并建立 id2seq 索引
id2seq = {}

# 遍历基因组目录下的所有指定后缀的文件
for filename in os.listdir(genome_dir):
    if filename.endswith(file_suffix):
        filepath = os.path.join(genome_dir, filename)
        for record in SeqIO.parse(filepath, 'fasta'):
            id2seq[record.id] = str(record.seq)

# 步骤2：使用 pandas 读取 Orthogroups.txt 文件，并构建数据框 orth
orth = pd.read_csv(orthogroups_file, sep='\t|\s+', index_col=0, header=None)
orth.index = orth.index.str.rstrip(':')

# 步骤3：提取序列并保存为 .fasta 文件到 out 目录
for index, row in orth.iterrows():
    if row.count() > min_count and row.count() < max_count:
        outfile = os.path.join(output_dir, f"{index}.fasta")
        with open(outfile, 'w') as f:
            for ortholog_id in row.dropna():
                if ortholog_id in id2seq:
                    f.write(f">{ortholog_id}\n{id2seq[ortholog_id]}\n")

# 定义处理函数
def process_file(filename):
    infile = os.path.join(output_dir, filename)

    # 使用 muscle 进行多序列比对
    result = subprocess.run(["muscle", "-align", infile, "-output", infile + ".align.fasta"])
    if result.returncode == 0:
        print("\nmuscle工作已完成！\n")
    else:
        print("\nmuscle工作失败！\n")

    # 使用 trimal 进行序列剪切
    result = subprocess.run(["trimal", "-in", infile + ".align.fasta",
                             "-out", infile + ".align.trimal.fasta", "-automated1"])
    if result.returncode == 0:
        print("\ntrimal工作已完成！\n")
    else:
        print("\ntrimal工作失败！\n")

    # 使用 iqtree2 构建系统发育树
    result = subprocess.run(["iqtree2", "-s", infile + ".align.trimal.fasta",
                             "-B", "1000", "-bnni", "-T", "AUTO"])
    if result.returncode == 0:
        print("\niqtree2工作已完成！\n")
    else:
        print("\niqtree2工作失败！\n")

    # 复制生成的树文件并重命名
    treefile = infile + ".align.trimal.fasta.treefile"
    if os.path.exists(treefile):
        shutil.copy(treefile, infile + ".align.trimal.fasta.treefile.nwk")
        print("\n文件复制成功！\n")
        print("\n最终文件：" + infile + ".align.trimal.fasta.treefile.nwk")
    else:
        print("\n错误：树文件不存在！")
    print("\n完成！！！")

# 步骤4：对生成的 .fasta 文件进行多序列比对、剪切和构树操作，并将结果放到 out 目录
if __name__ == '__main__':
    fasta_files = [filename for filename in os.listdir(output_dir) if filename.endswith('.fasta')]

    # 创建进程池，根据命令行参数设定最大进程数
    with Pool(processes=num_processes) as pool:
        pool.map(process_file, fasta_files)

    
# -*- coding: UTF-8 -*-
# FileName  : muscle2iqtree.py
# Time      : 2022/3/24
# Author    : xian

import argparse
import subprocess
import shutil

def run_muscle(infile):
    subprocess.run(["muscle", "-align", infile, "-output", infile + ".align.fasta"])
    print("\nmuscle work completed!\n")

def run_iqtree(infile):
    subprocess.run(["iqtree2", "-s", infile + ".align.fasta", "-B", "1000", "-bnni", "-redo", "-T", "AUTO"])
    print("\niqtree work completed!\n")

def copy_treefile(infile):
    shutil.copy(infile + ".align.fasta.treefile", infile + ".align.fasta.treefile.nwk")
    print("\nFinal file: " + infile + ".align.fasta.treefile.nwk")

def main(infile):
    # Run muscle
    run_muscle(infile)
    
    # Run iqtree2
    run_iqtree(infile)

    # Copy treefile
    copy_treefile(infile)
    print("\nDone!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run muscle and iqtree2 on a fasta file.")
    parser.add_argument("infile", type=str, help="Input fasta file name")
    args = parser.parse_args()
    main(args.infile)

#### 脚本使用 muscle5 + iqtree2,支持win和linux
#### 超快 bootstrap 检验: 参数-B (v1.x版的等价参数为-bb) ,推荐最小复制数 1000
#### 模型冲突的情况下，快速BS会高估BS值，推荐加上参数-bnni
#### 标准 bootstrap 检验: 参数-b，推荐最小复制数 100，比快速要慢很多
#### 类SH近似似然比检验(SH-like approximate likelihood ratio test): 参数-alrt，推荐最小复制数 1000
#### 以上-B、-b可与-alrt同时使用,会在树枝上产生两个检验值。
#### 默认识别未完成任务，如需覆盖原任务，需加上参数-redo
#### 模型选择：-m MF 运行 ModelFinder 找到最合适的模型不做树重建
#### 模型选择：-m MFP 运行 ModelFinder Plus 找到最合适的模型，后续用此模型构树等
#### -m MFP 模型选择参数可不加，默认使用-m MFP选模型，也可使用-m 指定模型
#### -T (v1.x版的等价参数为-nt)对应的是CPU线程数，AUTO时软件自选，将AUTO改为4即设置4线程
#### -ntmax 设置最大调用内核数
#### 移除列中间隙占比90%以上的列(-gt), 除非移除后剩余序列的长度小于60%(-cons)
#### trimal -in <inputfile> -out <outputfile> -gt 0.9 -cons 60 
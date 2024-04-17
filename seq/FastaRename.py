# -*- coding:utf-8 -*-
# @FileName :FastaRename.py
# @Time     :2021/11/11 19:25
# @Author   :Xian

import re
from Bio.Seq import Seq
from Bio import SeqIO
import argparse

def OldFa2NewFa(infa,ingff,ty='fasta'):
    Oid2Nid = {}
    with open(ingff, 'r+') as gff:
        for line in gff:
            line = line.strip()
            if re.match('^#', line): continue
            a = line.split("\t")
            # print(a[6],a[1])
            Oid2Nid[a[6]] = a[1]
    SeqList = []
    for seq_record in SeqIO.parse(infa,ty):
        # print(seq_record.description)
        b = re.split('[\=\s+]',seq_record.description) ## 每次均需要核对
        selectId = b[0]  ## 每次均需要核对
        if selectId not in Oid2Nid.keys(): 
            print( "Wrong select : " + selectId )
            continue
        print (seq_record.id,Oid2Nid[selectId])
        seq_record.id = Oid2Nid[selectId]
        seq_record.description = selectId
        SeqList.append(seq_record)
    return SeqList
if __name__ == "__main__":
    # 定义命令行解析器对象
    parser = argparse.ArgumentParser(description='Modify fasta sequence according to gff')

    # 添加命令行参数
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")
    parser.add_argument("fasta", type=str, help="fasta sequence filepath ")
    parser.add_argument("gff", type=str, help="gff filepath ")
    parser.add_argument("-o", "--outfasta", type=str, default='out.fasta',
                        help="output file name (default: out.fasta)")
    parser.add_argument("-t", "--type", type=str, default='fasta', help="output file type (default: fasta)")

    # 解析命令行参数
    args = parser.parse_args()

    # 执行
    seqList = OldFa2NewFa(args.fasta, args.gff)
    SeqIO.write(seqList, args.outfasta, args.type)

    # CdsList = OldFa2NewFa("Vvinifera_145_Genoscope.12X.cds.fa","Vvi.new.gff")
    # SeqIO.write(CdsList,"Vvi.new.cds.fa","fasta")
    # PepList = OldFa2NewFa("Vvinifera_145_Genoscope.12X.protein.fa", "Vvi.new.gff")
    # SeqIO.write(PepList, "Vvi.new.pep.fa", "fasta")


# -*- coding:utf-8 -*-
# @FileName :FastaRename.py
# @Time     :2021/11/11 19:25
# @Author   :Xian

import re
from Bio.Seq import Seq
from Bio import SeqIO

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
        print (seq_record.id,Oid2Nid[seq_record.id])
        seq_record.id = Oid2Nid[seq_record.id]
        seq_record.description = seq_record.id
        SeqList.append(seq_record)
    return SeqList
if __name__ == "__main__":
    CdsList = OldFa2NewFa("Vvinifera_145_Genoscope.12X.cds.fa","Vvi.new.gff")
    SeqIO.write(CdsList,"Vvi.new.cds.fa","fasta")
    PepList = OldFa2NewFa("Vvinifera_145_Genoscope.12X.protein.fa", "Vvi.new.gff")
    SeqIO.write(PepList, "Vvi.new.pep.fa", "fasta")


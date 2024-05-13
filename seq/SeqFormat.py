# -*- coding:UTF-8 -*-
# FileName  :seqFormat.py
# Time      :2021/3/30
# Author    :xian

import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from Bio.SeqUtils import GC123
import sys
import argparse

def cdsTrans(cdsSeq, pepSeq):
    '''
    :param cdsSeq: 输入cds序列文件名
    :param pepSeq: 输出氨基酸序列文件名
    :return: 返回氨基酸序列文件
    '''
    with open(pepSeq, 'w') as outfasta:
        for seq_record in SeqIO.parse(cdsSeq, "fasta"):
            print(seq_record.id)
            # print(seq_record.seq)
            # print(repr(seq_record.seq))
            # print(len(seq_record))
            seq_record.seq = seq_record.seq.translate(table=1)
            # 定义密码子表使用table参数（https: // www.ncbi.nlm.nih.gov / Taxonomy / Utils / wprintgc.cgi）
            # print(seq_record.seq)
            SeqIO.write(seq_record, outfasta, "fasta")

def gcFilter(naSeq, gcTable):
    '''
    :param naSeq: 输入核酸序列文件名
    :param gcTable: 输出序列gc含量列表名
    :return: 返回gc含量列表
    '''
    with open(gcTable, 'w') as outtable:
        outtable.write("seq_id\tGCtotal\tGC1\tGC2\tGC3 \n")
        for seq_record in SeqIO.parse(naSeq, "fasta"):
            # print(seq_record.id, GC(seq_record.seq)) ### 计算每条序列GC含量
            # print(seq_record.id, GC123(seq_record.seq))  ### 计算每条序列total,first,second,third位的GC含量,存储于元组tuple中
            gcCount = repr(GC123(seq_record.seq)[0]) + "\t" + repr(GC123(seq_record.seq)[1]) + "\t" + repr(
                GC123(seq_record.seq)[2]) + "\t" + repr(GC123(seq_record.seq)[3])
            print(seq_record.id,gcCount)
            # print(seq_record.seq)
            # print(repr(seq_record.seq))
            outlist = seq_record.id
            outtable.write(seq_record.id +"\t"+ gcCount +"\n")

def seqDict(seq, form = "fasta"):
    '''
    用途：将给定的生物序列放入字典中（内存存储，文件过大可使用Bio.SeqIO.index()
    用法：用一个变量接收即可得到id为key，序列SeqRecord为value的的字典
    SeqRecord 包含id,name,seq,description等属性，支持Bio.seq相关操作(https://biopython-cn.readthedocs.io/zh_CN/latest/cn/chr04.html)
    :param seq:字符串格式序列文件名
    :param form:字符串格式序列格式如 fasta、genebank等
    :return:序列id为key，序列的SeqRecord对象为value
    '''
    #n = seq.find('.')
    #seq_dict = seq[0:seq.find('.')]
    seq_dict = SeqIO.to_dict(SeqIO.parse(seq, str(form)))
    #return seq_dict.keys()
    print("\nSequence dictionary has been created!\n")
    return seq_dict

def readId(idfile):
    with open(idfile, 'r') as ids:
        content = ids.read()
        linelist = content.splitlines()
    print("\nThe id file has been read!\n")
    return linelist

def seqExtract(allSeq, idfile, extractSeq):
    '''
    :param allSeq: 包含所有序列的文件，默认fasta格式
    :param idfile: 待提取的id文件，所有id单列排列
    :param extractSeq: 提取到的序列文件，默认fasta格式
    :return: 提取到的序列文件
    '''
    seq_dict = seqDict(allSeq)
    idList = readId(idfile)
    seq_records = []
    for id in idList:
        print(seq_dict[id].id)
        seq_records.append(seq_dict[id])
    SeqIO.write(seq_records, extractSeq, "fasta")

def readCutfile(cutfile):
    with open(cutfile, 'r') as cutfs:
        content = cutfs.read()
        cutlist = []
        linelist = content.splitlines()
        for line in linelist:
            line.strip()
            cutarr = re.split("\t", line)
            print(cutarr[0])
            cutlist.append(cutarr)
    print("\nThe cut file has been read!\n")
    return cutlist

def cutSeq(allSeq, cutfile, cutSeq):
    '''
    :param allSeq: 包含所有序列的文件，默认fasta格式
    :param cutfile: 包含切取信息的文件，依次为 原id、起始位点、终止位点、新id，以table分隔
    :param cutSeq: 切取到序列文件，默认fasta文件
    :return: 切取到的序列文件
    '''
    seq_dict = seqDict(allSeq)
    cutList = readCutfile(cutfile)
    seq_records = []
    for cut in cutList:
        # print(cut[0])
        print(seq_dict[cut[0]].id)
        # print(seq_dict[cut[0]].seq)
        cutStart = int(cut[1]) - 1
        cutEnd = int(cut[2])
        # print(type(cutEnd))
        my_seq = seq_dict[cut[0]].seq

        # 起始位点大于终止位点，需要反向截取
        if cutStart > cutEnd:
            cutStart = int(cut[2]) - 1
            cutEnd = int(cut[1])
            my_seq = my_seq[::-1] ## 序列取反,非互补

        # 终止位点过长
        if cutEnd > len(seq_dict[cut[0]].seq):
            print(f'\n\tWARNING: The value of CutEnd is greater than the length of {cut[0]}!\n')
            break

        # 构造 SeqRecord
        my_record = SeqRecord(Seq(my_seq),
                              id = cut[3],
                              description = "\"Subseq of " + cut[0] + "\""
                              )
        my_record.seq = my_record.seq[cutStart:cutEnd]
        seq_records.append(my_record)
    SeqIO.write(seq_records, cutSeq, "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "用途：操作核苷酸和氨基酸序列")
    parser.add_argument("--cdstrans", action = "store_true",
                        help = "Translating nucleic acid coding sequnece into amino acid sequence(need -c -p)")
    parser.add_argument("-c", "--cdsSeq", type = str, help = "Input coding sequence")
    parser.add_argument("-p", "--pepSeq", type = str, default = "out.pep.fasta",
                        help = "Output amino acid sequence (default: out.pep.fasta)")

    parser.add_argument("--gc", action = "store_true", help = "Calculate the gc content of nucleic acid sequences(need -n -g)")
    parser.add_argument("-n", "--naSeq", type = str, help = "Input nucleic acid sequence")
    parser.add_argument("-g", "--gcTable", type = str, default = "out.gctable.txt",
                        help = "Summary table of gc content (default: out.gctable.txt).")

    parser.add_argument("--extract", action = "store_true", help = "Extract the sequence of a given id(need -s -i -e)")
    parser.add_argument("--cut", action="store_true",
                        help="Cut the sequence according to the given location information(need -s -l -m)")
    parser.add_argument("-s", "--Seq", type = str, help = "Input nucleic acid or amino acid sequence")
    parser.add_argument("-i", "--Id", type = str, help = "Imported id list in text file")
    parser.add_argument("-e", "--extractedSeq",type = str, default = "out.extractedSeq.fasta",
                        help = "Extracted sequence file (default: out.extractedSeq.fasta)." )
    parser.add_argument("-l", "--Loc", type=str, help="Imported cut location in text file")
    parser.add_argument("-m", "--miniSeq", type=str, default="out.cutSeq.fasta",
                        help="Cutted sequence file (default: out.cutSeq.fasta).")

    args = parser.parse_args()
    if args.cdstrans:
        cdsTrans(args.cdsSeq, args.pepSeq)
    if args.gc:
        gcFilter(args.naSeq, args.gcTable)
    if args.extract:
        seqExtract(args.Seq, args.Id, args.extractedSeq)
    if args.cut:
        cutSeq(args.Seq, args.Loc, args.miniSeq)



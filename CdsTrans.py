# -*- coding:utf-8 -*-
'''
用途：翻译cds序列使用标准密码子表,定义密码子表使用table参数（https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi）
用法： python CdsTrans.py [cdsfile] [proteinfile]
作者：xian
日期：2020年4月24日
'''

from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import GC123
import sys

def get_accession(record):
    """
    "Given a SeqRecord, return the accession number as a string.
    利用管道符拆分，提取第四部分并返回,可作为SeqIO.to_dict的参数
    e.g. "gi|2765613|emb|Z78488.1|PTZ78488" -> "Z78488.1"
    """
    parts = record.id.split("|")
    assert len(parts) == 5 and parts[0] == "gi" and parts[2] == "emb"
    return parts[3]

def SeqDict(seq,form):
    '''
    用途：将给定的生物序列放入字典中（内存存储，文件过大可使用Bio.SeqIO.index()）
    用法：用一个变量接收即可得到id为key的字典，对应value为序列record格式，支持Bio.seq相关操作
    :param seq:字符串格式序列文件名
    :param form:字符串格式序列格式如 fasta、genebank等
    :return:
    '''
    #n = seq.find('.')
    #seq_dict = seq[0:seq.find('.')]
    seq_dict = SeqIO.to_dict(SeqIO.parse(seq, str(form)))
    #return seq_dict.keys()
    return seq_dict

if __name__ == '__main__':
    # sl_dict =SeqDict('sl.4.0.cds.fa','fasta')
    # print(sl_dict.keys())
    with open(sys.argv[2], 'w') as outfasta:
        for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
            #print(seq_record.id,GC(seq_record.seq)) ### 计算每条序列GC含量
            #print(seq_record.id, GC123(seq_record.seq))  ### 计算每条序列total,first,second,third位的GC含量
            # print(seq_record.seq)
            # print(repr(seq_record.seq))
            seq_record.seq = seq_record.seq.translate(table=1)
            SeqIO.write(seq_record, outfasta, "fasta")


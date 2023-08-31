#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: xian
@version: 1.0.0
@license: Apache Licence
@file: ete3tree.py
@time: 2023/6/11 2:24
"""

from ete3 import PhyloTree
from Bio import AlignIO
from ete3 import TreeStyle

# 读取树文件
tree = PhyloTree("InR_tree.nwk")

# 读取比对文件
# alignment = Alignment("InR_alignment.fa")

# 读取序列文件
# sequences = SeqIO.parse("sequences.fasta", "fasta")

# 将比对结果添加到树上
tree.link_to_alignment(alignment = "InR_alignment.fa",  alg_format="fasta")

# # 将序列标签添加到树上
# for node in tree.iter_descendants("postorder"):
#     if not node.is_leaf():
#         continue
#     seq_id = node.name
#     sequence = next(sequences)
#     node.add_features(sequence=sequence.seq)

# 定义一个TreeStyle对象来控制树的显示方式
ts = TreeStyle()

# 显示树和比对结果
tree.render("tree_and_alignment.pdf", tree_style=ts)


# # Load a tree and link it to an alignment.
# t = PhyloTree("InR_tree.nwk")
# print(t)
#
# t.link_to_alignment(alignment="InR_alignment.fa", alg_format="fasta")
# print(t)

if __name__ == '__main__':
    pass

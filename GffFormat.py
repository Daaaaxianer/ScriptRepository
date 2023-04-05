# -*- coding:utf-8 -*-
# @FileName :GffFormat.py
# @Time     :2021/11/4 22:48
# @Author   :Xian

import re

# try:
#     infile = open('Vvinifera_145_Genoscope.12X.gene.gff3', 'r')
#     print(infile.read())
#     print(infile.readlines()) ### readlines()等价于 for i in infile:得到的列表
#     print(infile.readline())
# finally:
#     if f:
#         f.close()

###
with open('Vvinifera_145_Genoscope.12X.gene.gff3', 'r') as gff3:
    GffList = gff3.readlines()

SelList = []
for m in GffList:
    line = m.split("\t")
    if re.match('^#',line[0]):
        continue
    if re.match('(?!mRNA)', line[2]):
        continue
    if re.search(r'\D+$',line[0]):
        line[0] = "0"
    line[0] = re.sub(r'chr','',line[0])
    group = re.split(r'=|;|\.',line[8])
    # print(line[0]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[6]+"\t"+group[1])
    LineList = [line[0],line[3],line[4],line[6],group[1]]
    SelList.append(LineList)

NewList= sorted(SelList, key = lambda x: (int(x[0]),int(x[1])),reverse = False) ## 按1、2列升序排列

ChrNum  = 0
GeneNum = 0
GeneBp = 0
ChrOrder = {}
ChrBp = {}
with open('Vvi.new.gff', 'w+') as gff:
    for i in range(len(NewList)):
    # for i in NewList:
        if(int(NewList[i][0]) > ChrNum):
            ChrNum = int(NewList[i][0])
            GeneNum = 1
        else:
            GeneNum = GeneNum + 1
        ChrBp[ChrNum] = int(NewList[i][2])
        ChrOrder[ChrNum] = GeneNum
        GeneO = "{:0>5d}".format(GeneNum) #### 格式化
        NewId = "VVi"+NewList[i][0]+"g"+ GeneO
        NewList[i].insert(1,NewId)
        NewList[i].insert(5,str(GeneNum))
        gff.write("\t".join(NewList[i])+"\n")

with open('Vvi.lens', 'w+') as lens:
    for k in ChrOrder.keys():
        lens.write(str(k) +"\t"+str(ChrBp[k])+"\t"+ str(ChrOrder[k]) +"\n")
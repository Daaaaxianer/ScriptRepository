#_*_ coding:utf-8 _*_
import re

infile = open('Vvinifera_145_Genoscope.12X.gene.gff3', 'r')
outfile1 = open('new.gff','w')
outfile2 = open('len.lens','w')

G2C={} ### old gene id to its chromosome
G2C["1"]="vv"+"1"
G2S={} ### old gene id to its startpoint
G2E={} ### old gene id to its endpoint
G2D={} ### old gene id to its direction
G2N={} ### old gene id to its newid
G2O={} ### old gene id to its order
C2L={} ### chromosome to its length
oids=[]  ### list of old gene id

for line in infile:	
	if re.match('^#', line): #### 去除##开头的行
		continue
	line = line.strip() #### strip 默认去除字符串两端\n \r \t 和' '
	a = line.split("\t") #### tab拆分当前行, 并存入列表
	if re.match('(?!mRNA)', a[2]): #### match 不是mRNA返回None
		continue
	if (a[0][-1] == 'm') or (a[0][-1] == 'n'):
		continue
	# print(a) ### 输出a列表看一下格式
	# ['chr1', 'phytozomev10', 'mRNA', '7323830', '7327622', '.', '+', '.', listb]
	chrab = a[0][3:] #### 染色体编号
	startpoint = a[3] #### 起始
	endpoint = a[4] #### 终止
	direction = a[6] #### 方向

	b = re.split('=|:|;',a[8]) #### python 内置split只能单分隔符拆分，故借用re.split分割当前行最后一列
	# print(b) ### 输出b列表看一下格式
	# ['ID', 'GSVIVT01000534001.Genoscope12X', 'Name', 'GSVIVT01000534001', ...]

	#### 提取出下面5个排序
	oid = b[1] #### 旧ID
	oids.append(oid)  #### 旧ID
	G2C[oid] = "" + chrab #### 染色体编号
	G2S[oid] = startpoint #### 起始
	G2E[oid] = endpoint #### 终止
	G2D[oid] = direction #### 方向
	# print ('chromosome is:'+ G2C[oid] +'\nstartpoint is:'+ G2S[oid] + '\nendpoint is:'+ G2E[oid] +'\ndirection is:' + G2D[oid] +'\n') ### 输出旧id,起始终止位点和转录方向
	# break ### 输出一个, 检查是否有问题

#### 对旧名字列表按G2S和G2C进行排序,注意优先级问题

#sortoids = sorted(oids,key = lambda x:(int(G2S[x])), reverse = False)
#sortoids = sorted(sortoids,key = lambda x:(G2C[x]), reverse = False)
#### 上例为匿名函数lambda,也可自定义函数排序

def by_chr(x):
	return G2C[x]

def by_startpoint(x):
	return int(G2S[x])

## sort / sorted
sortoids = sorted(oids,key = by_startpoint , reverse = False) #### 排序结果存入sortoids列表
sortoids = sorted(sortoids,key = by_chr , reverse = False) #### False 递增

#### 排序后生成新名字 写入新gff
#### 注意总结遍历列表的几种方法

lastchr = '' ### 上一个染色体的编号
num = 0 ### default count number

for i, val in enumerate(sortoids):
	# print('oldid:' + val + '\nlocats in chrosome:' + G2C[val] + '\nwith startpoint:' + G2S[val] + '\n')
	# break
	# GSVIVT01012261001.Genoscope12X vv1 10731

	if G2C[val] != lastchr:
		print('>>> New chromosome' + G2C[val] + '...\n')

		#### chromosome length
		if G2C[val][-2:] != 00 and lastchr != '':
			C2L[lastchr[-2:]]=num
		num = 0
		num += 1 ### num = num + 1
		G2O[val] = num



	if G2C[val] == lastchr:
		print('>>> Chromosome' + G2C[val] + '...\n')
		num += 1   ### c += a 等效于 c = c + a
		G2O[val] = num

	numlen = len(str(num))
	m = 4 - numlen
	zero = "0" * m
	nid = G2C[val] + 'g' + zero + str(num)
	G2N[val] = nid

	lastchr = G2C[val] ### marked lastchr

	outfile1.write(G2C[val]+"\tVvi"+G2N[val]+"\t"+G2S[val]+"\t"+G2E[val]+"\t"+G2D[val]+"\t"+str(G2O[val])+"\t"+val+"\n")

C2L[G2C[val][-2:]] = num


for key in C2L:
	print(key+"\t"+str(C2L[key])+"\n")
	outfile2.write(key+"\t"+str(C2L[key])+"\n")

infile.close
outfile1.close
outfile2.close
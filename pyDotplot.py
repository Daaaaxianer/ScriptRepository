#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: xian
@version: 1.0.0
@license: GPL Licence
@file: pyDotplot.py
@time: 2023/4/22 22:00
"""

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
from matplotlib.patches import Rectangle
import pandas as pd
import numpy as np
import argparse

# mplstyle.use('fast')

## 输出 matplotlibrc 配置文件路径
print(mpl.matplotlib_fname())

# 查看内置的ttf字体
# a = sorted([f.name for f in mpl.font_manager.fontManager.ttflist])
# for i in a:
#     print(i)

# 指定默认字体
# sci 一般要求使用无衬线字体：Arial或Helvetica
# 支持中文的字体："Songti SC", "Wawati TC", "STHeiti", "Arial Unicode MS"
mpl.rcParams['font.sans-serif'] = ['Arial']
mpl.rcParams['axes.unicode_minus'] = False # 解决保存图像是负号’-'显示为方块的问题

## 定义 Figure 基础参数
# figsize 默认使用 inch 作为单位,可以使用 Figure.dpi 属性来访问图形的 DPI，默认值为 100。
figWidth = 12
figHeight = 10
# fig, ax = plt.subplots(1, 1, figsize=(figWidth, figHeight), layout='constrained')
fig, ax = plt.subplots(1, 1, figsize=(figWidth, figHeight))
ax.invert_yaxis() # 将y轴翻转
# fig.set_tight_layout('tight') ## 紧凑布局
# fig.set_tight_layout({"pad": 0.2})  ## 手动设置标准化布局参数
plt.subplots_adjust(left=0.02, right=1.02, top=0.98, bottom=-0.02) ## 手动调整边距

ax.set_axis_off()  # 关闭轴的可见性，不显示轴线、刻度和标签
# ax.set_xticks([])  # 设置x轴刻度为空列表，不显示x轴刻度标签
# ax.set_yticks([])  # 设置y轴刻度为空列表，不显示y轴刻度标签

# ax.set_title("Grape vs. Tomato genome dotplot", fontsize=25, fontweight=800, y=0.98)
# 设置图形的标题，标题的字体大小为25，字体粗细为800，标题的垂直位置为0.99（0到1之间）

## 设置pandas显示选项，不省略中间的列
pd.set_option("display.max_columns", None)

# 读入 lens 文件，构建 DataFrame
# chrnum chr_bp chr_order total_bp total_order
def readLens(file):
    """
    读取标准化的 lens 文件存储为 dataframe
    lens文本文件应包含如下列：
    chrnum chr_bp chr_order total_bp total_order cum_bp cum_order

    :param file: 读入lens文件名
    :return: chrnum为index的dataframe
    """
    # lens—>DataFrame ，其中chrnum为index行名
    lens = pd.read_csv(file, sep="\t", header=None, names=["chrnum", "chr_bp", "chr_order"])

    # 约束数据类型
    lens['chrnum'] = lens['chrnum'].astype(str)
    lens['chr_bp'] = lens['chr_bp'].astype(int)
    lens['chr_order'] = lens['chr_order'].astype(int)

    # 增加两列分别为bp和order的总和
    lens['total_bp'] = lens['chr_bp'].sum()
    lens['total_order'] = lens['chr_order'].sum()

    # 设置chrnum为行索引
    lens.set_index('chrnum', inplace=False)

    # 增加累积的cum_bp和cum_order列 (cumulative bp和cumulative order) ,可选是否包含当前行的bp或order
    lens['cum_bp'] = lens['chr_bp'].cumsum()
    # lens['cum_bp'] = lens['chr_bp'].cumsum() - lens['chr_bp']
    lens['cum_order'] = lens['chr_order'].cumsum()
    # lens['cum_order'] = lens['chr_order'].cumsum() - lens['chr_order']

    return lens

def readGff(gff):
    """
    读取标准化的 gff 文件存储为 dataframe

    读入 gff文本文件应包含如下列：
    chr id initial end direction order oldid

    存储的 dataframe 包含如下列：
    id chr initial end direction order oldid average

    :param gff: 读入gff文件名
    :return: id 为index，增加 平均bp(average)列 的新 dataframe
    """
    gff = pd.read_csv(gff, sep="\t", header=None, names=["chr", "id", "initial", "end", "direction", "order", "oldid"])
    gff['chr'] = gff['chr'].astype(str)
    gff['id'] = gff['id'].astype(str)
    gff['initial'] = gff['initial'].astype(int)
    gff['end'] = gff['end'].astype(int)
    gff['direction'] = gff['direction'].astype(str)
    gff['order'] = gff['order'].astype(int)
    gff['oldid'] = gff['oldid'].astype(str)

    # 增加平均bp位置列
    # 计算 gff DataFrame 中 “initial” 和 “end” 两列数据的平均值，并将该平均值存储为gff dataframe 中的 average 列
    average = gff[["initial", "end"]].apply(np.mean, axis=1)
    gff['average'] = round(average).astype(int)

    #设置id为行索引
    gff.set_index('id', inplace=True)
    return gff


def read_blast_outfmt6(file_path):
    """
    读取 NCBI-blast生成的 outfmt6文件，存储为 dataframe

    读入 blast文件 应为 >blast-2.10+ 输出的 outfmt6 文本文件，包含如下列：
    query_id id subject_id percent_identity alignment_length mismatches gap_opens query_start query_start query_end subject_start subject_end evalue bit_score

    存储的 dataframe 包含如下列：
    query_id id subject_id percent_identity alignment_length mismatches gap_opens query_start query_start query_end subject_start subject_end evalue bit_score

    :param file_path: blast 文本文件 in outfmt 6
    :return:
    """
    blast_columns = ["query_id", "subject_id", "percent_identity", "alignment_length",
                     "mismatches", "gap_opens", "query_start", "query_end", "subject_start",
                     "subject_end", "evalue", "bit_score"]
    colums_type = ["str", "str", "float", "int64", "int64", "int64", "int64", "int64", "int64", "int64", "float", "float"]

    type_dict = dict(zip(blast_columns, colums_type))

    # 读入blast文件，指定列名及其数据格式
    blast_df = pd.read_csv(file_path, sep="\t", header=None, names=blast_columns,
                           dtype= type_dict)
    return blast_df


def process_blast_results(blast_df, vgff, hgff):
    '''
    处理 NCBI Blast 结果文件
    存储的 dataframe 包含如下列：
    query_id subject_id percent_identity alignment_length mismatches gap_opens query_start query_start query_end subject_start subject_end evalue bit_score
     count total_count dot_color qchr schr

    :param blast_df: 读入的 NCBI Blast 数据框
    :return: 处理后的 NCBI Blast DataFrame
    '''
    # 新增“出现次数”列
    # cumcount(): 对于每个组中的元素，返回元素在该组中的序列号（从0开始）。
    # 例如，如果有三条query_id都为q1的数据，则它们在分组后的DataFrame中分别具有序号0、1和2。
    # 把序列号加1，方便从1开始计数。
    blast_df['count'] = blast_df.groupby('query_id').cumcount() + 1

    # 新增“总出现次数”列
    # 这里使用了map的 Series映射 进行 数据替换
    total_count = blast_df['query_id'].value_counts()
    blast_df['total_count'] = blast_df['query_id'].map(total_count)

    # 新增dot_color列
    # 创建数字与颜色的映射字典，1映射为red，2映射为blue，3-5均映射为gray，其他映射为空。
    # 这里使用了map的 字典dict映射 进行 数据替换
    # 三色方案：#bb3b2c, #435888, #8aa4b5
    # color_dict = {1: "#F44336", 2: "#1E88E5", 3: "lightgray", 4: "lightgray", 5: "lightgray", 6: "lightgray"}
    # color_dict = {1: "#E53935", 2: "#207dff", 3: "lightgray", 4: "lightgray", 5: "lightgray", 6: "lightgray"}
    color_dict = {1: "#bb3b2c", 2: "#207dff", 3: "lightgray", 4: "lightgray", 5: "lightgray", 6: "lightgray"}
    blast_df['dot_color'] = blast_df['count'].map(color_dict)

    # 基于gff 在blast中新增qchr、schr、qorder、sorder、qaverage、saverage列
    # 使用了map的 Series映射 进行 数据替换
    qchr_dict = vgff['chr']
    schr_dict = hgff['chr']
    qorder_dict = vgff['order']
    sorder_dict = hgff['order']
    qaverage_dict = vgff['average']
    saverage_dict = hgff['average']

    blast_df['qchr'] = blast_df['query_id'].map(qchr_dict).astype("str")
    blast_df['schr'] = blast_df['subject_id'].map(schr_dict).astype("str")
    blast_df['qorder'] = blast_df['query_id'].map(qorder_dict).fillna(0).astype("int64")
    blast_df['sorder'] = blast_df['subject_id'].map(sorder_dict).fillna(0).astype("int64")
    blast_df['qaverage'] = blast_df['query_id'].map(qaverage_dict).fillna(0).astype("int64")
    blast_df['saverage'] = blast_df['subject_id'].map(saverage_dict).fillna(0).astype("int64")

    # 删除按顺序超过第6次出现的行query_id
    blast_df = blast_df[blast_df['count'] <= 3]

    return blast_df

def draw_dotplot(vsp, vgff, vlens, hsp, hgff, hlens, blast, pos="order"):
    '''
    绘制 genome dotplot
    :param vsp: 垂直 Vertical 排列的 物种名称
    :param vgff: 垂直 Vertical 排列的 gff
    :param vlens: 垂直 Vertical 排列的 lens
    :param hsp: 水平 Horizontal 排列的 物种名称
    :param hgff: 水平 Horizontal 排列的 gff
    :param hlens: 水平 Horizontal 排列的 lens
    :param blast: ncbi Blast+ outfmt6
    :param pos: 位置表示方式默认为order，可选：order 或 bp
    :return:
    '''

    # 读取lens文件
    ver_lens = readLens(vlens)
    hor_lens = readLens(hlens)

    # 读取gff文件
    ver_gff = readGff(vgff)
    hor_gff = readGff(hgff)

    # 读取ncbi blast
    blast_df = read_blast_outfmt6(blast)
    print(blast_df.info())
    print("# " * 30)
    blast_df = process_blast_results(blast_df, ver_gff, hor_gff)

    # 基于lens 在blast中新增qcumbp、scumbp、qcumorder、scumorder列
    # 使用了map的 Series映射 进行 数据替换
    qcumbp_dict = pd.Series(ver_lens['cum_bp'].values - ver_lens['chr_bp'].values, index=ver_lens['chrnum'].values)
    scumbp_dict = pd.Series(hor_lens['cum_bp'].values - hor_lens['chr_bp'].values, index=hor_lens['chrnum'].values)

    qcumorder_dict = pd.Series(ver_lens['cum_order'].values - ver_lens['chr_order'].values, index=ver_lens['chrnum'].values)
    scumorder_dict = pd.Series(hor_lens['cum_order'].values - hor_lens['chr_order'].values, index=hor_lens['chrnum'].values)

    blast_df['qcum_bp'] = blast_df['qchr'].map(qcumbp_dict).fillna(0).astype("int64")
    blast_df['scum_bp'] = blast_df['schr'].map(scumbp_dict).fillna(0).astype("int64")

    blast_df['qcum_order'] = blast_df['qchr'].map(qcumorder_dict).fillna(0).astype("float64")
    blast_df['scum_order'] = blast_df['schr'].map(scumorder_dict).fillna(0).astype("float64")


    # 计算垂直和水平的单位量 并 计算染色体边界
    if pos == "order":
        ver_total = ver_lens['chr_order'].sum()
        hor_total = hor_lens['chr_order'].sum()
        verUnit = figHeight / ver_total
        horUnit = figWidth / hor_total

        ver_lens["yLoc"] = ver_lens["cum_order"] * verUnit
        hor_lens["xLoc"] = hor_lens["cum_order"] * horUnit

        # 计算query 和 subject 的位置
        blast_df['qloc'] = blast_df.apply(lambda row: (row['qcum_order'] + row['qorder']) * verUnit, axis=1)
        blast_df['sloc'] = blast_df.apply(lambda row: (row['scum_order'] + row['sorder']) * horUnit, axis=1)

    elif pos == "bp":
        ver_total = ver_lens['chr_bp'].sum()
        hor_total = hor_lens['chr_bp'].sum()
        verUnit = figHeight / ver_total
        horUnit = figWidth / hor_total

        ver_lens["yLoc"] = ver_lens["cum_bp"] * verUnit
        hor_lens["xLoc"] = hor_lens["cum_bp"] * horUnit

        # 计算query 和 subject 的位置
        blast_df['qloc'] = blast_df.apply(lambda row: (row['qcum_bp'] + row['qaverage']) * verUnit, axis=1)
        blast_df['sloc'] = blast_df.apply(lambda row: (row['scum_bp'] + row['saverage']) * horUnit, axis=1)

    # 删除或保留blast_df中的 lens不包含的染色体基因
    blast_df = blast_df[blast_df['qchr'].isin(ver_lens['chrnum'])]
    # # blast_df = blast_df[~blast_df['qchr'].isin(ver_lens['chrnum'])]
    blast_df = blast_df[blast_df['schr'].isin(hor_lens['chrnum'])]
    # # blast_df = blast_df[~blast_df['schr'].isin(hor_lens['chrnum'])]

    print("lens " + "$" * 30)
    print(ver_lens.info())
    # print(hor_lens.info())
    print("lens " + "$" * 30)
    print("\n")

    print("gff " + "$" * 30)
    print(ver_gff.info())
    # print(hor_gff.head(10))
    print("gff " + "$" * 30)
    print("\n")

    print("blast " + "$" * 30)
    # print(blast_df.head(10))
    print(blast_df.info())
    print("blast " + "$" * 30)
    print("\n")

    # 绘制染色体分隔线
    verlist = ver_lens["yLoc"].tolist()
    horlist = hor_lens["xLoc"].tolist()
    # verlist.insert(0, 0)
    # horlist.insert(0, 0)

    ax.hlines(verlist, 0, figWidth, lw=0.5, color='k', linestyle='-', alpha=0.8)
    ax.vlines(horlist, 0, figHeight, lw=0.5, color='k', linestyle='-',alpha=0.8)

    # 根据原始 Series 创建等长的元素全为 -0.1 的 Series 对象
    ver_x = pd.Series([-0.1] * len(ver_lens["yLoc"]))
    hor_y = pd.Series([-0.1] * len(hor_lens["xLoc"]))

    # 分别绘制垂直和水平方向染色体编号的文本
    for i in range(len(ver_lens["chrnum"])):
        ax.text(ver_x[i], ver_lens["yLoc"][i] - 0.5 * ver_lens["chr_order"][i] * verUnit , ver_lens["chrnum"][i], fontsize= 14, rotation=0, ha='right', va='center')

    for i in range(len(hor_lens["chrnum"])):
        ax.text(hor_lens["xLoc"][i] - 0.5 * hor_lens["chr_order"][i] * horUnit, hor_y[i]-0.15, hor_lens["chrnum"][i], fontsize= 14, rotation=0, ha='center', va='top')

    # 分别绘制垂直和水平方向的物种名称
    ax.text(-0.80, 0.5 * verUnit * ver_total, vsp, fontsize= 24, fontweight= 600, rotation= 90, ha= 'left', va= 'center')
    ax.text(0.5 * horUnit * hor_total, -0.45, hsp, fontsize= 24, fontweight= 600, rotation= 0, ha= 'center', va= 'center')

    # 绘制相似性的点
    # 按相似性由低到高绘制
    ax.scatter(blast_df[blast_df["count"] >= 3]['sloc'], blast_df[blast_df["count"] >= 3]['qloc'],  s=0.8, facecolors=blast_df[blast_df["count"] >= 3]['dot_color'], edgecolors='none')
    ax.scatter(blast_df[blast_df["count"] == 2]['sloc'], blast_df[blast_df["count"] == 2]['qloc'],  s=0.8, facecolors=blast_df[blast_df["count"] == 2]['dot_color'], edgecolors='none')
    ax.scatter(blast_df[blast_df["count"] == 1]['sloc'], blast_df[blast_df["count"] == 1]['qloc'],  s=0.8, facecolors=blast_df[blast_df["count"] == 1]['dot_color'], edgecolors='none')

    # ax.scatter(blast_df[blast_df.count == 2]['sloc'], blast_df[blast_df.count == 2]['qloc'],  s=1.0, facecolors=blast_df[blast_df.count == 2]['dot_color'], edgecolors='none')
    # ax.scatter(blast_df[blast_df.count == 1]['sloc'], blast_df[blast_df.count == 1]['qloc'],  s=1.0, facecolors=blast_df[blast_df.count == 1]['dot_color'], edgecolors='none')

    # ax.scatter(blast_df['sloc'], blast_df['qloc'],  s=1.0, facecolors=blast_df['dot_color'], edgecolors='none')

    # 绘制矩形外框
    # 创建一个矩形对象
    rectangle = Rectangle((0, 0), horlist[-1], verlist[-1], edgecolor='black', facecolor='none', linewidth=1.0)
    # 将矩形对象添加到图形中
    ax.add_patch(rectangle)

    # 将数据框写入csv文件
    blast_df.to_csv('processed_blast.csv', index=False)


if __name__ == "__main__":
    # 定义命令行解析器对象
    parser = argparse.ArgumentParser(description='Visualization Program of Genome Dotplot')
    # 添加命令行参数
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")
    parser.add_argument("query_name", type=str, help="Query specie name in vertical order")
    parser.add_argument("subject_name", type=str, help="Subject specie names in vertical order")
    parser.add_argument("query_gff", type=str, help="GFF file of query specie")
    parser.add_argument("subject_gff", type=str, help="GFF file of subject specie")
    parser.add_argument("query_lens", type=str, help="LENS file of query specie")
    parser.add_argument("subject_lens", type=str, help="LENS file of subject specie")
    parser.add_argument("blastout", type=str, help="Ncbi blast result file in outfmt6 format")
    parser.add_argument("-p", "--pos", type=str, default='order',choices=['order','bp'], help="Optional location type: 'order' or 'bp'(default: order) ")
    parser.add_argument("-o", "--output", type=str, default='Genome.dotplot.png', help="Figure name with format suffix.  Png, jpg, svg, pdf, eps, tiff, bmp and other matplotlib formats are supported (default: Genome.dotplot.png)")

    # 解析命令行参数
    args = parser.parse_args()

    # 执行绘图
    draw_dotplot(args.query_name, args.query_gff, args.query_lens, args.subject_name, args.subject_gff, args.subject_lens, args.blastout, args.pos)
    # draw_dotplot("Grape", "vv.12x.gff", "vv.12x.lens", "Tomato", "sl.4.0.all.gff", "sl.4.0.all.lens", "vv.sl.pep.1e-5.mts20.outfmt6",pos="order")

    # 保存图像
    # fig.savefig('grape.vs.tomato.genome.dotplot.png', dpi=300)
    fig.savefig(args.output, dpi=300)
    # plt.clf() #清除当前图像
    plt.show()
    sys.exit(0)
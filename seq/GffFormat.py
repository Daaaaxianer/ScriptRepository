# -*- coding:UTF-8 -*-
# FileName  :GffFormat.py
# Time      :2023/9/5 
# Author    :yuxian

import pandas as pd
import argparse
import numpy as np
import re
# from SeqFormat import readId

## 设置pandas显示宽度，不省略中间的列
# pd.set_option("display.max_columns", None)
pd.set_option('display.max_columns', 30)
# pd.set_option('display.max_rows', 30)
pd.set_option('display.max_colwidth', 200)
pd.set_option('display.width', 300)

def process_gff3(gff3_file):
    '''
    读取GFF3文件并存储为DataFrame
    其中提取最长转录本代表该位点基因
    :param gff3_file: 标准GFF3文件
    :return: 仅剩最长转录本的DataFrame
    '''
    # 1. 读取GFF3文件并存储为DataFrame
    gff3_df = pd.read_csv(gff3_file, sep="\t", header=None,
                          names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "group"])

    # 2. 使用布尔索引筛选出不以"#"开头的行并删除
    gff3_df = gff3_df[~gff3_df['seqname'].str.startswith("#")]

    # 3. 删除某些无效染色体seqname的行
    values_to_remove = ['ChrSy', 'ChrUn', 'chrUn']

    gff3_df = gff3_df[~gff3_df['seqname'].isin(values_to_remove)]
    gff3_df = gff3_df[~(gff3_df['seqname'].str.endswith(('_random',)))] # 可根据实际情况增减

    # gff3_df = gff3_df[~gff3_df['seqname'].isin(['ChrSy', 'ChrUn', 'chrUn'])]

    # 4. 只保留feature为"mRNA"的行，并重置索引
    gff3_df = gff3_df[gff3_df['feature'] == 'mRNA'].reset_index(drop=True)

    # 5. 从group列 name字段中 提取 gene loc 和转录本 name
    # gff3_df['name'] = gff3_df['group'].str.extract(r'([\w.]+)\.')
    gff3_df['loc'] = gff3_df['group'].str.split("Parent=").str[1].str.split(".").str[0] ## 每次均需要核对
    gff3_df['name'] = gff3_df['group'].str.split("Name=").str[1].str.split(";").str[0] ## 每次均需要核对

    # 6. 比较每个gene的转录本的长度差异，仅保留长度差异最大的转录本
    gff3_df['length_diff'] = gff3_df['end'] - gff3_df['start']
    gff3_df = gff3_df.loc[gff3_df.groupby('loc')['length_diff'].idxmax()].reset_index(drop=True)

    # 将start和end列的类型设置为整数
    gff3_df['start'] = gff3_df['start'].astype(int)
    gff3_df['end'] = gff3_df['end'].astype(int)
    gff3_df['length_diff'] = gff3_df['length_diff'].astype(int)

    print(gff3_df)

    return gff3_df

def process_dataframe(df, letters):
    '''
    将GFF3文件转换为标准GFF格式
    :param df: 输入GFF3
    :param letters: 缩写字母
    :return: 标准gff格式
    '''

    # 提取seqname列末尾的数字并创建新列chr
    df['chr'] = df['seqname'].str.extract(r'(\d+)$').astype(int)

    # 按照chr和start列进行升序排序
    df = df.sort_values(by=['chr', 'start'], ascending=[True, True])

    # 新增order列，根据chr进行分组后的行计数
    df['order'] = df.groupby('chr').cumcount() + 1

    # 新增lensbp列和lensorder列
    df['lensbp'] = df.groupby('chr')['end'].transform('max')
    df['lensorder'] = df.groupby('chr')['order'].transform('max')

    # 新增id列，自定义三个字母 + chr列的值 + g + order列的数值组成的字符串
    # df['id'] = (df['order'].apply(lambda x: f"{x:05d}")
    #             .apply(lambda x: f"{letters}{df['chr'].values[0]}g{x}"))

    df['id'] = (letters + df['chr'].astype(str) + 'g' + df['order'].astype(str).str.zfill(5))

    # 分组chr与其对应的lensbp和lensorder值输出到单独的文件
    lens_mapping = df[['chr', 'lensbp', 'lensorder']].drop_duplicates().reset_index(drop=True)
    lens_mapping.to_csv(f"{letters}.lens.txt", sep='\t', header=False, index=False)

    # 修改列顺序为chr, id, start, end, strand, order, name
    df = df[['chr', 'id', 'start', 'end', 'strand', 'order', 'name']]

    # 输出数据框到文件
    df.to_csv(f"{letters}.gff", header=False, sep='\t', index=False)

    return df

def id2bed(idfile, gff_file, bed_file):
    '''
    :param idfile: 待提取id文件，单列
    :param gff_file: 经 GffFormat.py gff3 处理后的 自定义Gff文件。
    :return: bed file for MG2C (http://mg2c.iask.in/mg2c_v2.1/)
    '''
    idList = readId(idfile)
    gffDf = readGff(gff_file)
    bed_df = gffDf.loc[gffDf["Id"].isin(idList), ["Oid","start","end","chrId"]]
    # 输出数据框到文件
    bed_df.to_csv(f"{bed_file}", header=False, sep='\t', index=False)

def readGff(gff_file):
    '''
    :param gff_file: custom gff file from gff3
    :return: gff_df
    '''
    # 读取cunstom gff文件并存储为DataFrame
    gff_df = pd.read_csv(gff_file, sep="\t", header=None,
                          names=["chrId", "Id", "start", "end", "strand", "order", "Oid"])
    return gff_df

def readId(idfile):
    with open(idfile, 'r') as ids:
        content = ids.read()
        linelist = content.splitlines()
    print("\nThe id file has been read!\n")
    return linelist


if __name__ == '__main__':
    # 定义命令行解析器对象
    parser = argparse.ArgumentParser(description='Flexible conversion involving gff3, gff, bed')

    # 添加命令行参数
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")

    # 功能1: Convert standard gff3 to custom gff
    parser.add_argument("--gff3gff", action="store_true",
                        help="Convert standard gff3 to custom gff format(need -g -s)")
    parser.add_argument("-g", "--gff3", type=str, help="Normal filepath of GFF3")
    parser.add_argument("-s", "--species", type=str, default='species',
                        help="Optional Species abbreviation (default: species)")

    # 功能2: Read id file and extract bed for MG2C from custom gff
    parser.add_argument("--id2bed", action="store_true",
                        help="extract bed for MG2C from custom gff(need -i -b)")
    parser.add_argument("-i", "--id", type=str,help="Id file to be extracted")
    parser.add_argument("-f", "--gff", type=str, help="custom gff file")
    parser.add_argument("-b", "--bed", type=str, default='out.bed', help="Extracted bed file")


    # 解析命令行参数
    args = parser.parse_args()

    # 执行操作
    if args.gff3gff:
        gff3_pd = process_gff3(args.gff3)
        gff_df = process_dataframe(gff3_pd, args.species)
        # print(gff_df)

    if args.id2bed:
        id2bed(args.id, args.gff, args.bed)

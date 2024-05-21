# -*- coding:UTF-8 -*-
# FileName  :gff_filter.py
# Time      :2024/5/21 
# Author    :xian

import pandas as pd
import argparse
import numpy as np
import re
import os
import shutil


## 设置pandas显示宽度，不省略中间的列
# pd.set_option("display.max_columns", None)
pd.set_option('display.max_columns', 30)
# pd.set_option('display.max_rows', 30)
pd.set_option('display.max_colwidth', 200)
pd.set_option('display.width', 300)


def readGff(gff_file):
    '''
    :param gff_file: custom gff file from gff3 :"chrId", "Id", "start", "end", "strand", "order", "Oid"
    :return: gff_df
    '''
    gff_df = pd.read_csv(gff_file, sep='\t', header=None, names=["chrId", "Id", "start", "end", "strand", "order", "Oid"])

    print(gff_df)
    return gff_df

def readGffFilter(gffFilter):
    '''
    :param gffFilter: gff filter file :Chr  Start   End Name
    :return: gffFilter_df
    '''
    gff_filter_df = pd.read_csv(gffFilter, sep='\t', header=None, names=["Chr", "Start", "End", "Name"])
    
    print(gff_filter_df)
    return gff_filter_df

def process_gff_data(gff_df, gff_filter_df):
    output_directory = "output_ids"
    
    # 删除已存在的输出目录并重新创建
    if os.path.exists(output_directory):
        shutil.rmtree(output_directory)
    
    os.makedirs(output_directory)

    # 按照条件筛选出符合要求的 Id
    selected_ids = []
    for _, row_filter in gff_filter_df.iterrows():
        ids = gff_df[(gff_df['chrId'] == row_filter['Chr']) & (gff_df['start'] >= row_filter['Start']) & (gff_df['end'] <= row_filter['End'])]['Id'].tolist()
        selected_ids.append(','.join(ids))

    gff_filter_df['selId'] = selected_ids

    print(gff_filter_df)

    # 将中间df输出到文件
    gff_filter_df.to_csv("temp.csv", index=False)

    # 按照Name进行分组
    grouped = gff_filter_df.groupby('Name')

    for name, group in grouped:
        # 检查 selId 列是否为空，并且是否有非空白字符串
        if not group['selId'].empty and any(sel_id.strip() for sel_id in group['selId']):
            all_ids = []
            # 遍历每个分组，提取并拆分 selId 列中的 IDs
            for sel_ids in group['selId']:
                if sel_ids.strip():  # 检查 selId 是否存在任何字符串
                    all_ids.extend(sel_ids.split(','))

            # 输出每个 Name 分组的 ID 到单独的文本文件中
            with open(f"{output_directory}/{name}.txt", 'w') as file:
                file.write("\n".join(all_ids))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process GFF data.')
    parser.add_argument('-gff', help='Path to the GFF file')
    parser.add_argument('-fgff', help='Path to the GFF filter file')

    args = parser.parse_args()

    if args.gff and args.fgff:
        gff_df = readGff(args.gff)
        gff_filter_df = readGffFilter(args.fgff)
        process_gff_data(gff_df, gff_filter_df)
    else:
        print("Please provide both GFF file and GFF filter file paths.")

# -*- coding:UTF-8 -*-
# FileName  :te_classification_and_counting.py
# Time      :2024/12/15 
# Author    :yuxian

import re
import argparse
import os


def count_class_family(file_path):
    # 初始化一个空字典来存储class/family及其出现次数
    class_family_counter = {}

    # 计数器，用于跳过前4行
    line_count = 0

    with open(file_path, 'r') as file:
        for line in file:
            line_count += 1
            # 跳过前4行
            if line_count <= 4:
                continue

            # 跳过以特定字符（例如#）开头的任何标题或注释行
            if line.startswith('#') or not line.strip():
                continue

            # 使用正则表达式根据一个或多个空白字符分割行成列
            columns = re.split(r'\s+', line.strip())

            # 定义class/family列的索引（根据实际文件格式调整）
            class_family_idx = 10

            try:
                # 从行中提取class/family信息
                class_family = columns[class_family_idx]

                # 输出当前行的class/family信息
                # print(f"Line {line_count}: Class/Family: {class_family}")

                # 更新该class/family的计数器
                if class_family in class_family_counter:
                    class_family_counter[class_family] += 1
                else:
                    class_family_counter[class_family] = 1
            except IndexError:
                print(f"警告: 第 {line_count} 行中的列数不符合预期。跳过该行。")
                continue

    return class_family_counter


def write_sorted_counts_to_file(counter, output_file_path):
    # 按class/family字母顺序对计数器进行排序
    sorted_items = sorted(counter.items())

    # 将排序后的项写入输出文件
    with open(output_file_path, 'w') as output_file:
        for class_family, count in sorted_items:
            output_file.write(f"{class_family}\t{count}\n")


def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="统计并输出te文件中每个class/family的数量")
    parser.add_argument('input_file', help='输入的te文件路径')
    parser.add_argument('-o', '--output', help='输出文件路径 (默认为输入文件路径+“_classification_and_counting.txt”)')

    args = parser.parse_args()

    input_file_path = args.input_file
    if args.output:
        output_file_path = args.output
    else:
        # 如果没有提供输出文件名，则默认在输入文件的末尾附加"_classification_and_counting.txt"
        output_file_path = f"{os.path.splitext(input_file_path)[0]}_classification_and_counting.txt"

    # 统计class/family出现次数，并将排序后的计数写入输出文件
    counter = count_class_family(input_file_path)
    write_sorted_counts_to_file(counter, output_file_path)

    print(f"统计结果已写入 {output_file_path}")


if __name__ == '__main__':
    main()
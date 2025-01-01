# -*- coding: utf-8 -*-
# FileName  : 2.codonBias_plot.py
# Time      : 2024/12/31 01:29
# Author    : xian

import argparse
from CodonU import vizualizer as viz
import os
import logging

def main(nuc_file, prot_file, organism_name, output_dir, genetic_table_num=1):
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    try:
        # 打开核酸序列文件和蛋白质序列文件
        with open(nuc_file, 'r') as nuc_handle, open(prot_file, 'r') as prot_handle:
            # 1.绘制 ENc Plot
            logging.info(f"Generating ENc Plot for {organism_name}...")
            viz.plot_enc(
                handle=nuc_handle,
                genetic_table_num=genetic_table_num,
                organism_name=organism_name,
                save_image=True,
                folder_path=output_dir
            )
            logging.info(f"ENc Plot saved to: {output_dir}")
            nuc_handle.seek(0)  # 重置核酸文件指针到文件开头

            # 2.绘制 Neutrality Plot
            logging.info(f"Generating Neutrality Plot for {organism_name}...")
            viz.plot_neutrality(
                handle=nuc_handle,
                min_len_threshold=200,
                organism_name=organism_name,
                save_image=True,
                folder_path=output_dir
            )
            logging.info(f"Neutrality Plot saved to: {output_dir}")
            nuc_handle.seek(0)  # 重置核酸文件指针到文件开头

            # 3.绘制 Parity Rule 2 Plot
            logging.info(f"Generating Parity Rule 2 Plot for {organism_name}...")
            viz.plot_pr2(
                handle=nuc_handle,
                min_len_threshold=200,
                organism_name=organism_name,
                save_image=True,
                folder_path=output_dir
            )
            logging.info(f"Parity Rule 2 Plot saved to: {output_dir}")
            nuc_handle.seek(0)  # 重置核酸文件指针到文件开头

            # 4.绘制 Codon Count 分析 (每个密码子)
            logging.info(f"Generating Codon Count analysis for each codon for {organism_name}...")
            viz.plot_ca_codon_freq_codon(
                handle=nuc_handle,
                genetic_table_num=genetic_table_num,
                organism_name=organism_name,
                save_image=True,
                folder_path=output_dir
            )
            logging.info(f"Codon Count analysis for each codon saved to: {output_dir}")
            nuc_handle.seek(0)  # 重置核酸文件指针到文件开头

            # # 5.绘制 Codon Count 分析 (每个基因) ### need repair
            # logging.info(f"Generating Codon Count analysis for each gene for {organism_name}...")
            # viz.plot_ca_codon_freq_gene(
            #     handle=nuc_handle,
            #     genetic_table_num=genetic_table_num,
            #     organism_name=organism_name,
            #     save_image=True,
            #     folder_path=output_dir
            # )
            # logging.info(f"Codon Count analysis for each gene saved to: {output_dir}")
            # nuc_handle.seek(0)  # 重置核酸文件指针到文件开头

            # 6.绘制 Codon RSCU 分析 (每个密码子)
            logging.info(f"Generating Codon RSCU analysis for each codon for {organism_name}...")
            viz.plot_ca_codon_rscu_codon(
                handle=nuc_handle,
                genetic_table_num=genetic_table_num,
                organism_name=organism_name,
                save_image=True,
                folder_path=output_dir
            )
            logging.info(f"Codon RSCU analysis for each codon saved to: {output_dir}")
            nuc_handle.seek(0)  # 重置核酸文件指针到文件开头

            # # 7.绘制 Codon RSCU 分析 (每个基因) ## need repair
            # logging.info(f"Generating Codon RSCU analysis for each gene for {organism_name}...")
            # viz.plot_ca_codon_rscu_gene(
            #     handle=nuc_handle,
            #     genetic_table_num=genetic_table_num,
            #     organism_name=organism_name,
            #     save_image=True,
            #     folder_path=output_dir
            # )
            # logging.info(f"Codon RSCU analysis for each gene saved to: {output_dir}")

            # 8.绘制 AA Count 分析 (每个氨基酸)
            prot_handle.seek(0)  # 重置蛋白质文件指针到文件开头
            logging.info(f"Generating AA Count analysis for each AA for {organism_name}...")
            viz.plot_ca_aa_freq_aa(
                handle=prot_handle,
                genetic_table_num=genetic_table_num,
                organism_name=organism_name,
                save_image=True,
                folder_path=output_dir
            )
            logging.info(f"AA Count analysis for each AA saved to: {output_dir}")

            # # 9.绘制 AA Count 分析 (每个基因) 使用 'GRAVY' ### need repair
            # nuc_handle.seek(0)  # 重置核酸文件指针到文件开头
            # logging.info(f"Generating AA Count analysis for each AA using 'GRAVY' for {organism_name}...")
            # viz.plot_ca_aa_freq_gene(
            #     handle=nuc_handle,
            #     genetic_table_num=genetic_table_num,
            #     scale='gravy',
            #     organism_name=organism_name,
            #     save_image=True,
            #     folder_path=output_dir
            # )
            # logging.info(f"AA Count analysis for each gene using 'GRAVY' saved to: {output_dir}")
            #
            # # 10.绘制 AA Count 分析 (每个基因) 使用 'Aromaticity' ## need repair
            # nuc_handle.seek(0)  # 重置核酸文件指针到文件开头
            # logging.info(f"Generating AA Count analysis for each gene using 'Aromaticity' for {organism_name}...")
            # viz.plot_ca_aa_freq_gene(
            #     handle=nuc_handle,
            #     genetic_table_num=genetic_table_num,
            #     scale='aroma',
            #     organism_name=organism_name,
            #     save_image=True,
            #     folder_path=output_dir
            # )
            # logging.info(f"AA Count analysis for each gene using 'Aromaticity' saved to: {output_dir}")

        print("所有图表已成功生成并保存到指定目录。")

    except Exception as e:
        logging.error(f"Error during analysis: {e}")
        print(f"发生错误: {e}")


if __name__ == "__main__":
    # 设置日志记录
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    parser = argparse.ArgumentParser(description="绘制多个生物信息学分析图")

    # 添加位置参数
    parser.add_argument("nuc_file", help="核酸序列文件路径")
    parser.add_argument("prot_file", help="蛋白质序列文件路径")
    parser.add_argument("organism_name", help="物种名称")
    parser.add_argument("output_dir", help="输出图像保存的文件夹路径")

    # 添加可选参数
    parser.add_argument("-g", "--genetic_table_num", type=int, default=1, help="遗传密码表编号，默认为1")

    args = parser.parse_args()

    main(args.nuc_file, args.prot_file, args.organism_name, args.output_dir, args.genetic_table_num)
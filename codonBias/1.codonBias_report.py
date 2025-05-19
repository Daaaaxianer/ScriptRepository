# -*- coding:UTF-8 -*-
# FileName  : 1.codonBias_report.py
# Time      : 2024/12/31 00:35
# Author    : xian


import argparse
from CodonU import analyzer as an
import pandas as pd
from Bio.SeqIO import parse, write
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import logging
import glob
import shutil  

# 设置日志配置
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("analysis.log"),
        logging.StreamHandler()
    ]
)

def process_sequences(nuc_file, prot_file, output_dir, genetic_code_num=1):
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 获取文件名前缀
    nuc_prefix = os.path.splitext(os.path.basename(nuc_file))[0] if nuc_file else None
    prot_prefix = os.path.splitext(os.path.basename(prot_file))[0] if prot_file else None

    # 清理旧的报告文件
    report_dir = os.path.join(output_dir, 'Report')
    clean_old_reports(report_dir)

    # 处理核苷酸序列
    if nuc_file:
        logging.info(f"Processing nucleotide sequences from {nuc_file}...")
        calculate_nucleotide_metrics(nuc_file, nuc_prefix, genetic_code_num, output_dir)

    # 处理蛋白质序列
    if prot_file:
        logging.info(f"Processing protein sequences from {prot_file}...")
        calculate_protein_metrics(prot_file, prot_prefix, genetic_code_num, output_dir)

def clean_old_reports(report_dir):
    """删除旧的报告文件"""
    try:
        # 确保报告目录存在
        os.makedirs(report_dir, exist_ok=True)

        # 查找并删除已有的报告文件 (例如 report_*.txt)
        report_pattern = os.path.join(report_dir, 'report_*.txt')
        existing_reports = glob.glob(report_pattern)
        if existing_reports:
            for report in existing_reports:
                try:
                    os.remove(report)  # 删除旧的报告文件
                    logging.info(f"Deleted existing report file: {report}")
                except Exception as e:
                    logging.error(f"Failed to delete {report}: {e}")
    except Exception as e:
        logging.error(f"Error cleaning old reports: {e}")

def calculate_nucleotide_metrics(in_file, file_prefix, genetic_code_num, output_dir):
    try:
        # 计算RSCU
        rscu = an.calculate_rscu(in_file, genetic_code_num, gene_analysis=True, save_file=False)
        rscu_df = pd.DataFrame(rscu).T
        rscu_csv_path = os.path.join(output_dir, f"{file_prefix}_RSCU.csv")
        rscu_df.to_csv(rscu_csv_path, index_label="Codon")
        logging.info(f"The RSCU score file can be found at: {rscu_csv_path}")

        # 计算CAI
        cai = an.calculate_cai(in_file, genetic_code_num, min_len_threshold=200, gene_analysis=True, save_file=False)
        cai_df = pd.DataFrame(cai).T
        cai_csv_path = os.path.join(output_dir, f"{file_prefix}_CAI.csv")
        cai_df.to_csv(cai_csv_path, index_label="Gene")
        logging.info(f"The CAI score file can be found at: {cai_csv_path}")

        # 计算CBI
        cbi = an.calculate_cbi(in_file, genetic_code_num, gene_analysis=True, save_file=False)
        cbi_df = pd.DataFrame(cbi).T
        cbi_csv_path = os.path.join(output_dir, f"{file_prefix}_CBI.csv")
        cbi_df.to_csv(cbi_csv_path, index_label="Gene")
        logging.info(f"The CBI score file can be found at: {cbi_csv_path}")

        # 计算ENc
        enc = an.calculate_enc(in_file, genetic_code_num, gene_analysis=True, save_file=False)
        enc_df = pd.DataFrame(list(enc.items()), columns=["Name", "Value"])
        enc_df.set_index("Name",inplace=True)
        enc_csv_path = os.path.join(output_dir, f"{file_prefix}_ENc.csv")
        enc_df.to_csv(enc_csv_path, index_label="Gene")
        logging.info(f"The ENc score file can be found at: {enc_csv_path}")

        # 生成核酸分析报告
        report_dir = os.path.join(output_dir, 'Report')
        os.makedirs(report_dir, exist_ok=True)  # 确保报告目录存在
        an.generate_report(in_file, 'nuc', genetic_code_num, min_len_threshold=200, res_folder_path=report_dir)
        logging.info(f"Nucleotide analysis report generated in the '{report_dir}' folder.")

    except Exception as e:
        logging.error(f"Error processing nucleotide sequences: {e}")

def calculate_protein_metrics(in_file, file_prefix, genetic_code_num, output_dir):
    try:
        # 定义一个函数来移除蛋白质序列末尾的终止密码子
        def remove_stop_codon(record):
            seq = str(record.seq)
            if seq.endswith('*'):
                seq = seq[:-1]
            return SeqRecord(Seq(seq), id=record.id, description=record.description)

        # 读取并处理蛋白质序列
        with open(in_file, 'r') as f:
            filtered_records = [remove_stop_codon(record) for record in parse(f, "fasta")]
            if not filtered_records:
                logging.warning(f"No valid protein sequences found in {in_file}.")
                return

        # 创建临时文件以保存处理后的序列
        temp_fasta_path = os.path.join(output_dir, f"{file_prefix}_temp.fasta")
        with open(temp_fasta_path, 'w') as temp_fasta:
            write(filtered_records, temp_fasta, "fasta")

        # 计算Aromaticity
        aroma = an.calculate_aromaticity(temp_fasta_path, gene_analysis=True, save_file=False, folder_path=output_dir)
        aroma_df = pd.DataFrame(list(aroma.items()), columns=["Protein_name", "Aroma_score"])
        aroma_df.set_index("Protein_name", inplace=True)
        aroma_csv_path = os.path.join(output_dir, f"{file_prefix}_Aromaticity.csv")
        aroma_df.to_csv(aroma_csv_path)
        logging.info(f"The Aromaticity score file can be found at: {aroma_csv_path}")
        # logging.info(f"The Aromaticity score file can be found at: {os.path.join(output_dir, f'{file_prefix}_Aromaticity.xlsx')}")

        # 计算GRAVY
        gravy = an.calculate_gravy(temp_fasta_path, gene_analysis=True, save_file=False, folder_path=output_dir)
        gravy_df = pd.DataFrame(list(gravy.items()), columns=["Protein_name", "Gravy_score"])
        gravy_df.set_index("Protein_name", inplace=True)
        gravy_csv_path = os.path.join(output_dir, f"{file_prefix}_GRAVY.csv")
        gravy_df.to_csv(gravy_csv_path)
        logging.info(f"The GRAVY score file can be found at: {gravy_csv_path}")
        # logging.info(f"The GRAVY score file can be found at: {os.path.join(output_dir, f'{file_prefix}_GRAVY.xlsx')}")

        # 生成蛋白分析报告
        report_dir = os.path.join(output_dir, 'Report')
        os.makedirs(report_dir, exist_ok=True)  # 确保报告目录存在
        an.generate_report(temp_fasta_path, 'aa', genetic_code_num, min_len_threshold=66, res_folder_path=report_dir)
        logging.info(f"Protein analysis report generated in the '{report_dir}' folder.")

        # 删除临时文件
        os.remove(temp_fasta_path)
        logging.info(f"Temporary file {temp_fasta_path} has been deleted.")

    except Exception as e:
        logging.error(f"Error processing protein sequences: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Codon usage analysis.")
    parser.add_argument("nucleotide_file", help="Path to the nucleotide FASTA file")
    parser.add_argument("protein_file", help="Path to the protein FASTA file")
    parser.add_argument("output_directory", help="Path to the output directory")
    parser.add_argument("-g","--genetic_table_num", type=int, default=1, help="Genetic code number (default: 1)")

    args = parser.parse_args()

    # 将 'None' 字符串转换为 None 类型
    nuc_fasta = args.nucleotide_file if args.nucleotide_file.lower() != 'none' else None
    prot_fasta = args.protein_file if args.protein_file.lower() != 'none' else None

    process_sequences(nuc_fasta, prot_fasta, args.output_directory, args.genetic_table_num)

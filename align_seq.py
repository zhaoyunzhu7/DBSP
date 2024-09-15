from Bio import AlignIO
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from collections import Counter

def align_seqs(input_fasta, file_path):
    # 读取包含多个序列的fasta文件
    # input_fasta = "D:\\Python\\PythonProject\\clu_ass_graph\\before_align_sequences.fa"
    sequences = SeqIO.parse(input_fasta, "fasta")

    # 创建一个ClustalW多序列比对进程
    clustalw_cline = ClustalwCommandline("D:\\MSA\\clustalw\\clustalw2.exe", infile=input_fasta)
    clustalw_cline()

    # 读取比对后的结果
    alignment = AlignIO.read(file_path, "clustal")
    # alignment2 = AlignIO.read("D:\\Python\\PythonProject\\clu_ass_graph\\result\\msa_sequences_{:.2f}.aln".format(error_rate), "clustal")

    # 初始化一个空字符串来存储最终的正确序列
    correct_sequence = ""

    # 循环遍历每个位置
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]  # 获取当前位置的碱基列
        counts = Counter(column)  # 统计当前位置不同碱基的出现次数
        most_common_base = counts.most_common(1)[0][0]  # 获取出现次数最多的碱基
        correct_sequence += most_common_base

    # 输出最长序列的碱基序列（去除破折号）
    correct_sequence = str(correct_sequence).replace("-", "")
    # print("Correct Sequence:")
    # print(correct_sequence)
    return correct_sequence


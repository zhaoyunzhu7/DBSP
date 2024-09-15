import os
from collections import defaultdict
from graph import build_dbg, compute_merged_sequence, read_fasta
from sequence_alignment_diagram import *
from align_seq import align_seqs
from recovery import data_recovery_rate



def comparision_seq(ss_list, G, k):
    """
    :param ss_list: 骨干序列的k-mer列表
    :param G: 所有序列的DBG图
    :param k: k-mer的参数
    :return: 返回一条组装后的序列
    """
    container = []  # container存储DBG中选择路径的节点信息，即一个一个k-mer
    container.append("CCTGCAGAGTA")
    i = 0  # 指示骨干序列中每个k-mer的位置
    # drift = 0  # 骨干序列上允许偏移的跨度
    drift = [0]
    while len(container) < 151 - k + 1:
        # print("iiiiiiiiiiiiiiiiiiiiiiiiii")
        try:
            current_node = container[-1]  # 当前节点为路径的最后一个节点
            next_nodes = list(G.successors(current_node))
            # print("current_node:", current_node,"next_nodes:",next_nodes)
            # print("##################################",len(next_nodes))
            # print()
            if len(next_nodes) == 0:
                # 删除container的最后一个元素
                container.pop()
                while len(next_nodes) == 0:
                    container.pop()
                    i -= 1
                current_node = container[-1]  # 当前节点为路径的最后一个节点
                next_nodes = list(G.successors(current_node))
                sorted_neighbors = sorted(next_nodes, key=lambda x: G[current_node][x]['frequency'],
                                          reverse=True)
                container.append(sorted_neighbors[0])
                # container.append(ss_list[i + 1])
                # continue
            elif len(next_nodes) == 1:
                next_node = next_nodes[0]  # 将当前节点的下一个节点添加到container中，并将其设为新的当前节点
                container.append(next_node)
                i += 1

            elif len(next_nodes) > 1:
                flag = False
                # for j in range(-drift, drift+1):
                for j in drift:
                    if 0 <= i + j < len(ss_list) and current_node == ss_list[i + j]:
                        # print("current_node:", current_node, "j:", j,"i:", i)
                        if 0 <= i + j + 1 < len(ss_list) and ss_list[i + j + 1] in next_nodes:
                            # print("next_node:", ss_list[i + j + 1])
                            container.append(ss_list[i + j + 1])
                            i += 1
                            flag = True
                            break
                        else:
                            sorted_neighbors = sorted(next_nodes, key=lambda x: G[current_node][x]['frequency'],
                                                      reverse=True)
                            container.append(sorted_neighbors[0])
                            # container.append(ss_list[i + 1])
                            i += 1
                            flag = True
                            break
                    else:
                        # print(i)
                        continue
                if not flag:
                    # print('xxx')
                    sorted_neighbors = sorted(next_nodes, key=lambda x: G[current_node][x]['frequency'],
                                              reverse=True)
                    container.append(sorted_neighbors[0])
                    # container.append(ss_list[i + 1])
                    i += 1

                continue

        except Exception as e:
            # print("error:", e)
            container = ss_list
            break

    # print(container)
    ss_list_merge = ss_list[0]
    for i in ss_list[1:]:
        ss_list_merge += i[-1]
    # print("聚类后的序列：", ss_list_merge)
    # 合并前后缀并还原为原始DNA序列
    assem_sequence = container[0]
    for segment in container[1:]:
        assem_sequence += segment[-1]
    # print("组装后的序列：", assem_sequence)
    # print("container",len(container))
    # ss_list

    return assem_sequence


if __name__ == "__main__":
    error_rates = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17,
                   0.18, 0.19, 0.2]
    number_of_droplets = 1  # 选择前多少条
    reps = 1  # 文件数目（实验次数）
    yes_file_sequences = read_fasta("D:\\Python\\PythonProject\\clu_ass_graph\\yes.fa")
    # 数据恢复率，只算一条的(yes.fa中只有一条)
    yes_file_sequences_string = ''.join(yes_file_sequences)
    k = 11
    assem_sequences_all = []  # List to store all assem sequences
    all_fp_MSA_sequence_recovery_rate = 0
    all_BSDBG_recovery_rate = 0
    all_BSDBG_no_beam_recovery_rate = 0
    average_BSDBG_recovery_rate = 0
    average_fp_MSA_sequence_recovery_rate = 0
    average_BSDBG_no_beam_recovery_rate = 0
    beam_width = 10
    for error_rate in error_rates:
        directory = f"D:\\Python\\PythonProject\\clu_ass_graph\\data\\error{error_rate}"
        files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
        for filename in files:
            if filename.endswith(".fa"):
                print("文件名：", filename)
                # 组合文件路径
                filepath = os.path.join(directory, filename)
                filepath_aln = "D:\\Python\\PythonProject\\clu_ass_graph\\data\\error{}\\mix_eDNAs_all_error{}_iteration1.aln".format(error_rate, error_rate)
                filepath_msa_sequence = align_seqs(filepath, filepath_aln)

                # assem_sequences = []  # 清空 assem_sequences 列表
                # 然后进行您的循环迭代和序列添加操作

                ss, ss_no_beam = compute_merged_sequence(filepath, beam_width, error_rate)
                ss_list = [ss[i:i + k] for i in range(len(ss) - k + 1)]
                ss_no_beam_list = [ss_no_beam[i:i + k] for i in range(len(ss_no_beam) - k + 1)]
                filepath_msa_sequence_list = [filepath_msa_sequence[i:i + k] for i in
                                              range(len(filepath_msa_sequence) - k + 1)]
                # print("filepath:",filepath)
                print("ss:", ss)
                print("nb:", ss_no_beam)
                print("fp:", filepath_msa_sequence)
                sequences = read_fasta(filepath)
                dbg_nodes, edges = build_dbg(sequences, k)
                G = nx.DiGraph()
                for kmer in dbg_nodes:
                    G.add_node(kmer)
                for (kmer1, kmer2), freq in edges.items():
                    G.add_edge(kmer1, kmer2, frequency=freq)
                    # print(f"{kmer1} -> {kmer2} : {freq}")
                BSDBG_assem_sequence = comparision_seq(ss_list, G, k)
                BSDBG_no_beam_assem_sequence = comparision_seq(ss_no_beam_list, G, k)
                fp_MSA_sequence = comparision_seq(filepath_msa_sequence_list, G, k)
                # print("组装后的序列：", assem_sequence)
                # yes_file_sequences为原始序列yes.fa


                # BSDBG_recovery_rate = data_recovery_rate(yes_file_sequences, assem_sequence)
                # msa_BSDBG_recovery_rate = data_recovery_rate(yes_file_sequences, msa_assem_sequence)
                BSDBG_recovery_rate = data_recovery_rate(yes_file_sequences_string, BSDBG_assem_sequence)
                BSDBG_no_beam_recovery_rate = data_recovery_rate(yes_file_sequences_string,
                                                                 BSDBG_no_beam_assem_sequence)
                fp_MSA_sequence_recovery_rate = data_recovery_rate(yes_file_sequences_string, fp_MSA_sequence)
                # print(yes_file_sequences_string)
                print("assem_sequence", BSDBG_assem_sequence)
                print("msa_a_sequence", BSDBG_no_beam_assem_sequence)
                print("fp_MSA_sequence", fp_MSA_sequence)

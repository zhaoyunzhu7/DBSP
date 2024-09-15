from sequence_beam_search import *
from lev_path import print_edit_operations
from sequence_alignment_diagram import *
import Levenshtein
from find_clusting_graph_max_node import *
from align_seq import align_seqs
from collections import defaultdict
# 基于编辑距离（Levenshtein距离）将相似的序列聚类成图，并找到每个图中连接边数最多的节点。
# 读取eDNAs_all.fa文件并将序列存储到一个列表中
# 读取FASTA文件中的序列
def write_before_align(path_comparision, path_fasta_file):
    with open(path_fasta_file, "w") as output_file:
        for i, sequence in enumerate(path_comparision, start=1):
            output_file.write(f">{i}\n{sequence}\n")

# 计算Levenshtein距离
def levenshtein_distance(seq1, seq2):
    return Levenshtein.distance(seq1, seq2)


def build_graph(sequences):  # 计算序列之间的最小编辑距离，构建图
    graph = nx.Graph()
    # 添加节点
    for i, sequence in enumerate(sequences):
        graph.add_node(i, sequence=sequence)
    # 添加边，连接每个序列与编辑距离最小的序列
    for i in range(len(sequences)):
        min_distance = float('inf')
        min_distance_sequence = None
        for j in range(len(sequences)):
            if i != j:
                dist = distance(sequences[i], sequences[j])
                if dist < min_distance:
                    min_distance = dist
                    min_distance_sequence = j
        if min_distance_sequence is not None:
            graph.add_edge(i, min_distance_sequence, weight=min_distance)
    return graph

def build_dbg(sequences, k):
    dbg_nodes = set()
    edges = defaultdict(int)
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            dbg_nodes.add(kmer)
            if i < len(seq) - k:
                next_kmer = seq[i + 1:i + k + 1]
                edges[(kmer, next_kmer)] += 1

    # 去除频率小于阈值的边
    filtered_edges = {(kmer1, kmer2): freq for (kmer1, kmer2), freq in edges.items() if freq >= 2}

    return dbg_nodes, filtered_edges

def find_max_degree_node_and_sequence(graph, sequences):
    degrees = dict(graph.degree())
    max_degree_node = max(degrees, key=degrees.get)
    max_degree_sequence = sequences[max_degree_node]
    return max_degree_node, max_degree_sequence

def read_fasta(filename):
    sequences = []
    with open(filename, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(str(record.seq))
    return sequences
def compute_merged_sequence(filepath, beam_width, error_rate):
    sequences = read_fasta(filepath)
    cluster_graph = build_graph(sequences)
    max_degree_node, max_degree_sequence = find_max_degree_node_and_sequence(cluster_graph, sequences)
    max_degree_node_neighbors = []#保存最大度节点以及相邻节点序列
    max_degree_node_neighbors.append(max_degree_sequence)
    sequence_graph = create_sequence_graph(max_degree_sequence)
    for neighbor in cluster_graph.neighbors(max_degree_node):
        neighbor_sequence = cluster_graph.nodes[neighbor]['sequence']
        # max_degree_node_neighbors.append(neighbor_sequence)
        re_pos, del_pos, ins_pos = print_edit_operations(max_degree_sequence, neighbor_sequence)
        print("ins_pos:",ins_pos)
        # 更新图的信息
        sequence_graph = update_sequence_graph(max_degree_sequence, sequence_graph, re_pos, ins_pos, del_pos)
        max_degree_node_neighbors.append(neighbor_sequence)
    for node, data in sequence_graph.nodes(data=True):
        print(f"Node: {node}, Data: {data}")

    # visualize_sequence_graph(sequence_graph)  # 可视化
    beam_search_sequences = beam_search(sequence_graph, beam_width)
    # write_before_align(before_align_sequences, "D:\\Python\\PythonProject\\clu_ass_graph\\result\\before_align_sequences.fa")
    # msa_correct_alignment_sequence = align_seqs("D:\\Python\\PythonProject\\clu_ass_graph\\result\\before_align_sequences.fa")#没有束搜索，熔融
    BSDBG_fasta = "D:\\Python\\PythonProject\\clu_ass_graph\\result\\BSDBG_sequences_{:.2f}.fa".format(error_rate)
    BSDBG_no_beam_fasta = "D:\\Python\\PythonProject\\clu_ass_graph\\result\\BSDBG_no_beam_sequences_{:.2f}.fa".format(error_rate)
    write_before_align(beam_search_sequences, BSDBG_fasta)
    write_before_align(max_degree_node_neighbors, BSDBG_no_beam_fasta)
    BSDBG_aln = "D:\\Python\\PythonProject\\clu_ass_graph\\result\\BSDBG_sequences_{:.2f}.aln".format(error_rate)
    BSDBG_no_beam_aln = "D:\\Python\\PythonProject\\clu_ass_graph\\result\\BSDBG_no_beam_sequences_{:.2f}.aln".format(error_rate)
    BSDBG_correct_alignment_sequence = align_seqs(BSDBG_fasta, BSDBG_aln)
    BSDBG_no_beam_correct_alignment_sequence = align_seqs(BSDBG_no_beam_fasta, BSDBG_no_beam_aln)
    return BSDBG_correct_alignment_sequence, BSDBG_no_beam_correct_alignment_sequence
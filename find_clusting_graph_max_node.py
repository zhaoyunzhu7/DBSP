"""
通过最小编辑距离构建聚类图，并找到每个簇中的最大度节点以及它的相邻节点

"""
from Levenshtein import distance
import networkx as nx
from Bio import SeqIO

# 读取1.fa文件中的序列
def read_sequences_from_file(file_path):
    sequences = []
    with open(file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(str(record.seq))
    return sequences


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

# 找到图中度最大的节点及其对应的序列
def find_max_degree_node_and_sequence(graph, sequences):
    degrees = dict(graph.degree())
    max_degree_node = max(degrees, key=degrees.get)
    max_degree_sequence = sequences[max_degree_node]
    return max_degree_node, max_degree_sequence

if __name__ == "__main__":
    # 读取1.fa文件中的序列
    file_path = "mix_eDNAs_all_iteration1.fa"
    sequences = read_sequences_from_file(file_path)  # list
    # print(sequences)
    # 计算序列之间的最小编辑距离，构建图
    sequence_graph = build_graph(sequences)  # 得到所有序列的全连接图
    print(sequence_graph)
    # 找到图中度最大的节点及其对应的序列
    max_degree_node, max_degree_sequence = find_max_degree_node_and_sequence(sequence_graph, sequences)
    # 输出结果
    print("具有最大度的节点为:", max_degree_node, "度为:", sequence_graph.degree(max_degree_node))
    print("具有最大度的的序列为:", max_degree_sequence)

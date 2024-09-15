"""
在序列比对图上进行波束搜索，得到k条高质量的序列

"""

from sequence_alignment_diagram import *


def beam_search(graph, k):
    # 初始化beam，每个元素是一个路径和对应的权值
    start_node = (-1, 'root')
    beam = [([start_node], 0)]
    new_beam = beam
    neighbor_dict = graph.adj
    beam_final = []

    while new_beam:
        # 进行波束搜索
        if k == 0:
            break
        # for _ in range(graph.number_of_nodes()):
        new_beam = []

        # 对于当前beam中的每条路径
        for path, score in beam:
            # 获取当前路径的最后一个节点
            current_node = path[-1]
            # if current_node == (-1, 'end'):
            #     break

            # 对于当前节点的每个邻居
            for neighbor in graph.successors(current_node):

                if neighbor == (-1, 'end'):
                    new_beam.append((path + [neighbor], 100))
                    break
                # 计算新路径的权值
                new_score = score + graph[current_node][neighbor].get('weight', 1)

                # 将新路径添加到新的beam中
                new_beam.append((path + [neighbor], new_score))

        # 根据权值对新的beam进行排序，保留前k个路径
        new_beam.sort(key=lambda x: x[1], reverse=True)
        if new_beam == []:
            break
        beam = new_beam[:k]
        # print("beam", beam)

        for b in beam:
            if b[0][-1] == (-1, 'end'):
                # print(b)
                beam_final.append(b)
                k -= 1

    sequences = []
    for path in beam_final:
        sequence = ''
        for seq in path[0]:
            if seq[0] != -1:
                sequence += seq[1]
        sequences.append(sequence)
        # print(sequences)
    return sequences


if __name__ == '__main__':
    seq = 'TCTGCAAGCGT'
    G = create_sequence_graph(seq)
    print(G)
    visualize_sequence_graph(G)

    re_pos = {0: 'C'}
    ins_pos = {6: 'G', 3: 'G'}
    del_pos = {8: 'C'}
    G_1 = update_sequence_graph(seq, G, re_pos, ins_pos, del_pos)
    visualize_sequence_graph(G_1)

    re_pos1 = {0: 'C'}
    ins_pos1 = {6: 'G', 3: 'G'}
    del_pos1 = {8: 'C'}
    G_2 = update_sequence_graph(seq, G, re_pos1, ins_pos1, del_pos1)
    visualize_sequence_graph(G_2)

    width = 16
    sequences = beam_search(G_2, width)
    print(sequences)

    # 连续错误目前还不行，后面我在改，你试试目前不加连续错误时候的效果吧。

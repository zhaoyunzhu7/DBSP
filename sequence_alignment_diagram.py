# """
# 通过聚类后的序列构建序列比对图
#
import networkx as nx
import matplotlib.pyplot as plt


def visualize_sequence_graph(G):  # 可视化
    plt.figure(figsize=(10, 10))
    pos = nx.circular_layout(G)

    # labels = nx.get_node_attributes(G, 'base')
    edge_labels = nx.get_edge_attributes(G, 'weight')

    nx.draw(G, pos, with_labels=True, node_size=700, node_color="skyblue", font_size=10, font_color="black",
            font_weight="bold", arrowsize=20)
    # nx.draw_networkx_labels(G, pos, labels=labels)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

    plt.show()


def create_sequence_graph(backbone):
    G = nx.DiGraph()
    G.add_node((-1, 'root'))

    for i, base in enumerate(backbone):
        if not G.has_node((i, base)):
            G.add_node((i, base))

        # 添加边
        if i < len(backbone) - 1:
            next_base = backbone[i + 1]
            if not G.has_node((i + 1, next_base)):
                G.add_node((i + 1, next_base))

            G.add_edge((i, base), (i + 1, next_base), weight=1)

    G.add_edge((-1, 'root'), (0, backbone[0]), weight=1)
    G.add_node((-1, 'end'))
    G.add_edge((len(backbone) - 1, backbone[-1]), (-1, 'end'), weight=1)
    return G


def update_sequence_graph(backbone, G, re_pos, ins_pos, del_pos):
    # global curr_node
    backbone_list = []
    for i, base in enumerate(backbone):
        backbone_list.append((i, base))

    # 处理发生错误的边
    for index, base in re_pos.items():
        if not G.has_node((index, base)):
            G.add_node((index, base))

            for node in backbone_list:
                if node[0] == int(index):
                    current_node = node

            pre_node = list(G.predecessors(current_node))[0]
            last_node = list(G.successors(current_node))[0]
            G.add_edge(pre_node, (index, base), weight=1)
            G.add_edge((index, base), last_node, weight=1)
        else:
            pre_node = list(G.predecessors((index, base)))[0]
            last_node = list(G.successors((index, base)))[0]
            G[pre_node][(index, base)]['weight'] += 1
            G[(index, base)][last_node]['weight'] += 1

    for index, base in ins_pos.items():
        for node in backbone_list:
            if node[0] == int(index):
                current_node = node
        pre_node = list(G.predecessors(current_node))[0]

        if not G.has_node((-index, base)):  # 为了在图中区分插入错误位置和主干DNA序列的位置，选择使用负数表示插入错误的位置。
            G.add_node((-index, base))  # 如果插入错误发生在DNA序列的位置 index 处，那么将这个位置表示为 -index。
            G.add_edge(pre_node, (-index, base), weight=1)
            G.add_edge((-index, base), current_node, weight=1)
        else:
            if not G.has_edge(pre_node, (-index, base)):
                G.add_edge(pre_node, (-index, base), weight=1)
            else:
                G[pre_node][(-index, base)]['weight'] += 1
            if not G.has_edge((-index, base), current_node):
                G.add_edge((-index, base), current_node, weight=1)
            else:
                G[(-index, base)][current_node]['weight'] += 1

    for index, base in del_pos.items():
        for node in backbone_list:
            if node[0] == int(index):
                current_node = node
        pre_node = list(G.predecessors(current_node))[0]
        last_node = list(G.successors(current_node))[0]

        if G.has_edge(pre_node, last_node):
            G[pre_node][last_node]['weight'] += 1
        else:
            G.add_edge(pre_node, last_node, weight=1)

    # 对没发生错误位置的边加权重
    list_1 = []
    list_1.append((-1, 'root'))
    list_1.extend(backbone_list)
    list_1.append((-1, 'end'))

    for i, node_i in enumerate(list_1):
        if i != len(list_1) - 1:
            if node_i[0] not in ins_pos and list_1[i + 1][0] not in ins_pos:
                if node_i[0] not in re_pos and list_1[i + 1][0] not in re_pos:
                    if node_i[0] not in del_pos and list_1[i + 1][0] not in del_pos:
                        G[node_i][list_1[i + 1]]['weight'] += 1

            elif node_i[0] in ins_pos and list_1[i + 1][0] not in ins_pos:
                G[node_i][list_1[i + 1]]['weight'] += 1

    return G


if __name__ == '__main__':
    seq = 'TCTGCAAGCGT'
    G = create_sequence_graph(seq)
    print(G)
    visualize_sequence_graph(G)

    re_pos = {0: 'C'}
    ins_pos = {6: 'G'}
    del_pos = {8: 'C'}
    G_1 = update_sequence_graph(seq, G, re_pos, ins_pos, del_pos)
    visualize_sequence_graph(G_1)

    re_pos1 = {0: 'A'}
    ins_pos1 = {6: 'T'}
    del_pos1 = {8: 'C'}
    G_2 = update_sequence_graph(seq, G_1, re_pos1, ins_pos1, del_pos1)
    visualize_sequence_graph(G_2)

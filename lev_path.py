from typing import List
import numpy as np
import copy


def levenshtein_distance(s1, s2):
    len_s1 = len(s1)
    len_s2 = len(s2)

    # 创建一个矩阵并初始化第一行和第一列（二维数组，行数和列数分别为len_s1+1和len_s2+1）
    matrix = np.zeros((len_s1 + 1, len_s2 + 1), dtype=int)
    for i in range(len_s1 + 1):
        matrix[i][0] = i
    for j in range(len_s2 + 1):
        matrix[0][j] = j

    # 填充矩阵
    for i in range(1, len_s1 + 1):
        for j in range(1, len_s2 + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            matrix[i][j] = min(matrix[i - 1][j] + 1, matrix[i][j - 1] + 1, matrix[i - 1][j - 1] + cost)

    return matrix


def get_edit_paths(matrix: List[List[int]], source: List[str], target: List[str]):
    # 用来保存第一条编辑路径
    first_path = []

    def _edit_path(i, j, optimal_path):
        nonlocal first_path  # 使得在函数内部可以修改外部变量

        if not first_path:  # 如果已经找到了第一条路径，就不再继续递归
            if i > 0 and j > 0:
                diagonal = matrix[i - 1][j - 1]
                vertical = matrix[i - 1][j]
                horizontal = matrix[i][j - 1]
                current = matrix[i][j]

                flag = False
                minimal = min(diagonal, min(vertical, horizontal))

                if diagonal == minimal:
                    new_i = i - 1
                    new_j = j - 1
                    path_ = copy.deepcopy(optimal_path)
                    if diagonal == current - 1:
                        path_.append(f"Replace | {new_i + 1} | {target[new_j]}")
                        _edit_path(new_i, new_j, path_)
                    elif source[new_i] == target[new_j]:
                        flag = True
                        path_.append(f"Keep | {new_i + 1}")
                        _edit_path(new_i, new_j, path_)

                if not flag:
                    if vertical == minimal:
                        new_i = i - 1
                        new_j = j
                        path_ = copy.deepcopy(optimal_path)
                        path_.append(f"Delete | {new_i + 1}")
                        _edit_path(new_i, new_j, path_)

                    if horizontal == minimal:
                        new_i = i
                        new_j = j - 1
                        path_ = copy.deepcopy(optimal_path)
                        path_.append(f"Insert | {new_j + 1} | {target[new_j]}")
                        _edit_path(new_i, new_j, path_)

            else:
                while i > 0:
                    i -= 1
                    optimal_path.append(f"Delete | {i + 1}")

                while j > 0:
                    j -= 1
                    optimal_path.append(f"Insert | {j + 1} | {target[j]}")

                first_path = list(reversed(optimal_path))

    # 获取行数和列数
    row_len = len(matrix) - 1
    col_len = len(matrix[0]) - 1
    # 调用
    _edit_path(row_len, col_len, optimal_path=[])
    return [first_path]


# s1 = "TCTGCAAGGCT"
# s2 = "CCTGGAGAGGT"
# matrix = levenshtein_distance(s1, s2)
# print(matrix)
# paths = get_edit_paths(matrix, source=list(s1), target=list(s2))
# for path in paths:
#      print(path)
def print_edit_operations(s1, s2):
    matrix = levenshtein_distance(s1, s2)
    paths = get_edit_paths(matrix, source=list(s1), target=list(s2))
    re_pos = {}
    del_pos = {}
    ins_pos = {}
    for path in paths:
        source_str = list(s1)
        # target_str = list(s2)

        for step in path:
            action, *rest = step.split(" | ")
            if action == "Replace":
                position = int(rest[0]) - 1
                new_char = rest[1]
                source_str[position] = new_char
                re = position  # （针对正确序列发生错误的位置）替换位置
                # print(f"Replace '{s1[position]}' at position {re} with '{new_char}': {''.join(source_str)}")
                re_pos[re] = new_char
            elif action == "Delete":
                position = int(rest[0]) - 1
                del_char = source_str[position]
                de = position  # 删除位置
                # del source_str[position]
                # print(f"Delete character at position {de}: {''.join(source_str)}")
                del_pos[de] = del_char
                # print(del_char)
            elif action == "Insert":
                position = int(rest[0]) - 1
                new_char = rest[1]
                source_str.insert(position, new_char)
                ins = position  # 插入位置
                # print(f"Insert '{new_char}' at position {ins}: {''.join(source_str)}")
                ins_pos[ins] = new_char

    return re_pos, del_pos, ins_pos


if __name__ == '__main__':
    s1 = "CCTGCAGAGTAGGCGGGGACCATACGTTTTCCAACTTGTTAAGTAGACGATCATAGGCTGCTATGGTGCACTATCGGCCCCTCTTCAAGTTGCCTCCAGTTCATGTTGTCACATGGTGGGACCAGTTTGCAGAGGTGATGGACGCGCACCG"
    s2 = "CCTGCAGAGTAGGCGGGGGACCATACGTTTTCCCAACTTGTTAAGTAGACGATCATAGGCCTGCTATGGGTGCACTATCGGCCCCTCTTCAAGTTGCCTCCAGTTCATGTTTGTCACATGGTTGGGACCAGTTTGCAGAGGGTGATGGGACGCGCACCG"

    re_pos, del_pos, ins_pos = print_edit_operations(s1, s2)
    print("re_pos:", re_pos)
    print("del_pos:", del_pos)
    print("ins_pos:", ins_pos)

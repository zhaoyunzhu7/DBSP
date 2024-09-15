from graph import levenshtein_distance


def data_recovery_rate(original_sequence, sequence):
    distance = levenshtein_distance(original_sequence, sequence)
    # print(original_sequence[0])
    # print(sequence)
    recovery_rate = (len(original_sequence) - distance) / len(original_sequence)

    return recovery_rate


if __name__ == "__main__":

    # 读取原始数据和恢复的DNA序列
    original_data_DBGPS = read_fasta("D:\\Python\\PythonProject\\clu_ass_graph\\yes_no_primer_index.fa")  # 用实际的原始数据替换
    original_data_BSDBG = read_fasta("D:\\Python\\PythonProject\\clu_ass_graph\\yes_no_primer_index.fa")
    DBGPS_recovered_sequences = read_fasta("D:\\Python\\PythonProject\\clu_ass_graph\\result\\assem_sequences_DBGPS.fa")
    DBGPS_recovery_rate = data_recovery_rate(original_data_DBGPS, DBGPS_recovered_sequences)
    BSDBG_recovered_sequences = read_fasta("D:\\Python\\PythonProject\\clu_ass_graph\\assem_sequences.fa")
    BSDBG_recovery_rate = data_recovery_rate(original_data_BSDBG, BSDBG_recovered_sequences)
    print(f"DBGPS数据恢复率：{DBGPS_recovery_rate:.2f}%")
    print(f"BSDBG数据恢复率：{BSDBG_recovery_rate:.2f}%")
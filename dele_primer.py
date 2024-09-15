def remove_primer_from_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequence = ""
        for line in infile:
            if line.startswith(">"):
                # 如果是FASTA文件的标题行，写入标题
                if sequence:
                    # 如果已经有序列数据，删除引物后写入文件
                    trimmed_sequence = sequence[32:-32]
                    outfile.write(f">{title}\n{trimmed_sequence}\n")
                title = line.strip()
                sequence = ""
            else:
                sequence += line.strip()

        # 处理最后一个序列
        if sequence:
            trimmed_sequence = sequence[32:-32]
            outfile.write(f">{title}\n{trimmed_sequence}\n")


input_file = "encode.fasta"
output_file = "data/encode_dele_primer.fasta"
remove_primer_from_fasta(input_file, output_file)

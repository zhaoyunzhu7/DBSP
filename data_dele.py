import os

# 指定目标目录
directory_path = "result"

# 获取目录下的所有文件
files = os.listdir(directory_path)

# 遍历并删除文件
for file in files:
    file_path = os.path.join(directory_path, file)
    if os.path.isfile(file_path):
        os.remove(file_path)
        print(f"已删除文件: {file_path}")

print(f"目录 '{directory_path}' 下的所有文件已删除")

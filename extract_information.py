import os

def extract_subfolder_names(folder_path):
    result = []
    # 遍历目标文件夹中的所有条目
    for entry in os.listdir(folder_path):
        # 获取条目的完整路径
        full_path = os.path.join(folder_path, entry)
        # 检查是否是目录
        if os.path.isdir(full_path):
            # 查找所有下划线的位置
            underscores = [i for i, char in enumerate(entry) if char == '_']
            # 确保至少有2个下划线
            if len(underscores) >= 2:
                # 获取第二个下划线的位置
                second_underscore_pos = underscores[1]
                # 截取并保存结果
                result.append(entry[:second_underscore_pos])
    return result

# 使用示例
if '__name__' == '__main__':
    folder_path = "/root/autodl-tmp/crossdocked_v1.1_rmsd1.0"  # 替换为你的文件夹路径
    result_list = extract_subfolder_names(folder_path)



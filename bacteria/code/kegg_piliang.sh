#!/bin/bash

# 设置输入和输出目录
input_dir="/home/wtao/projects/czh/P24052701/bacteria/outdir/Prokka"
output_dir="/home/wtao/projects/bacteria/kegg_res"  # 替换为你的输出目录

# 遍历目录中的每个 .faa 文件
for folder in "$input_dir"/*; do
    if [ -d "$folder" ]; then
        folder_name=$(basename "$folder")
        faa_file="$folder/$folder_name.faa"
        
        # 构造输出文件名
        output_file="$output_dir/${folder_name}_kegg.txt"
        
        # 执行命令
        exec_annotation -o "$output_file" --cpu 20 --format mapper -E 1e-5 "$faa_file"
    fi
done

echo "批量处理完成！"
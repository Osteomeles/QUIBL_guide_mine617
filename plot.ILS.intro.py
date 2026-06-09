/data/01/user187/0.data_all/01.Talpidae/11.diff_tree/01.3.4d_windows/04.quibl/summary_with_only_sygansu/summaryFileanalysis/run1/each_species_pair_use_totalIntroProp_totalILSProp/plot_ils_intro/plot.ILS.intro.py
#!/usr/bin/env python3
"""
绘制基于自定义物种顺序的渐渗支持比例热图（下三角）和ILS支持比例热图（上三角）
输入：
  1. 渐渗数据文件：两列，制表符分隔，第一列为"物种A_物种B"，第二列为支持比例（0-1）
  2. ILS数据文件：两列，制表符分隔，第一列为"物种A_物种B"，第二列为支持比例（0-1）
  3. 物种顺序文件：每行一个物种名，按行顺序定义行列顺序
输出：
  完整热图（下三角渐渗，上三角ILS），颜色从白色(0)到莫兰迪哑光紫(1)
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches

# ============================================================
# 用法示例：
#   python plot_intro_ils_heatmap.py intro_data.txt ils_data.txt species_order.txt output_prefix
#
# 渐渗数据文件格式示例（intro_data.txt）：
#   speciesA_speciesB    0.75
#   speciesA_speciesC    0.32
#   speciesB_speciesC    0.98
#
# ILS数据文件格式示例（ils_data.txt）：
#   speciesA_speciesB    0.25
#   speciesA_speciesC    0.68
#   speciesB_speciesC    0.02
#
# 物种顺序文件示例（species_order.txt）：
#   speciesA
#   speciesB
#   speciesC
# ============================================================

def parse_args():
    if len(sys.argv) < 5:
        sys.exit("""
用法: python plot_intro_ils_heatmap.py <intro_data_file> <ils_data_file> <species_order_file> <output_prefix>

  intro_data_file   : 渐渗数据文件，两列制表符分隔，第一列为"物种A_物种B"，第二列为支持比例(0-1)
  ils_data_file     : ILS数据文件，两列制表符分隔，第一列为"物种A_物种B"，第二列为支持比例(0-1)
  species_order_file: 每行一个物种名，定义热图行列顺序（从上到下，从左到右）
  output_prefix     : 输出文件前缀（将生成 .pdf 和 .png 文件）
        """)
    return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

def read_species_order(order_file):
    """读取物种顺序文件，返回物种列表（保持顺序）"""
    species_list = []
    with open(order_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:  # 跳过空行
                species_list.append(line)
    if not species_list:
        sys.exit(f"错误: 物种顺序文件 {order_file} 为空或格式不正确")
    print(f"读取到 {len(species_list)} 个物种，顺序如下:")
    for i, sp in enumerate(species_list):
        print(f"  {i+1}. {sp}")
    return species_list

def read_data(data_file, species_list, data_type=""):
    """
    读取数据文件，构建矩阵
    返回: n x n 的numpy矩阵（填充数值，其余为NaN）
    data_type: 用于提示信息
    """
    n = len(species_list)
    # 创建物种名到索引的映射
    species_to_idx = {sp: i for i, sp in enumerate(species_list)}
    
    # 初始化矩阵为NaN
    mat = np.full((n, n), np.nan)
    
    # 记录读取了多少条数据
    data_count = 0
    missing_species = set()
    
    with open(data_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 2:
                print(f"警告 ({data_type}): 第{line_num}行格式不正确（需要至少2列），跳过: {line}")
                continue
            
            species_pair = parts[0].strip()
            try:
                value = float(parts[1].strip())
                # 确保值在0-1范围内
                if value < 0 or value > 1:
                    print(f"警告 ({data_type}): 第{line_num}行数值超出[0,1]范围，跳过: {line}")
                    continue
            except ValueError:
                print(f"警告 ({data_type}): 第{line_num}行数值转换失败，跳过: {line}")
                continue
            
            # 解析物种对
            if '_' not in species_pair:
                print(f"警告 ({data_type}): 第{line_num}行物种对格式不正确（需要下划线连接），跳过: {line}")
                continue
            
            sp1, sp2 = species_pair.split('_', 1)  # 只分割第一个下划线
            
            # 检查物种是否在顺序列表中
            if sp1 not in species_to_idx:
                missing_species.add(sp1)
                continue
            if sp2 not in species_to_idx:
                missing_species.add(sp2)
                continue
            
            # 获取索引
            idx1 = species_to_idx[sp1]
            idx2 = species_to_idx[sp2]
            
            # 填充矩阵（对称填充，但后续下三角和上三角会分开使用）
            mat[idx1, idx2] = value
            mat[idx2, idx1] = value
            data_count += 1
    
    if missing_species:
        print(f"警告 ({data_type}): 以下物种在数据文件中出现但不在物种顺序列表中: {missing_species}")
    
    if data_count == 0:
        sys.exit(f"错误 ({data_type}): 未从数据文件 {data_file} 中读取到任何有效数据")
    
    print(f"成功读取 {data_count} 条数据记录 ({data_type})")
    
    return mat

def plot_heatmap(intro_mat, ils_mat, species_list, output_prefix, max_val=1.0, min_val=0.0):
    """
    绘制组合热图：下三角渐渗，上三角ILS
    intro_mat: n x n 矩阵，渐渗数据
    ils_mat: n x n 矩阵，ILS数据
    species_list: 物种列表（顺序）
    output_prefix: 输出文件前缀
    max_val: 颜色映射最大值（默认1）
    min_val: 颜色映射最小值（默认0）
    """
    n = len(species_list)
    
    # 创建组合矩阵：下三角用intro_mat，上三角用ils_mat，对角线设为NaN
    combined_mat = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if i > j:  # 下三角，渐渗
                combined_mat[i, j] = intro_mat[i, j]
            elif i < j:  # 上三角，ILS
                combined_mat[i, j] = ils_mat[i, j]
            # 对角线保持NaN
    
    # 创建莫兰迪风格的哑光紫色渐变（基于R脚本的色系：白色→灰色→purple4）
    # 使用更柔和的紫色调
    colors = ["#FFFFFF", "#E8E0F0", "#D1C4E0", "#B09EC8", "#8F78B0", "#6E5298", "#4D2C80", "#3A1F66"]
    # 进一步柔化，降低饱和度，模拟莫兰迪色系
    cmap = LinearSegmentedColormap.from_list("morandi_purple", colors, N=256)
    
    # 掩盖NaN值（对角线）
    masked = np.ma.masked_invalid(combined_mat)
    
    # 计算图形大小（根据物种数量自适应）
    fig_w = max(8, n * 0.6)
    fig_h = max(7, n * 0.55)
    
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    
    # 绘制热图
    im = ax.imshow(
        masked,
        cmap=cmap,
        vmin=min_val,
        vmax=max_val,
        interpolation="none",
        aspect="equal"
    )
    
    # 在每个非NaN的格子中写入数值（三位小数）
    for i in range(n):
        for j in range(n):
            val = combined_mat[i, j]
            if np.isfinite(val):
                # 格式化数值：保留三位小数，0显示为"0"
                if val == 0:
                    text = "0"
                else:
                    text = f"{val:.3f}"
                
                # 根据背景颜色深度决定文字颜色
                if val > 0.6:  # 阈值可调整
                    text_color = "white"
                else:
                    text_color = "black"
                
                ax.text(
                    j, i, text,
                    ha="center",
                    va="center",
                    fontsize=7 if n > 20 else 8,
                    color=text_color
                )
                
                # 对于非0的数值框，添加红色边框（加粗）
                if val != 0:
                    rect = plt.Rectangle((j-0.5, i-0.5), 1, 1, 
                                        fill=False, edgecolor="red", linewidth=1.5)
                    ax.add_patch(rect)
    
    # 设置坐标轴
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(species_list, rotation=90, fontsize=9)
    ax.set_yticklabels(species_list, fontsize=9)
    
    # 设置坐标轴范围
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    
    # 设置标题和标签
    ax.set_title("Genomic support patterns", fontsize=13)
    ax.set_xlabel("Species")
    ax.set_ylabel("Species")
    
    # 隐藏边框
    for sp in ax.spines.values():
        sp.set_visible(False)
    
    # 添加黑色边框（粗线）分隔每个格子
    # 绘制所有网格线（包括对角线位置）
    # 垂直网格线
    for x in range(n + 1):
        ax.axvline(x - 0.5, color="black", linewidth=1.2, linestyle="-")
    # 水平网格线
    for y in range(n + 1):
        ax.axhline(y - 0.5, color="black", linewidth=1.2, linestyle="-")
    
    ax.tick_params(which="minor", bottom=False, left=False)
    
    # 添加颜色条
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Support ratio (proportion of windows)", fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    
    # 在图例位置添加标注说明下三角和上三角的含义
    # 在左下角和右上角添加小标注
    ax.text(0.05, 0.05, "Introgression", transform=ax.transAxes, 
            fontsize=9, color="#4D2C80", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))
    ax.text(0.95, 0.95, "ILS", transform=ax.transAxes, 
            fontsize=9, color="#4D2C80", fontweight="bold",
            ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))
    
    plt.tight_layout()
    
    # 保存图片
    out_pdf = f"{output_prefix}_heatmap.pdf"
    out_png = f"{output_prefix}_heatmap.png"
    plt.savefig(out_pdf, format="pdf", dpi=300, bbox_inches="tight")
    plt.savefig(out_png, format="png", dpi=300, bbox_inches="tight")
    
    print(f"PDF 已保存至: {out_pdf}")
    print(f"PNG 已保存至: {out_png}")
    
    return fig, ax

def main():
    intro_file, ils_file, order_file, output_prefix = parse_args()
    
    # 1. 读取物种顺序
    species_list = read_species_order(order_file)
    
    # 2. 读取渐渗数据
    print("\n--- 读取渐渗数据 ---")
    intro_mat = read_data(intro_file, species_list, "渐渗")
    
    # 3. 读取ILS数据
    print("\n--- 读取ILS数据 ---")
    ils_mat = read_data(ils_file, species_list, "ILS")
    
    # 4. 绘制热图
    plot_heatmap(intro_mat, ils_mat, species_list, output_prefix, max_val=1.0, min_val=0.0)
    
    # 5. 输出矩阵文件（便于检查）
    # 渐渗矩阵文件
    intro_matrix_out = f"{output_prefix}_intro_matrix.tsv"
    with open(intro_matrix_out, 'w') as out:
        out.write("taxon\t" + "\t".join(species_list) + "\n")
        for i, taxon in enumerate(species_list):
            row_vals = []
            for j in range(len(species_list)):
                if np.isfinite(intro_mat[i, j]):
                    row_vals.append(f"{intro_mat[i, j]:.6f}")
                else:
                    row_vals.append("NA")
            out.write(taxon + "\t" + "\t".join(row_vals) + "\n")
    print(f"渐渗矩阵文件已保存至: {intro_matrix_out}")
    
    # ILS矩阵文件
    ils_matrix_out = f"{output_prefix}_ils_matrix.tsv"
    with open(ils_matrix_out, 'w') as out:
        out.write("taxon\t" + "\t".join(species_list) + "\n")
        for i, taxon in enumerate(species_list):
            row_vals = []
            for j in range(len(species_list)):
                if np.isfinite(ils_mat[i, j]):
                    row_vals.append(f"{ils_mat[i, j]:.6f}")
                else:
                    row_vals.append("NA")
            out.write(taxon + "\t" + "\t".join(row_vals) + "\n")
    print(f"ILS矩阵文件已保存至: {ils_matrix_out}")
    
    print("\n=== 运行完成 ===")

if __name__ == "__main__":
    main()

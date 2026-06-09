#!/usr/bin/env python3
"""
判断三联体的外群是否符合物种树拓扑
用法: python check_outgroup.py --input triplets.txt --tree species_tree.nwk --output result.txt
"""

import argparse
from ete3 import Tree
import sys

def read_triplets_file(filename):
    """读取三联体文件，每行格式: triplet: A_B_C outgroup: X"""
    triplets = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # 解析 "triplet: A_B_C outgroup: X" 格式
            parts = line.split()
            if len(parts) >= 4:
                triplet_str = parts[1]  # A_B_C
                outgroup = parts[3]     # X
                species = triplet_str.split('_')
                if len(species) == 3:
                    triplets.append({
                        'species': species,
                        'outgroup': outgroup,
                        'line': line
                    })
    return triplets

def get_lca_distance(tree, species1, species2):
    """获取两个物种在树上的距离（枝长之和）"""
    node1 = tree & species1
    node2 = tree & species2
    return tree.get_distance(species1, species2)

def check_outgroup_correct(tree, species, candidate_outgroup):
    """
    判断 candidate_outgroup 是否是三个物种中离另外两个最远的
    返回: (is_correct, reason, distances_dict)
    """
    a, b, c = species
    out = candidate_outgroup
    
    # 获取所有成对距离
    dist_ab = get_lca_distance(tree, a, b)
    dist_ac = get_lca_distance(tree, a, c)
    dist_bc = get_lca_distance(tree, b, c)
    
    # 找出候选外群
    # 外群应该满足：外群到另外两个物种的距离之和 > 另外两个物种之间的距离
    if out == a:
        dist_out_to_others = dist_ab + dist_ac
        dist_between_others = dist_bc
        is_correct = dist_out_to_others > dist_between_others
        reason = f"{out}到另外两个的距离和({dist_out_to_others:.6f}) vs 另外两个之间距离({dist_between_others:.6f})"
    elif out == b:
        dist_out_to_others = dist_ab + dist_bc
        dist_between_others = dist_ac
        is_correct = dist_out_to_others > dist_between_others
        reason = f"{out}到另外两个的距离和({dist_out_to_others:.6f}) vs 另外两个之间距离({dist_between_others:.6f})"
    elif out == c:
        dist_out_to_others = dist_ac + dist_bc
        dist_between_others = dist_ab
        is_correct = dist_out_to_others > dist_between_others
        reason = f"{out}到另外两个的距离和({dist_out_to_others:.6f}) vs 另外两个之间距离({dist_between_others:.6f})"
    else:
        is_correct = False
        reason = f"外群 {out} 不在三个物种 {a},{b},{c} 中"
    
    distances = {
        f"{a}-{b}": dist_ab,
        f"{a}-{c}": dist_ac,
        f"{b}-{c}": dist_bc
    }
    
    return is_correct, reason, distances

def find_true_outgroup(tree, species):
    """
    根据物种树找出三个物种中真正的离最远的外群
    返回真正的outgroup名称
    """
    a, b, c = species
    dist_ab = get_lca_distance(tree, a, b)
    dist_ac = get_lca_distance(tree, a, c)
    dist_bc = get_lca_distance(tree, b, c)
    
    # 外群应该是使得"外群到另外两个的距离和"最大的那个物种
    sum_a = dist_ab + dist_ac
    sum_b = dist_ab + dist_bc
    sum_c = dist_ac + dist_bc
    
    max_sum = max(sum_a, sum_b, sum_c)
    if max_sum == sum_a:
        return a
    elif max_sum == sum_b:
        return b
    else:
        return c

def main():
    parser = argparse.ArgumentParser(description='判断三联体外群是否符合物种树')
    parser.add_argument('-i', '--input', required=True, help='输入的三联体文件')
    parser.add_argument('-t', '--tree', required=True, help='物种树文件 (newick格式)')
    parser.add_argument('-o', '--output', default='outgroup_check_result.txt', help='输出文件')
    parser.add_argument('-v', '--verbose', action='store_true', help='显示详细信息')
    
    args = parser.parse_args()
    
    # 读取物种树
    try:
        tree = Tree(args.tree, format=1)
        print(f"成功读取物种树: {args.tree}")
    except Exception as e:
        print(f"读取物种树失败: {e}")
        sys.exit(1)
    
    # 读取三联体文件
    triplets = read_triplets_file(args.input)
    print(f"读取到 {len(triplets)} 个三联体")
    
    # 检查每个三联体
    results = []
    correct_count = 0
    
    for i, trip in enumerate(triplets, 1):
        species = trip['species']
        outgroup = trip['outgroup']
        
        # 检查物种是否都在树中
        missing = []
        for s in species:
            try:
                tree & s
            except:
                missing.append(s)
        
        if missing:
            is_correct = False
            true_outgroup = "N/A"
            reason = f"物种 {missing} 不在物种树中"
        else:
            is_correct, reason, distances = check_outgroup_correct(tree, species, outgroup)
            true_outgroup = find_true_outgroup(tree, species)
        
        results.append({
            'triplet': '_'.join(species),
            'species': species,
            'given_outgroup': outgroup,
            'true_outgroup': true_outgroup,
            'is_correct': is_correct,
            'reason': reason
        })
        
        if is_correct:
            correct_count += 1
    
    # 输出结果到控制台
    print("\n" + "="*80)
    print("检查结果")
    print("="*80)
    
    for r in results:
        status = "✓ 正确" if r['is_correct'] else "✗ 错误"
        print(f"{status} | 三联体: {r['triplet']} | 给定外群: {r['given_outgroup']} | 正确外群应为: {r['true_outgroup']}")
        if args.verbose and not r['is_correct']:
            print(f"    原因: {r['reason']}")
    
    print("\n" + "="*80)
    print(f"总结: {correct_count}/{len(results)} 个三联体的外群正确 ({correct_count/len(results)*100:.1f}%)")
    print("="*80)
    
    # 写入输出文件
    with open(args.output, 'w') as f:
        f.write("triplet\tgiven_outgroup\ttrue_outgroup\tis_correct\treason\n")
        for r in results:
            f.write(f"{r['triplet']}\t{r['given_outgroup']}\t{r['true_outgroup']}\t{r['is_correct']}\t{r['reason']}\n")
    
    print(f"\n详细结果已保存到: {args.output}")
    
    # 输出所有错误的行（方便检查）
    errors = [r for r in results if not r['is_correct']]
    if errors:
        print("\n" + "-"*40)
        print("外群错误的三联体（原格式，可用于后续处理）:")
        print("-"*40)
        for r in errors:
            # 输出与原文件相同的格式
            print(f"triplet: {r['triplet']} outgroup: {r['given_outgroup']}")
    
    return results

if __name__ == "__main__":
    main()

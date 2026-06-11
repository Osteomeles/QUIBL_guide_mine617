# QUIBL-使用QuIBL过程中的经验记录
自用经验，可能有不准确的地方

# 参考REFERENCE：

- Genomic architecture and introgression shape a butterfly radiation

  - 下载文章的附件里面有详细的QUIBL解读

  - QUIBL的github（https://github.com/miriammiyagi/QuIBL）里面还有很多后续分析的脚本可以使用（下面这里就用到了）

- Extensive Genome-Wide Phylogenetic Discordance Is Due to Incomplete Lineage Sorting and Not Ongoing Introgression in a Rapidly Radiated Bryophyte Genus

  - 下载文章的附件里面有QUIBL结果附表，可以看到作者是怎么处理QUIBL结果文件的

- 其他有QUIBL分析的文章，可以参考，看看其他人是如何构造输入文件和分析结果的

  - Incomplete lineage sorting and phenotypic evolution in marsupials
  - Genomic insights into the secondary aquatic transition of penguins

# 输入文件INPUT FILE：

- 将系统发育树中的物种分为四个物种一组的若干组合，其中一个是固定的外群物种，另外三个物种在剩余物种中随机组合得到，这里生成若干个四物种组合，每个四物种组合都单独运行一次QUIBL
- ```
  python3 01_comb_4spec.py
  ```
  这个脚本把所有需要随机组合的物种名称和固定外群写在脚本里面，就可以自动生成所有的组合写入out_four_species_array.txt文件
  这个脚本来自：https://wu-tz.github.io/2024/04/09/QuIBL%E6%96%B9%E6%B3%95%E6%A3%80%E9%AA%8C%E5%9F%BA%E5%9B%A0%E6%B8%90%E6%B8%97%E5%92%8C%E4%B8%8D%E5%AE%8C%E5%85%A8%E8%B0%B1%E7%B3%BB%E5%88%86%E9%80%89/
  里面还有一些关于QUIBL运行的解读或者脚本，可以进一步阅读了解
  这里我每个四物种组合都有一个单独的文件夹，放各自的输入文件，配置文件和输出文件
```
  while read -r a b c d
    do
    CHR="${a}"_"${b}"_"${c}"
    mkdir $CHR
  done < out_four_species_array.txt
```
## 树文件：

我这里用的是划窗口（5kb大小的窗口，窗口大小参考的是Genomic architecture and introgression shape a butterfly radiation，并且窗口再小的话我的服务器有点无法承受）得到的若干窗口树，每行一个nwk格式的树文件

由于我的窗口树太多了，我随机抽样了5000棵树作为输入（这里可以参考不同的文章抽样的大小，树越多，运行的时间成本越大；这里我后续做了多次重复来避免结果的随机性
 
```
    shuf -n 5000 window.tree > 5000.tree
```
随后我需要把1.sub.tree针对不同的四物种组合提取各自的子树作为输入
```
    while read -r a b c d
    do
    CHR="${a}_${b}_${c}"
    while read i
    do
    echo "$i" | /data/00/user/user187/miniconda3/bin/nw_prune - -v $a $b $c Outgroup >> $CHR/1.sub.tree
    done < 5000.tree
    done < out_four_species_array.txt
```
这里Outgroup是我外群的名字

## 配置文件sampleInputFile.txt
模板：
- /path_to_file/推荐使用绝对路径
- treefile: /path_to_file/1.sub.tree指定输入的树文件也就是上面的1.sub.tree
- OutputPath: /path_to_file/Output.1.csv指定输出文件路径，这里我推荐每个四物种组合都有一个单独的文件夹，放各自的输入文件，配置文件和输出文件，故我使用的输出文件名称都是一致的
- totaloutgroup: Outgroup指定外群物种名称
```
[Input]
treefile: /path_to_file/1.sub.tree
numdistributions: 2
likelihoodthresh: 0.01
numsteps: 50
gradascentscalar: 0.5
totaloutgroup: Outgroup
multiproc: True
maxcores:10
[Output]
OutputPath: /path_to_file/Output.1.csv
```
写一个循环给每个四物种组合生成配置文件放入各自的文件夹

# 运行
- 需要按照QUIBL的github里面的教程构建环境，这里我的环境在lzu
```
conda activate /data/00/user/user187/miniconda3/envs/QuIBL_env
```
把所有组合的运行命令都写入一个sh文件，提交任务，用pNormal就好，c用1，m用1，运行需要的时间成本是比较高的；命令行示例：
```
python /data/00/user/user187/00.apps/QuIBL/QuIBL.py ./path_to_file/sampleInputFile.txt
```
# 分析结果文件&可视化
首先，下面有一些脚本需要判断拓扑是否和物种树一致，物种树定义为count数最多的拓扑；为避免count数最多的拓扑实际并不是物种树拓扑的情况出现，这里首先进行一个检查，是否count最多的就是物种树一致拓扑：
- 列出count最多的拓扑
```
find ./ -name "*.Output.1.csv" | while read i                                                                 
do
awk -F',' 'NR==1 {next} {if($NF>max) {max=$NF; line=$0}} END {split(line, a, ","); print "triplet:", a[1], "outgroup:", a[2]}' $i
done >> check.list
```
- 看count最多拓扑是否就是物种树一致拓扑，这里需要输入一个nwk格式的物种树（我用Astral从窗口树里面推断的物种树）；输出结果很直白可以解答我们的疑惑
```
python3 check_outgroup.py -i check.list -t species_tree.nwk -o result.txt
```

先把所有的物种组合的结果都整合为一个：

每个物种组合的分析都会得到一个结果文件Output.csv

共有这些列：triplet,outgroup,C1,C2,mixprop1,mixprop2,lambda2Dist,lambda1Dist,BIC2Dist,BIC1Dist,count
- 具体意义参考：https://github.com/miriammiyagi/QuIBL

把所有组合的结果文件只保留首个文件的列名cat到一起得到一个总结文件output.1.all.txt（这里我把所有结果文件都软连接到了一个allfile文件夹下面，因为名称是相同的故我软连接的时候将其改名在其前面加了物种组合的前缀）
```
(head -n 1 /path_to_file/allfile/run1/spe1_spe2_spe3.Output.1.csv && for f in /path_to_file/allfile/run1/*.Output.1.csv; do tail -n +2 "$f"; done) > output.1.all.txt
```
不知道为什么我得到的这个文件有时候会有格式问题，我会：
```
dos2unix output.1.all.txt
```

①首先用QUIBL提供的一个脚本可以得到“有百分之多少的不一致位点（或全基因组所有位点）支持渐渗”
- 去QUIBL的github里面下载相应脚本summaryFileAnalysis.R
- 可以得到：显著支持渐渗的tree的数量占所有窗口树数量的比例；显著支持渐渗的tree的数量占与物种树不一致的窗口树的数量的比例
```
Rscript summaryFileAnalysis.R -i output.1.all.txt -o output > result.txt
```

②此外这里试图得到的是Extensive Genome-Wide Phylogenetic Discordance Is Due to Incomplete Lineage Sorting and Not Ongoing Introgression in a Rapidly Radiated Bryophyte Genus的Fig.5A

得到output.1.all.txt后计算三个数值，diffBIC,totalILSProp,totalIntroProp，把他们加载txt文件后面得到一个新文件
- diffBIC=BIC2-BIC1，当diffBIC大于10，说明ILSonly模型显著更优；当diffBIC小于-10，说明有渐渗的模型显著优于ILSonly模型
- totalILSProp = (mixprop1*count)/num_alltree
  - 这里num_alltree是所有的窗口树，因为我一共抽样5000棵树，故我直接用的5000
- totalIntroProp = (mixprop2*count)/num_alltree
- 得到output.1.all.plus.txt得把新得到的列列名0改为对应的diffBIC/totalILSProp/totalIntroProp
```
awk -F "," '{print $1","$2","$3","$4","$5","$6","$7","$8","$9","$10","$11","$9-$10","$5*$11/5000","$6*$11/5000}' output.1.all.txt >> output.1.all.plus.txt
```

去得到所有的物种对，随后手动vim把姐妹枝的物种对去掉
```
python3 get_species_pairs.py output.1.all.plus.txt species.pairs.txt
```
去掉与物种树一致的拓扑行
```
while read -r a b                                                                                                            
do
CHR="${a}"_"${b}"  
mkdir $CHR                                                                           
done < species.pairs.txt  

head -n 1 output.1.all.plus.txt > head

while read -r a b c d                          
do
awk -F "," -v val="$b" -v val2="$d" '$1 == val && $2 != val2' output.1.all.plus.txt >> output.1.all.plus.dis.txt.tmp
done < ../check.list

cat head output.1.all.plus.dis.txt.tmp > output.1.all.plus.dis.txt

while read -r a b                                                                                                            
do
CHR="${a}"_"${b}"  
python3 filter_triplet.py output.1.all.plus.dis.txt $a $b ./$CHR/output.1.all.plus.$a.$b.txt
done < species.pairs.txt
```
把第十二列diffBIC小于-10也就是显著支持渐渗的拓扑行筛选出来，然后把这些行的第十四列totalIntroProp求平均，得到每个物种对的渐渗位点平均总比例
```
while read -r a b
do
    CHR="${a}_${b}"
    cd $CHR
    avg=$(awk -F "," '{if($12 < -10) print $14}' output.1.all.plus.$a.$b.txt | grep -v "totalIntroProp" | awk '{sum+=$1; count++} END {if(count>0) print sum/count; else print "NA"}')
    cd ..
    echo "$CHR"$'\t'"$avg"
done < species.pairs.txt >> all.species_pairs.totalIntroProp.txt
```
把第十二列diffBIC大于10也就是显著支持ILSonly的拓扑行筛选出来，然后把这些行的第十三列totalILSProp求平均，得到每个物种对的ILS-only位点平均总比例
```
while read -r a b                              
do
    CHR="${a}_${b}"
    cd $CHR
    avg=$(awk -F "," '{if($12 > 10) print $13}' output.1.all.plus.$a.$b.txt | grep -v "totalILSProp" | awk '{sum+=$1; count++} END {if(count>0) print sum/count; else print "NA"}')   
    cd ..
    echo "$CHR"$'\t'"$avg"
done < species.pairs.txt >> all.species_pairs.totalILSProp.txt  
```
把NA变为0
```
sed -i 's|NA|0|g' all.species_pairs.totalILSProp.txt all.species_pairs.totalIntroProp.txt
```
可视化
```
mkdir plot_ils_intro
cd plot_ils_intro
python3 plot.ILS.intro.py all.species_pairs.totalIntroProp.txt all.species_pairs.totalILSProp.txt species.list ILS.intro.plot
```
- 这里species.list是一个一列的文件，每行一个物种名，是为了固定绘图时物种排列顺序
- 这里和Fig.5A不同的是，上三角和下三角都是QUIBL的输出

③这里做一个类似的图，只是totalILSProp、totalIntroProp两个比值的分母变成了全部不一致树，而不是所有的树
后面用到的species.pairs.txt和②里面的生成方式一样可以直接用output.1.all.txt去生成
首先去得到去除一致树的结果文件
```
while read -r a b c d                                                                                                        
do
awk -F "," -v val="$b" -v val2="$d" '$1 == val && $2 != val2' output.1.all.txt >> output.1.all.dis.txt.tmp
done < ../check.list

head -n 1 output.1.all.txt > head

cat head output.1.all.dis.txt.tmp > output.1.all.dis.txt
```
随后去给每行拓扑加三列同样是diffBIC,totalILSProp,totalIntroProp
这里diffBIC保持不变，totalILSProp,totalIntroProp的分母变成all_discordant_tree，也就是两个不一致树拓扑行的count的加和
```
python3 add_quibl_col.py output.1.all.dis.txt output.1.all.plus.dis.txt
```
接下来的做法就和之前一样了，首先把相同物种对的分到各自的文件夹下面
```
while read -r a b                                                                                                            
do
CHR="${a}"_"${b}"                                                                                 
mkdir $CHR                                                                                 
done < species.pairs.txt

while read -r a b                                                                                                            
do
CHR="${a}"_"${b}"                                                                                 
python3 filter_triplet.py output.1.all.plus.dis.txt $a $b ./$CHR/output.1.all.plus.$a.$b.txt
done < species.pairs.txt
```
后续直接复刻“把第十二列diffBIC小于-10也就是显著支持渐渗的拓扑行筛选出来，然后把这些行的第十四列totalIntroProp求平均，得到每个物种对的渐渗位点平均总比例”后的内容即可


# 注
这里作图的逻辑是根据Extensive Genome-Wide Phylogenetic Discordance Is Due to Incomplete Lineage Sorting and Not Ongoing Introgression in a Rapidly Radiated Bryophyte Genus的附表推理出来的，可能有不准确的地方

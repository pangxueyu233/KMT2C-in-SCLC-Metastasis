# COVID19

我负责淋巴细胞，雪兰做髓系细胞，赵蕾做上皮、内皮和间充质

我、雪兰、赵蕾都需要完成的部分：

1. 分群
   1. population percentage between ctrl, COVID19 and treat
   2. key marker of population
2. COVID VS CTRL 
   1. DGEs, GO/KEGG
   2. pesudotime of population
   3. dynamics gene --> GO/KEGG
3. COVID VS TREAT
   1. DGEs, GO/KEGG
   2. pesudotime of population
   3. dynamics gene --> GO/KEGG
4. CTRL  VS TREAT
   1. DGEs, GO/KEGG
   2. pesudotime of population
   3. dynamics gene --> GO/KEGG

5. 病人的phenotype的结合
   1. COVID19中高表达的细胞群基因或者细胞群是否病人也高或者含有
   2. TREAT之后高表达的细胞群基因或者细胞群是否下调，病人中，低表达这些基因和或者含有更少的这些细胞群的，是不是相对良好的性状。
   3. TREAT之后和CTRL比差异的细胞群或者基因是否扔呈现出易发病感染的状态，是否和初期患者相似，或者治疗后患者相似。

针对COVID vs Ctr，描述模型攻毒时，免疫系统应答，免疫激活，可能已经失效，但是需要证明活化通路有所活化：

CD8 T 细胞有激活和扩增，

barplot统计 ACE2和TMPRSS2/9在各个细胞群和各组间的差异

分析模板：

1. Harmony消除批次
2. 分群图(剔除Unknown)TSNE/openTSNE
3. 密度图展示三阶段比例+percentage统计图(barplot带P-value和所有的个体统计分布图)
4. 分群所用marker，DOTPLOT+投影图+热图+小提琴图+boxplot(带有p-value)
5. Global的热图,COVID19 vs CTRL DGEs, COVID19 vs TREAT DGEs,TREAT vs CTRL DGEs
6. 每一个亚群的热图,COVID19 vs CTRL DGEs, COVID19 vs TREAT DGEs,TREAT vs CTRL DGEs
7. GO/KEGG找到功能联系的通路，或者通过文章该有什么通路，对应的抽提GSEA的通路的所有基因，绘制热图，然后计算signature score，绘制热图
8. 用通路或者基因绘制Foldchange大图X轴 COVID19 vs CTRL,Y轴COVID19 vs TREAT
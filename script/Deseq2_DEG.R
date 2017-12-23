# ==========================================================
#
#      生信技能树转录组入门课程  -- 差异表达分析
#      •   使用DESeq2进行差异分析
#
# ==========================================================

# 在开始之前保证当前项目目录下面存在以下文件：
# • exprSet.txt (表达矩阵)
# • group_info.txt (设计矩阵)
# • run_DEG.R (Jimmy的差异分析软件,作为参考)
# 文件获取方法：
#   wget  https://raw.githubusercontent.com/jmzeng1314/my-R/master/DEG_scripts/tair/exprSet.txt
#   wget  https://raw.githubusercontent.com/jmzeng1314/my-R/master/DEG_scripts/tair/group_info.txt
#   wget https://raw.githubusercontent.com/jmzeng1314/my-R/master/DEG_scripts/run_DEG.R


# 载入DESeq2和edgeR包
# BiocInstaller::biocLite("DESeq2")
# BiocInstaller::biocLite("edgeR")
library(DESeq2)
library(edgeR)

# 将表达矩阵和设计矩阵导入R
# exprSet_file <-  "C:\\Users\\Shenmy\\Desktop\\2017-12-10\\exprSet.txt"
# group_info_file <- "C:\\Users\\Shenmy\\Desktop\\2017-12-10\\group_info.txt"
exprSet_file <-  "https://raw.githubusercontent.com/shenmengyuan/RNA_seq_Biotrainee/master/data/exprSet.txt"
group_info_file <- "https://raw.githubusercontent.com/shenmengyuan/RNA_seq_Biotrainee/master/data/group_info.txt"
exprSet <- read.table(exprSet_file, stringsAsFactors = F, header = T)
group_info <- read.table(group_info_file, stringsAsFactors = F, header = T)
head(group_info)
group_list <- group_info[,1]
geneLists <- row.names(exprSet)

# 数据过滤
keepGene <- rowSums(edgeR::cpm(exprSet)>0) >=2
table(keepGene);dim(exprSet)
dim(exprSet[keepGene,])
exprSet <- exprSet[keepGene,]
rownames(exprSet) <- geneLists[keepGene]

(colData <- data.frame(row.names=colnames(exprSet), group_list=group_list) )

# DESeq2包可以导入不同的定量软件的输出结果，并且有对应的函数可用
# 参考：https://www.bioconductor.org/help/workflows/rnaseqGene/
# function	        package	           framework      output	              DESeq2 input function
# summarizeOverlaps	GenomicAlignments	R/Bioconductor	SummarizedExperiment	DESeqDataSet
# featureCounts	    Rsubread	        R/Bioconductor	matrix	              DESeqDataSetFromMatrix
# tximport	        tximport	        R/Bioconductor	list of matrices	    DESeqDataSetFromTximport
# htseq-count	      HTSeq	            Python	        files	                DESeqDataSetFromHTSeq

### 构建DESeqDataSet对象
# 到这里可以看到差异分析的开始必须要表达矩阵和设计矩阵这两个文件：
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
### Pre-filtering
# dim(dds)
# 数据过滤也可以在这一步做
# dds <- dds[ rowSums(counts(dds)) > 1 ,]
# dim(dds)

dds <- DESeq(dds) 
plotDispEsts(dds, main="Dispersion plot") 

### Count data transformations
# 为下游的可视化和聚类分析做准备
## rlog()、rlogTransformation()以及与ves()、varianceStabilizingTransformation()的功能是一样的，
## 样本数超过100时，使用vst() 比 rlog() 快. 因为rlog需要拟合每个样本和每个基因的收缩项很耗费时间；

rld <- rlogTransformation(dds)

exprMatrix_rlog <- assay(rld) 
write.table(exprMatrix_rlog,file='DESeq2.exprMatrix_rlog.txt',quote = F,sep = '\t',row.names = T)

# 分别用Day1,2,3与Day0进行比较：
res0vs1 <- results(dds, contrast = c("group_list","Day1","Day0"))
summary(res0vs1,alpha = 0.01)
res0vs2 <- results(dds, contrast = c("group_list","Day2","Day0"))
summary(res0vs2,alpha = 0.01)
res0vs3 <- results(dds, contrast = c("group_list","Day3","Day0"))
summary(res0vs3,alpha = 0.01)

# 按照padj进行排序
res0vs1_order <- res0vs1[order(res0vs1$padj),] 
write.csv(res0vs1_order,"DESeq2_DEG_Day0-Day1_deseq2_results.csv",quote = F,row.names = T) # 把结果输出

# 筛选结果，查看上调和下调基因
res1SigUp <- subset(res0vs1, padj < 0.01 & log2FoldChange >0 ) 
res1SigDown <- subset(res0vs1, padj < 0.01 & log2FoldChange <0 ) 
res2SigUp <- subset(res0vs2, padj < 0.01 & log2FoldChange >0 ) 
res2SigDown <- subset(res0vs2, padj < 0.01 & log2FoldChange <0 ) 
res3SigUp <- subset(res0vs3, padj < 0.01 & log2FoldChange >0 ) 
res3SigDown <- subset(res0vs3, padj < 0.01 & log2FoldChange <0 ) 

sigUpGene <- unique(c(rownames(res1SigUp),rownames(res2SigUp),rownames(res3SigUp)))
length(sigUpGene)

sigDownGene <- unique(c(rownames(res1SigDown),rownames(res2SigDown),rownames(res3SigDown)))
length(sigDownGene)


# 三组比较的差异基因
# 使用UpSetR
# install.packages("UpSetR")
library(UpSetR)
listInput <- list(ah0vs1=c(rownames(res1SigUp),rownames(res1SigDown)),ah0vs2=c(rownames(res2SigUp),
                                                                               rownames(res2SigDown)),ah0vs3=c(rownames(res3SigUp),rownames(res3SigDown)))

png("UpSetR.png")                                            
upset(fromList(listInput), order.by = "freq",text.scale =2,point.size = 5, 
      line.size = 2,mainbar.y.label = "0vs1 0vs2  0vs3", 
      sets.x.label = "Three Time")
dev.off()

# 使用韦恩图
# install.packages("VennDiagram")
library(VennDiagram)
venn.plot<-venn.diagram(
  listInput,
  col = "transparent",
  fill = c("red", "green", "yellow"),
  alpha =c(0.5,0.5,0.5),
  scaled = TRUE,
  ext.text = TRUE,
  ext.line.lwd = 2,
  ext.dist = -0.15,
  ext.length = 0.9,
  ext.pos = -4,
  inverted = TRUE,
  cex = 2.5,
  cat.cex = 2.5,
  filename = "venn.tiff",
  main = "Venn DEG Diagram",
  sub = "0vs1 0vs2  0vs3",
  main.cex = 2,
  sub.cex = 1)

[TOC]

# 实战练习：差异表达分析与功能分析


##### Mengyuan Shen (2017/12/23)

## 0 写在前面

各位小伙伴在生信技能树的两天两夜基础培训课程后，相信大家已经对转录组有了初步的了解。现在我们要开始项目实战，从数据下载、比对、定量到差异分析与功能注释，下面我会一步一步演示给大家看。

## 1 环境配置软件安装

关于环境配置这里不在多说，在之前的课程里我们都有自己处理数据的计算机或者服务器了。现在只需要安装我们需要的软件即可。

```shell
# 创建一个专门安装软件的文件夹
mkdir Biosoft & cd Biosoft 
```

### 1.1 软件列表

-  Miniconda (软件管理器，可一键安装生信软件，类似各种软件管家；)
-  git (可以用来下载GitHub上的软件；管理分享自己代码；windows上安装可以用git bash学习Linux的一些基本操作；)
-  notepad++ (代码编辑器，类似的还有editplus、Sublime)
-  fastqc\RSeQC  (高通量测序质控软件)
-  salmon (不需要比对的定量软件)
-  subread(史上最快的转录组流程，比对+定量)
-  R\Rstudio (统计、画图；用于后续分析与功能分析)

本实战练习要求：在**服务器**上安装：Miniconda/fastqc/salmon/subread 。在**Windows**上安装：git/notepad++/R/Rstudio 。下面是几个软件的安装代码，以供参考：

### 1.2 Minconda

```shell
# 下载安装包（Linux版本）
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# 也可以在清华镜像下载：https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/ 
# wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 直接安装,一路yes即可
bash Miniconda3-latest-Linux-x86_64.sh

# 添加生信软件包下载频道
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

# 下面是清华的频道镜像：https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/
# 据说有时候不稳定，大家可以尝试下：
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
# conda config --set show_channel_urls yes
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/msys2/
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/menpo/
```

上面频道的添加大家也可以**洲更**写的指导说明，也可以自己尝试下，该软件的安装并不难。Minconda3会自带Python3.6，所以大家不用再特地安装Python了。

### 1.3 fastqc

- conda安装

  `conda install -c bioconda fastqc=0.11.5`

- 自己下载安装包安装
  ```shell
  mkdir fastqc &&  cd fastqc
  wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
  # curl: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 
  unzip fastqc_v0.11.5.zip
  # 中间要安装unzip和java，根据系统提示安装，一般是sudo apt install后加软件名。
  ```

### 1.4 subread

Version 1.5.3

跟上面fastqc软件一样，可以使用conda进行subread的安装：`conda install -c bioconda subread`，也可以使用下载安装包的方法进行安装，各有好处，使用安装包安装的好处是知道自己安装软件中有哪些应用，安装最新版本；直接用conda安装则无法判断，虽然他帮你把一些应用加到环境变量里。

### 1.5 Salmon

Salmon v0.8.2

`conda install -c bioconda salmon`

### 1.6 RSeQC

RSeQC 是依赖于python的，直接使用 pip 进行安装: `pip install RSeQC`

官网：[RSeQC: An RNA-seq Quality Control Package](http://rseqc.sourceforge.net/)

中文版可以看：[高通量测序质控及可视化工具包RSeQC](http://mp.weixin.qq.com/s/zKXhy6Vhli1IQdWG8E7uwg)

### 1.7 R包

```R
install.packages("tidyverse") ; library(tidyverse)
install.packages("optparse") ; library(optparse)
install.packages("UpSetR") ; library(UpSetR)
install.packages("rjson") ; library(rjson)

source("https://bioconductor.org/biocLite.R")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
biocLite("DESeq2") ; library(DESeq2)
biocLite('edgeR') ; library(edgeR)
biocLite('limma') ; library(limma)
biocLite("clusterProfiler") ; library(clusterProfiler)
biocLite("DOSE") ; library(DOSE)
biocLite("KEGG.db") ; library(KEGG.db)
biocLite("org.At.tair.db") ; library(org.At.tair.db)
biocLite("pheatmap") ; library(pheatmap)
biocLite("RColorBrewer") ; library(RColorBrewer)
biocLite("AnnotationHub") ; library(AnnotationHub)
biocLite('GenomicFeatures') ; library(GenomicFeatures)
biocLite("tximport") ; library(tximport)
```

一定要提取安装好，保证以上代码正常运行；



## 2 读文章拿到测序数据

数据来自于发表在Nature commmunication上的一篇文章 “[**Temporal dynamics of gene expression and histone marks at the Arabidopsis shoot meristem during flowerin**](https://www.nature.com/articles/ncomms15120)”。原文用RNA-Seq的方式研究在开花阶段，芽分生组织在不同时期的基因表达变化。

实验设计： 4个时间段（0,1,2,3），分别有4个生物学重复，一共有16个样品。

### 2.1 测序数据

```shell
# 创建一个文件夹，用来存放FASTQ文件（公司返回的原始数据）
mkdir -p rna_practice/data/fastq

# 获取样本信息
wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5130/E-MTAB-5130.sdrf.txt
tail -n +2 E-MTAB-5130.sdrf.txt | cut -f 32,36 |sort -u

# 下载数据
# fastq文件下载链接在第几列
# head -n1 E-MTAB-5130.sdrf.txt | tr '\t' '\n' | nl | grep URI
# 根据上述返回数字，获取文件第33列，然后下载fastq文件
# tail -n +2 E-MTAB-5130.sdrf.txt | cut -f 33 | xargs -i wget {}
# 也可以编写批量下载数据的shell脚本，如下：
perl -alne 'if(/.*(ftp:.*gz).*/){print "nohup wget $1 &"}' E-MTAB-5130.sdrf.txt >fq_data_download.sh
bash fq_data_download.sh
```
### 2.2 参考基因组

```powershell
# 创建文件夹用来放置参考基因组或注释文件
mkdir -p rna_practice/data/ref
# 下载cdna（转录本）、dna（基因组）、gff3、gtf文件（注释文件）
nohup wget ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz &
nohup wget ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.28.dna.genome.fa.gz &
nohup wget ftp://ftp.ensemblgenomes.org/pub/plants/release-28/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.28.gff3.gz &
nohup wget ftp://ftp.ensemblgenomes.org/pub/plants/release-28/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.28.gtf.gz &

# 使用nohup &可以将任务放置后台运行，可以关闭远程控制端
```


##  3 序列比对及基因定量

### 3.1 原文的流程

`TopHat -> SummarizeOverlaps -> Deseq2 -> AmiGO`，其中比对的参考基因组为**TAIR10 ver.24** ，并且屏蔽了ribosomal RNA regions (2:3471–9557; 3:14,197,350–14,203,988)。**Deseq2**只计算至少在一个时间段的FPKM的count > 1 的基因。

>Next generation sequencing (NGS) reads were mapped to Arabidopsis thaliana reference transcriptome TAIR10 ver. 24, with ribosomal RNA regions (2:3471–9557; 3:14,197,350–14,203,988) masked, using **TopHat 2.0.13 **(no-mixed alignments; up to 20 secondary alignments; no novel junctions)61. Counts of NGS reads covering transcripts were computed using the function **summarizeOverlaps62** in R. Expressed genes were defined as those having the value of FPKM>1 at least at one time point. Read counts were submitted to differential gene expression analysis in **Deseq2** (default parameters, FDR<0.05)63. Regularized logarithms of read count computed by Deseq2, denoted by rlog, were used for the analysis of relationships between gene expression level and histone modifications signal.



### 3.2 实战流程

![](D:\Biotrainee\rna_seq\RNA_seq_Biotrainee\pic\rnaseq_workflow.png)

#### 3.2.1 Salmon流程

**不需要比对，直接对转录水平进行定量。**

- 创建索引

  `salmon index -t Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz -i athal_index_salmon`

- 定量
```
#! /bin/bash
index=~/rna_seq_practice_2/data/ref/athal_index_salmon ## 指定索引文件夹
quant=~/rna_seq_practice_2/work/quant_salmon # 指定定量文件输出文件夹
for fn in ERR1698{194..209};
do
   sample=`basename ${fn}`
# basename命令用于去掉路径信息，返回纯粹的文件名，如果指定的文件的扩展名，则将扩展名也一并去掉。
   echo "Processin sample ${sampe}"
  salmon quant -i $index -l A \
       -1 ${sample}_1.fastq.gz \
       -2 ${sample}_2.fastq.gz \
       -p 5 -o $quant/${sample}_quant
done
# quant_salmon.sh
```
```nohup bash quant_salmon.sh &```

#### 3.2.2  subread流程

- 创建索引
```
gunzip Arabidopsis_thaliana.TAIR10.28.dna.genome.fa.gz
subread-buildindex -o athal_index_subread   Arabidopsis_thaliana.TAIR10.28.dna.genome.fa
```

- 比对
```
#! /bin/bash
index=~/rna_seq_practice_2/data/ref/athal_index_subread
map=~/rna_seq_practice_2/work/map
for fn in ERR1698{194..209};
do
    sample=`basename ${fn}`
    echo "Processin sample ${sampe}" 
    subjunc -i $index \
        -r ${sample}_1.fastq.gz \
        -R ${sample}_2.fastq.gz \
        -T 5 -o $map/${sample}_subjunc.bam
# 比对的sam自动转为bam，但是并不按照参考基因组坐标排序
done
# map_subjunc.sh
```
打开一个screen窗口运行`bash map_subjunc.sh `，保证当ssh连接断开后服务器也能在后台运行。

- 定量
```
featureCounts=~/anaconda2/bin
# gff3=~/rna_seq_practice_2/data/ref/Arabidopsis_thaliana.TAIR10.28.gff3.gz
gtf=~/rna_seq_practice_2/data/ref/Arabidopsis_thaliana.TAIR10.28.gtf
count=~/rna_seq_practice_2/work/quant_subjunc
nohup $featureCounts/featureCounts  -T 5 -p -t exon -g gene_name -a $gtf -o  $count/counts.txt   *.bam &
nohup $featureCounts/featureCounts  -T 5 -p -t exon -g gene_id -a $gtf -o  $count/counts_id.txt   *.bam &
```

#### 3.2.3 其他比对软件

```shell
# hisat
hisat -p 5 -x $hisat2_mm10_index -1 $fq1 -2 $fq2 -S $sample.sam 2>$sample.hisat.log
samtools sort -O bam -@ 5  -o ${sample}_hisat.bam $sample.sam
# bwa

bwa mem -t 5 -M  $bwa_mm10_index $fq1 $fq2 1>$sample.sam 2>/dev/null 
samtools sort -O bam -@ 5  -o ${sample}_bwa.bam $sample.sam
# bowtie

bowtie -p 5 -x $bowtie2_mm10_index -1 $fq1  -2 $fq2 | samtools sort  -O bam  -@ 5 -o - >${sample}_bowtie.bam
# star

## star软件载入参考基因组非常耗时，约10分钟，也比较耗费内存，但是比对非常快，5M的序列就两分钟即可
star --runThreadN  5 --genomeDir $star_mm10_index --readFilesCommand zcat --readFilesIn  $fq1 $fq2 --outFileNamePrefix  ${sample}_star 
## --outSAMtype BAM  可以用这个参数设置直接输出排序好的bam文件
samtools sort -O bam -@ 5  -o ${sample}_star.bam ${sample}_starAligned.out.sam
```
也可以使用`HTseq`进行计数。

## 4 差异表达分析

### 4.1 设计矩阵和表达矩阵

差异表达分析 

`Deseq2_DEG.R`

• 使用DESeq2进行差异分析 

```R
## 数据过滤
# dds <- dds[ rowSums(counts(dds)) > 1 ,]
# dim(dds)
dds <- DESeq(dds)
plotDispEsts(dds, main="Dispersion plot")
rld <- rlogTransformation(dds)
exprMatrix_rlog=assay(rld)
res <- results(dds, contrast=c("condition",'Day1','Day0'))
resOrdered <- res[order(res$padj),] 
res_Day1_Day0=as.data.frame(resOrdered)
```

### 4.2 一步法差异分析

原教程网址：<https://github.com/jmzeng1314/my-R/tree/master/DEG_scripts>

下载Jimmy的一步法差异分析脚本

```shell
# wget https://raw.githubusercontent.com/jmzeng1314/my-R/master/DEG_scripts/run_DEG.R
# wget https://raw.githubusercontent.com/jmzeng1314/my-R/master/DEG_scripts/tair/exprSet.txt
# wget https://raw.githubusercontent.com/jmzeng1314/my-R/master/DEG_scripts/tair/group_info.txt
# 下载以上表达矩阵、分组矩阵和代码，进行运行：
Rscript run_DEG.R -e exprSet.txt -g group_info.txt -c 'Day1-Day0' -s counts  -m DESeq2
```

## 5 功能分析

差异分析结果的功能注释：

`function_DEG.R`

• 导入差异分析结果，判断显著差异基因 

• 画个火山图看看挑选的差异基因合理与否 

• 显著差异基因的GO/KEGG注释 

• OrgDb类型注释数据学习，了解基因注释原理其实是ID转换

## 6 参考资料

- DESeq2官网说明书：[Analyzing RNA-seq data with DESeq2](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- 生信菜鸟团转录组实战系列：[一个植物转录组项目的实战](http://www.bio-info-trainee.com/2809.html)
- github 转录组教程：[Differential Gene Expression using RNA-Seq (Workflow)](https://github.com/twbattaglia/RNAseq-workflow)
- 生信技能树论坛转录组入门：[RNA-seq基础入门传送门](http://www.biotrainee.com/thread-1750-1-1.html)
- [Jimmy一步法差异分析](https://github.com/jmzeng1314/my-R/tree/master/DEG_scripts)
- [Bioconductor的workflows](http://www.bioconductor.org/help/workflows/)
  - [RNA-Seq workflow: gene-level exploratory analysis and differential expression](http://www.bioconductor.org/help/workflows/rnaseqGene/)
  - [RNA-seq analysis is easy as 1-2-3](http://www.bioconductor.org/help/workflows/RNAseq123/)
  - [Gene Expression Normalization Workflow](http://www.bioconductor.org/help/workflows/ExpressionNormalizationWorkflow/)
  - [Gene-level RNA-seq differential expression and pathway analysis](http://www.bioconductor.org/help/workflows/RnaSeqGeneEdgeRQL/) 
  - [From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline](https://bioconductor.org/help/workflows/RnaSeqGeneEdgeRQL/)
- 《RNA-seq Data Analysis-A Practical Approach》



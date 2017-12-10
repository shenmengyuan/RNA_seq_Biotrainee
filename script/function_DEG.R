# ==========================================================
#
#      生信技能树转录组入门课程  -- 差异分析结果的功能注释
#      •   导入差异分析结果，判断显著差异基因
#      •   画个火山图看看挑选的差异基因合理与否
#      •   显著差异基因的GO/KEGG注释
#      •   OrgDb类型注释数据学习，了解基因注释原理其实是ID转换
#
# ==========================================================

# 在开始之前保证当前项目目录下面存在以下文件：
# • Slamon_0vs1_deseq2_results.csv（差异分析结果）
#   也可以的使用https://raw.githubusercontent.com/jmzeng1314/my-R/master/DEG_scripts/tair/DESeq2_DEG.Day1-Day0.txt


################ 导入差异分析结果，判断显著差异基因 ################
# 为了保证下面课程顺利进行，直接采用GitHub仓库上的差异分析结果
prefix='Day1-Day0'
deg1 <- read.table('https://raw.githubusercontent.com/jmzeng1314/my-R/master/DEG_scripts/tair/DESeq2_DEG.Day1-Day0.txt')
dim(deg1)

# 删除有NA的行，出现NA有以下几个原因：
# Note on p-values set to NA: some values in the results table can be set to NA for one of the following reasons:
# 1. If within a row, all samples have zero counts, the baseMean column will be zero, and the log2 fold change estimates,
#    p value and adjusted p value will all be set to NA.
# 2. If a row contains a sample with an extreme count outlier then the p value and adjusted p value will be set to NA.
#    These outlier counts are detected by Cook’s distance.
# 3. If a row is filtered by automatic independent filtering, for having a low mean normalized count, 
#    then only the adjusted p value will be set to  NA. 

deg1 <- na.omit(deg1)
dim(deg1)
head(deg1)

# 很明显，这个时候用padj来挑选差异基因即可，不需要看foldchange了,padj小于0.05的基因数不多嘛
table(deg1$padj<0.05)
table(deg1$padj<0.01)

# 筛选差异基因，上调，下调
diff_geneList <- rownames(deg1[deg1$padj<0.05,])
up_geneList <- rownames(deg1[deg1$padj<0.05 & deg1$log2FoldChange >0,])
down_geneList <- rownames(deg1[deg1$padj<0.05 & deg1$log2FoldChange <0,])
length(diff_geneList)
length(up_geneList)
length(down_geneList)


################ 画个火山图看看挑选的差异基因合理与否 ###############
library(ggplot2)
library(stringr)

colnames(deg1)

# 这里我不准备用log2FoldChange来挑选差异基因，仅仅是用padj即可,所以log2FoldChange_Cutof设为0
log2FoldChange_Cutof = 0

# 会根据padj和log2FoldChange对基因进行分类：
# padj大于0.05为不显著为NOT;
# log2FoldChange大于0且padj小于0.05为上调UP；
# log2FoldChange小于0且padj小于0.05为下调DOWN；

deg1$change = as.factor(ifelse(deg1$padj < 0.05 & 
                                 abs(deg1$log2FoldChange) > log2FoldChange_Cutof,
                               ifelse(deg1$log2FoldChange > log2FoldChange_Cutof ,'UP','DOWN'),'NOT'))
head(deg1)
summary(deg1)

# 下面开始绘图,本质就是散点图：
this_tile <- paste0('Cutoff for log2FoldChange is ',round(log2FoldChange_Cutof,3),
                    '\nThe number of up gene is ',nrow(deg1[deg1$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(deg1[deg1$change =='DOWN',])
)
g_volcano <- ggplot(data=deg1, aes(x=log2FoldChange, y=-log10(padj), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile  ) + 
  theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
print(g_volcano)
ggsave(filename = paste0(prefix,"_volcano_plot.pdf"),g_volcano)


################ 显著差异基因的GO/KEGG注释 ################
# ============== 首先进行KEGG注释 ================
# source("https://bioconductor.org/biocLite.R")
# biocLite("clusterProfiler")
library(clusterProfiler)
library(org.At.tair.db)

# KEGG API官网：http://www.kegg.jp/kegg/docs/keggapi.html 
# 因为enrichKEGG()函数用的KEGG注释是通过API得到最新的数据，
# KEGG上除了有kegg自己的基因编号还有UniProt/NCBI-ProteinID/NCBI-geneID这些其他数据的ID,
# 如果有这些id,就可以获得注释信息。 
# 以AT1G03000基因为例：http://www.kegg.jp/dbget-bin/www_bget?ath:AT1G03000
# NCBI-GeneID: 839315
# NCBI-ProteinID: NP_171799
# MIPS: AT1G03000.1
# TAIR: AT1G03000
# Araport: AT1G03000
# UniProt: Q8RY16 A0A178WG59
# http://rest.kegg.jp/conv/genes/ncbi-geneid:839315（ncbi-geneid转成kegg-geneid） 
# http://rest.kegg.jp/conv/genes/ncbi-proteinid:NP_171799 (ncbi-proteinid转成kegg-geneid) 
# http://rest.kegg.jp/conv/genes/uniprot:Q8RY16(uniprot转成kegg-geneid)


?enrichKEGG
# 之前差异分析的差异基因名和kegg-geneid一样，所以ketTye默认为kegg，不需要修改；
# 富集分析的原理：http://mp.weixin.qq.com/s/JNVTcbKjtZbXf4XITLk0gg

diff.kk <- enrichKEGG(gene         = diff_geneList,
                      organism     = 'ath',
                      pvalueCutoff = 0.99,
                      qvalueCutoff=0.99
)

# 把geneID转成基因名(gene Symbol),这时候就要用OrgDb类型的数据库了，方法有很多种select函数也可以
kegg_diff_dt <- as.data.frame(setReadable(diff.kk,org.At.tair.db,keytype = 'TAIR'))
# head(kegg_diff_dt)
# View(kegg_diff_dt)

search_kegg_organism('ath', by='kegg_code')

up.kk <- enrichKEGG(gene         = up_geneList,
                    organism     = 'ath',
                    pvalueCutoff = 0.99,
                    qvalueCutoff=0.99
)

kegg_up_dt <- as.data.frame(setReadable(up.kk,org.At.tair.db,keytype = 'TAIR'))

down.kk <- enrichKEGG(gene         = down_geneList,
                      organism     = 'ath',
                      pvalueCutoff = 0.99,
                      qvalueCutoff=0.99
)

kegg_down_dt <- as.data.frame(setReadable(down.kk,org.At.tair.db,keytype = 'TAIR'))

write.csv(x = kegg_diff_dt,file = paste0(prefix,"_kegg_diff.csv"))
write.csv(x = kegg_up_dt,file = paste0(prefix,"_kegg_up.csv"))
write.csv(x = kegg_down_dt,file = paste0(prefix,"_kegg_down.csv"))


# ============== 可视化看看KEGG注释结果 ================
## KEGG patheay visulization: 
# 显著上调下调的差异基因对富集分析结果种的p值进行筛选，小于0.05的注释结果，显著富集，并分组为-1(下调)，1(上调)
down_kegg <- kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1

# 合并筛选出来的通路信息
dat <- rbind(up_kegg,down_kegg)
# View(dat)
colnames(dat)

# 画图看更直观
dat$pvalue <- -log10(dat$pvalue)
dat$pvalue <- dat$pvalue*dat$group 

dat <- dat[order(dat$pvalue,decreasing = F),]

g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
  geom_bar(stat="identity") + 
  scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
  scale_x_discrete(name ="Pathway names") +
  scale_y_continuous(name ="-log10P-value") +
  coord_flip() +
  ggtitle("Pathway Enrichment")
print(g_kegg)
ggsave(filename = paste0(prefix,"_kegg_plot.pdf"),g_kegg)

# 代谢通路图
# browseKEGG(diff.kk, 'ath00500')
# browseKEGG(diff.kk,'ath03010')


# ============== 接着进行GO注释 ================
for (onto in c('CC','BP','MF')){
  
  ego <- enrichGO(gene         = diff_geneList,
                  OrgDb         = org.At.tair.db, 
                  keyType = 'TAIR', 
                  ont           =  onto ,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.2,
                  qvalueCutoff  = 0.9)
  ego <- setReadable(ego, org.At.tair.db,keytype = 'TAIR')
  write.csv(as.data.frame(ego),paste0(prefix,"_",onto,".csv"))
  ego2 <- gofilter(ego,4)
  ego2=ego
  pdf(paste0(prefix,"_",onto,'_barplot.pdf'))
  p=barplot(ego2, showCategory=12)+scale_x_discrete(labels=function(x) str_wrap(x,width=20))
  print(p)
  dev.off()
  
  png(paste0(prefix,"_",onto,'_dotplot.png'))
  p=dotplot(ego2)
  print(p)
  dev.off()

  pdf(paste0(prefix,"_",onto,'_plotGOgraph.pdf'))
  plotGOgraph(ego2)
  dev.off()

  pdf(paste0(prefix,"_",onto,'_enrichMap.pdf'))
  enrichMap(ego2)
  dev.off()
}



################ OrgDb类型注释数据学习，了解基因注释原理其实是ID转 ###############
# 使用Bioconductor的注释包
# biocLite("org.At.tair.db")
# 可以在这个网址找Bioconductor有哪些单独的功能注释包：
# https://bioconductor.org/packages/3.5/data/annotation/
# 官网说明书：http://bioconductor.org/packages/release/data/annotation/manuals/org.At.tair.db/man/org.At.tair.db.pdf
library(org.At.tair.db)
ls('package:org.At.tair.db')

# 用columns()看看数据库里有哪些类型数据
columns(org.At.tair.db)
head(select(org.At.tair.db,columns="SYMBOL",keys = keys(org.At.tair.db)))

# 用everything查询下载下来的org.At.tair.sqlite数据库（234M）在以下路径：
# D:\Software\R-3.4.0\library\org.At.tair.db\extdata
# 直接用R对接该数据库进行了解其中的数据结构
# library(RSQLite)

# 连接数据库
# con <- dbConnect(SQLite(), "D:\\Software\\R-3.4.0\\library\\org.At.tair.db\\extdata\\org.At.tair.sqlite")

# 该数据库中有的表单
# dbListTables(con)

# 上面就是根据数据库里面存在哪些表。它们应该是通过X_id关联着的,以kegg和gene_info为例：
# head(dbReadTable(con, "enzyme"))
# head(dbReadTable(con, "go"))
# head(dbReadTable(con, "genes"))
# head(dbReadTable(con, "kegg"))
# kegg <- dbReadTable(con, "kegg")
# head(dbReadTable(con, "gene_info"))
# gene_infor <- dbReadTable(con, "gene_info")

# 所以进行我们差异基因的注释就是将symol通过X_id关联到go_id和path_id(ID转换)，
# 我们还是用包里的select()函数进行注释吧。

# 其实上一讲的AnnotaionHub除了TxDb类型的注释数据，
# 也能查询下载到拟南芥OrgDb类型注释数据
# library(AnnotationHub)
# ah <- AnnotationHub()
# ath <- query(ah,'thaliana')
# AnnotationHub with 5 records
# # snapshotDate(): 2017-04-25 
# # $dataprovider: UCSC, Inparanoid8, ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
# # $species: Arabidopsis thaliana
# # $rdataclass: TxDb, Inparanoid8Db, OrgDb
# # additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer,
# #   rdatadateadded, preparerclass, tags, rdatapath, sourceurl, sourcetype 
# # retrieve records with, e.g., 'object[["AH10456"]]' 
# 
# title                                     
# AH10456 | hom.Arabidopsis_thaliana.inp8.sqlite      
# AH52245 | TxDb.Athaliana.BioMart.plantsmart22.sqlite
# AH52246 | TxDb.Athaliana.BioMart.plantsmart25.sqlite
# AH52247 | TxDb.Athaliana.BioMart.plantsmart28.sqlite
# AH53758 | org.At.tair.db.sqlite 
# ath_org <- ath[['AH53758']]
# keytypes(ath_org)


## Then draw GO/kegg figures:
deg1$gene_id <- rownames(deg1)
id2symbol <- toTable(org.At.tairSYMBOL) 
deg1 <- merge(deg1,id2symbol,by='gene_id')
head(deg1)

kegg_go_annotion <- select(org.At.tair.db,keys = deg1$gene_id,columns = c("GENENAME","SYMBOL","GO","PATH"),
                           keytype = "TAIR")


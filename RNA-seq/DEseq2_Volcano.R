rm(list = ls())
options(stringsAsFactors = F)
library(cowplot)
library(ggplot2)
library(dplyr)
library(DESeq2)

getwd()
mycounts<-read.csv("sum_count_result.csv",row.names = 1)
mycounts_1<-mycounts[rowSums(mycounts) != 0,]
mymeta<-read.csv("mymeta.csv",stringsAsFactors = T)
id <- read.csv("mymeta.csv",row.names = 1)
id_name <-  rownames(id)
mycounts_2 <- mycounts_1[id_name]
colnames(mycounts_2) == mymeta$id


dds <- DESeqDataSetFromMatrix(countData=mycounts_2, 
                              colData=mymeta, 
                              design=~dex)

dds <- DESeq(dds)

res <- results(dds)

res_1<-data.frame(res)

res_1 %>% 
  mutate(group = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "UP",
    log2FoldChange <= -1 & padj <= 0.05 ~ "DOWN",
    TRUE ~ "NOT_CHANGE"
  )) -> res_2

table(res_2$group)

write.csv(res_2,file="diff_expr_result.csv",
          quote = F)
res_3 <- res_2[,c(2,6)]
colnames(res_3) <- c('log2FoldChange','padj')
write.csv(res_3,file= "DESeq2_diffExpression_logFC.csv",quote = F)

get1 <- res_2 %>% filter( group == "UP" )
get2 <- res_2 %>% filter( group == "DOWN" )
res_4 <- rbind(get1, get2)
write.csv(res_4,file= "DESeq2_diffExpression_logFC_up_and_down.csv",quote = F)


results  <- results(dds)

data <- data.frame(gene = row.names(results),
                   pval = -log10(results$padj), 
                   lfc = results$log2FoldChange)


data <- na.omit(data)

## If fold-change > 1 and pvalue > 1.3 (Increased significant)
## If fold-change < 1 and pvalue > 1.3 (Decreased significant)
data <- mutate(data, color = case_when(data$lfc > 1 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < -1 & data$pval > 1.3 ~ "Decreased",
                                       data$pval < 1.3 ~ "nonsignificant"))
table(data$color)
head(data)
library(ggplot2)
max(data$lfc)
min(data$lfc)
max(data$pval)
min(data$pval)
pdf("Volcano-plot.pdf",width = 6,height= 5)

vol <- ggplot(data, aes(x = lfc, y = pval, color = color))
options(repr.plot.width = 5, repr.plot.height =1)
vol +   
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 1, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "5",
                     values = c(Increased = "#0768AC", Decreased = "#CE3B2A", nonsignificant = "darkgray")) +
  theme_bw(base_size = 14) + 
  theme(legend.position = "right") + 
  xlab(expression(log[2]("R2Tg*/mock"))) + 
  ylab(expression(-log[10]("p-value"))) + 
  geom_hline(yintercept = 1.3, linetype=2, colour = "darkgrey") + geom_vline(xintercept=c(-1,1), linetype=2, color='darkgrey') + xlim(-8,8) + ylim(0,50) + theme(panel.grid=element_blank())  

dev.off()


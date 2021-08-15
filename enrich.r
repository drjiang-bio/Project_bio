rm(list = ls());gc()
library(tidyverse);library(ggplot2);library(cowplot)

###################
# 1. david数据处理
###################
# 载入david数据
getwd()
data_yz <- read.table('./result/enrich/DAVID/Functional Annotation chart_up.txt',
                     sep = '\t', header = T)
data_yf <- data.table::fread('./result/enrich/DAVID/Functional Annotation chart_down.txt',
                            data.table = F)
df_deg_l <- list('UP'=data_yz, 'DOWN'=data_yf)

# 提取数据
# GO
data_yz$Category %>% table
data_go <- lapply(df_deg_l, function(x) {
  filter(x, str_detect(Category, pattern = 'GOTERM'))
})
# kegg
data_kegg <- lapply(df_deg_l, function(x) {
  filter(x, str_detect(Category, pattern = 'KEGG'))
})

for (i in 1:2) {
  data_go[[i]]$TermID <- sapply(strsplit(data_go[[i]]$Term, '~'), function(x) x[1])
  data_go[[i]]$Term <- sapply(strsplit(data_go[[i]]$Term, '~'), function(x) x[2])
  data_go[[i]] <- arrange(data_go[[i]], Category, FDR)
}

for (i in 1:2) {
  data_kegg[[i]]$TermID <- sapply(strsplit(data_kegg[[i]]$Term, ':'), function(x) x[1])
  data_kegg[[i]]$Term <- sapply(strsplit(data_kegg[[i]]$Term, ':'), function(x) x[2])
  data_kegg[[i]] <- arrange(data_kegg[[i]], FDR) %>% .[1:20,]
}

# top10 GO
unique(data_go[[1]]$Category)
data_go_zi <- lapply(data_go, function(x) {
  lapply(unique(x$Category), function(y) {
    x %>% filter(., Category==y) %>% .[1:10,] 
  }) %>% do.call(rbind, .)
})


###################
# 2. 处理成画图所需格式
###################
up <- data_go_zi$UP
up <- up %>% arrange(desc(Category), desc(PValue))
up$Term <- factor(up$Term, levels = up$Term)
id_up <- levels(up$Term)
#
down <- data_go_zi$DOWN
down <- down %>% arrange(desc(Category), desc(PValue))
down$Term <- factor(down$Term, levels = down$Term)
id_down <- levels(down$Term)

n <- length(up$Term) - length(down$Term)
if (n > 0) {
  down$Term <- as.character(down$Term)
  for (i in paste('nn', as.character(seq(1, n, 1)), sep = '')) down <- rbind(list(i, NA, NA), down)
  down$Term <- factor(down$Term, levels = down$Term)
  id_down <- c(rep('', n), id_down)
} else if (n < 0) {
  up$Term <- as.character(up$Term)
  for (i in paste('nn', as.character(seq(1, abs(n), 1)), sep = '')) up <- rbind(list(i, NA, NA), up)
  up$Term <- factor(up$Term, levels = up$Term)
  id_up <- c(rep('', abs(n)), id_up)
}

library(RColorBrewer)
col <- colorRampPalette(brewer.pal(12,'Paired'))(12)
plot(1:12,rep(1,12),col= col,pch=16,cex=2)

#上调 GO
p_up <- ggplot(up, aes(Term, log(PValue, 10), fill=Category)) +
  geom_col(width = 0.7) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent')) +
  theme(axis.line.x = element_line(colour = 'black'), 
        axis.line.y = element_line(colour = 'transparent'), 
        axis.ticks.y = element_line(colour = 'transparent'),
        axis.text = element_text(face = 'plain', size = 14)) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold')) +
  theme(legend.position = 'bottom') +
  coord_flip() +
  geom_hline(yintercept = 0) +
  labs(x = '', y = '', title = 'UP') +
  scale_y_continuous(expand = c(0, 0), breaks = c(-15, -10, -5, 0), 
                     labels = as.character(c(-15, -10, -5, 0))) +     #这儿更改间距设置
  scale_x_discrete(labels = id_up) +
  scale_fill_manual(values=col[c(6,5,7)])
p_up
#下调 GO
p_down <- ggplot(down, aes(Term, -log(PValue, 10), fill=Category)) +
  geom_col(width = 0.7) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent')) +
  theme(axis.line.x = element_line(colour = 'black'), 
        axis.line.y = element_line(colour = 'transparent'), 
        axis.ticks.y = element_line(colour = 'transparent'),
        axis.text = element_text(face = 'plain', size = 14)) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold')) +
  theme(legend.position = 'bottom') +
  coord_flip() +
  geom_hline(yintercept = 0) +
  labs(x = '', y = '', title = 'DOWN') +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 10, 20, 30), 
                     labels = as.character(abs(c(0, 10, 20, 30)))) +       #这儿更改间距设置
  scale_x_discrete(labels = id_down, position = 'top') +
  scale_fill_manual(values=col[c(4,2,1)]) 
p_down

#合并输出pdf
pdf('./result/figure/go_top10_3.pdf', width = 22, height = 12)
plot_grid(p_up, p_down, nrow = 2, ncol = 2, rel_heights = c(9, 1), 
          labels = '-log10(p value)', 
          label_x = 0.5, label_y = 0, label_fontface = 'plain')
dev.off()


# KEGG
up <- data_kegg$UP %>% arrange(desc(PValue))
up$Term <- factor(up$Term, levels = up$Term)
id_up <- levels(up$Term)
#
down <- data_kegg$DOWN %>% arrange(desc(PValue))
down$Term <- factor(down$Term, levels = down$Term)
id_down <- levels(down$Term)

n <- length(up$Term) - length(down$Term)
if (n > 0) {
  down$Term <- as.character(down$Term)
  for (i in paste('nn', as.character(seq(1, n, 1)), sep = '')) down <- rbind(list(i, NA, NA), down)
  down$Term <- factor(down$Term, levels = down$Term)
  id_down <- c(rep('', n), id_down)
} else if (n < 0) {
  up$Term <- as.character(up$Term)
  for (i in paste('nn', as.character(seq(1, abs(n), 1)), sep = '')) up <- rbind(list(i, NA, NA), up)
  up$Term <- factor(up$Term, levels = up$Term)
  id_up <- c(rep('', abs(n)), id_up)
}

library(RColorBrewer)
col <- colorRampPalette(brewer.pal(12,'Paired'))(24)
plot(1:24,rep(1,24),col= col,pch=16,cex=2)

#上调 GO
p_up <- ggplot(up, aes(Term, log(PValue, 10))) +
  geom_col(fill=col[10], width = 0.7) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent')) +
  theme(axis.line.x = element_line(colour = 'black'), 
        axis.line.y = element_line(colour = 'transparent'), 
        axis.ticks.y = element_line(colour = 'transparent'),
        axis.text = element_text(face = 'plain', size = 14)) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold')) +
  theme(legend.position = 'left') +
  coord_flip() +
  geom_hline(yintercept = 0) +
  labs(x = '', y = '', title = 'UP') +
  scale_y_continuous(expand = c(0, 0), breaks = c(-15, -10, -5, 0), 
                     labels = as.character(c(-15, -10, -5, 0))) +     #这儿更改间距设置
  scale_x_discrete(labels = id_up)
p_up
#下调 GO
p_down <- ggplot(down, aes(Term, -log(PValue, 10))) +
  geom_col(fill=col[4],color = col[4], width = 0.7) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent')) +
  theme(axis.line.x = element_line(colour = 'black'), 
        axis.line.y = element_line(colour = 'transparent'), 
        axis.ticks.y = element_line(colour = 'transparent'),
        axis.text = element_text(face = 'plain', size = 14)) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold')) +
  theme(legend.position = 'right') +
  coord_flip() +
  geom_hline(yintercept = 0) +
  labs(x = '', y = '', title = 'DOWN') +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 10, 20, 30), 
                     labels = as.character(abs(c(0, 10, 20, 30)))) +       #这儿更改间距设置
  scale_x_discrete(labels = id_down, position = 'top')
p_down

#合并输出pdf
pdf('./result/figure/kegg_top20.pdf', width = 20, height = 8)
plot_grid(p_up, p_down, nrow = 2, ncol = 2, rel_heights = c(9, 1), 
          labels = '-log10(p value)', 
          label_x = 0.5, label_y = 0, label_fontface = 'plain')
dev.off()


##########调色板附录
library(RColorBrewer)
display.brewer.all(colorblindFriendly = T) #展示所有的色板
display.brewer.pal(11,'RdYIBu')# 展示'Accent'色板中8个颜色
col <- brewer.pal(8,'Accent')# 从Accent 方案中选出8个颜色，并赋值给col
plot(1:8,rep(1,8),col= col,pch=16,cex=2)#绘图看一下颜色


col1<- brewer.pal(8,'Dark2')
col2 <- brewer.pal(8,"Accent")
mycolor <- c(col1[5:8],col2[1:4])#挑选Dark2中后四个，accent中的前四个
plot(1:8,rep(1,8),col= mycolor,pch=16,cex=2)

#函数正常用法：从blue到red，生成10个渐变色
col <- colorRampPalette(c('blue','red'))(10)
plot(1:10,rep(1,10),col= col,pch=16,cex=2)

# 从只有8个颜色的Accent中生成16个颜色
col <- colorRampPalette(brewer.pal(8,'Accent'))(16)
plot(1:16,rep(1,16),col= col,pch=16,cex=2)

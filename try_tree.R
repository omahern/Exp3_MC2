phyo_frac
library(ggtree)
phyo_frac2=subset_taxa(phyo_frac, Strain !="529768d439b63fc1fd2991f7c6ea2071")
phyo_frac3=subset_taxa(phyo_frac2, Strain !="ee279c64d4b8ca7b8abd12858f46f0b7")
phyo_frac
phyo_frac2
phyo_frac3
tree=phy_tree(phyo_frac3)

plotTree(tree)


#write.csv(tax_table(phyo_frac3), 'tree_tips.csv')
## potential chimeras : 529768d439b63fc1fd2991f7c6ea2071
## remove cyano? ee279c64d4b8ca7b8abd12858f46f0b7

tree=phy_tree(phyo_frac3)
dd=read.csv(file='tree_tips2.csv',header=T)


p <- ggtree(tree, layout='rectangular') + geom_tiplab(size=3)
#p <- p %<+% dd + geom_tiplab(aes(color=incorp)) 
geom_tippoint(aes(colour = PATTERN))

p

dat <- data.frame(
  SPECIES = dd$Strain,
  PATTERN = sample(dd$incorp, size = length(dd$Strain), replace = TRUE))

setdiff(tree$tip.label, dd$Strain)

ggtree(tree, layout="rectangular") %<+%
  dat +
  geom_tiplab(offset = .3, hjust = .6, size=2) +
  geom_tippoint(aes(colour = PATTERN)) +
  theme(legend.position = "right")



dat <- data.frame(
  SPECIES = dd$Strain,
  PATTERN = sample(dd$incorp, size = length(dd$Strain), replace = TRUE),
  MAX2= sample(dd$MaxL2FChange, size = length(dd$Strain), replace = TRUE))

p <- ggtree(tree, layout='rectangular') 
p <- p %<+% dd + geom_tiplab(aes(color=incorp),size=2) +
  geom_tippoint(aes(size=MaxL2FChange, color=incorp), alpha=1)
p+theme(legend.position="right")


p <- ggtree(tree, layout='rectangular', aes(color=incorp)) + 
  geom_tiplab2(size=10, color='black')
p

dd=read.csv(file='tree_tips2.csv',header=T)

dat <- data.frame(
  SPECIES = dd$Strain,
  PATTERN = sample(dd$incorp, size = length(dd$Strain), replace = TRUE),
  MAX2= sample(dd$MaxL2FChange, size = length(dd$Strain), replace = TRUE),
  COL= sample(dd$cols, size = length(dd$Strain), replace = TRUE))

ggtree(tree) %<+% dd + aes(color=I(cols))

g=groupOTU(tree, subset)
ggtree(g, aes(color=group)) + scale_color_manual(values=c("black", "#008080"))

ggtree(tree) %<+% dd + geom_highlight(aes(colour=cols))

p=ggtree(tree)
p+geom_tiplab(geom="label", aes(color=))


ggtree(tree) %<+% dat + geom_tiplab(geom="image", aes(color=PATTERN), size=1)

p+geom_hilight(data=dd, mapping=aes(node=tip.label, fill=incorp))



d <- data.frame(node=1:(Nnode(tree)+length(tree$tip.label)), color = c(("black", 'green')))


d <- data.frame(node=1:(Nnode(tree)+length(tree$tip.label)))
# color = c(tip.colours, rep("black", Nnode(aves))))
edge.col[tip.branches] <- tip.colours

                
                
ggtree(tree) %<+% d + aes(color=I(color))


########



p3 <- ggtree(tree, layout="rectangular", branch.length='none') + geom_tiplab(size=2)
p3



read=read.csv(file="subsets_asvs_0.01.csv",header=T,row.names = 1)
data=data.frame(read[,2:6], row.names=read$OTU)
colnames(data)=c("9MC2","11MC2","12MC2","14MC2","15MC2")
phyo_hrsip_f=subset_samples(phyo_hrsip,Comm_Frac=="Frac" )
hrsip_tree=phy_tree(phyo_hrsip_f)
library(viridis)
library(phytools)
c2=c('gray96',(rev(magma(30))))
{phylo.heatmap(hrsip_tree, data, colors=c2, mar=c(1,5,1,1),lwd=2,
               fsize=c(0.7,1,1),split=c(1,0.3),standardize = F)}

ggsim <- ggtree(hrsip_tree) + geom_tiplab(size=2)

gheatmap(p=ggsim, data=data, low="white", high="navy", 
         color='black', colnames=TRUE, colnames_angle = 45)

phylo.heatmap(hrsip_tree, data, colors=c2, mar=c(1,5,1,1),lwd=2,
                 fsize=c(0.7,1,1),split=c(1,0.3),standardize = F)


ggsim <- ggtree(hrsip_tree, layout= 'rectangular',ladderize = TRUE) + 
  geom_tiplab(size=2.5, align = TRUE, hjust = -0.1) +     xlim_tree(0.5)
ggsim
gheatmap(p=ggsim, data=data, offset=0.3,
         font.size = 3, legend_title = "Log2 Fold Change",
         width=0.7,
         color='white', colnames=TRUE, colnames_angle = 45, 
         colnames_position = 'top') + scale_fill_continuous(type='viridis') +
   theme(legend.position="bottom")





p <- ggtree(hrsip_tree) + geom_tiplab(size=2) + theme_tree2()
p
pp <- p %>%
  gheatmap(data, offset=4, width=0.5, colnames=FALSE) %>%
  scale_x_ggtree()
pp + theme(legend.position="right")


read=read.csv(file="/Users/oliviaahern/Documents/GitHub/Exp3_MC2/asvs_subset_l2fchange.csv",header=T)
data1=data.frame(read[,13:18], row.names=read$ASVs)
max2=data.frame(read$ASVs, read[,6])
read=read.csv(file="/Users/oliviaahern/Documents/GitHub/Exp3_MC2/asvs_subset_l2fchange.csv",header=T)
my_subset <- subset(tax_table(phyo_frac), rownames(tax_table(phyo_frac)) %in% read$X.2)
x<-read.csv(file='/Users/oliviaahern/Documents/R/mc2/4Nov21/ASVs/asv-table_2.csv',header=TRUE,row.names=1)

OTU = otu_table(x, taxa_are_rows=T)
phyo_hrsip=phyloseq(OTU, map, my_subset, tree)
tree=phy_tree(phyo_hrsip)




p <- ggtree(tree) + geom_tiplab(size=2, align=TRUE)

p + geom_facet(panel = "Trait", data = max2, geom = geom_col, 
               aes(x = read...4.))

facet_plot(p, panel="Trait", data=max2, mapping=aes(x=max2$read...4., y=1:21), 
           geom=geom_col, orientation='y')
max2=data.frame(read[,6], row.names=read$ASVs)
colnames(max2)=c("value")
#max2=log10(max2)


p+ facet_plot(p, panel="Trait", data=max2, 
              mapping=aes(x=max2$value, y=1:21), fill='black',
           geom=geom_col, orientation='y', width=0.6) + theme_tree2()


data1=data.frame(read[,15:19], row.names=read$ASVs)
data1$pos <- 1:21
p + geom_facet(panel = "SNP", data = data1, geom = geom_point, mapping=aes(x=data1),
               shape = '|') +
  theme_tree2(legend.position=c(.05, .85))




?facet_plot

p + geom_facet(panel = "SNP", data = max2, geom = geom_col, mapping=aes(x=max2$read...4., y=1:21),
              orientation='y') +
  theme_tree2(legend.position=c(.05, .85))


p + geom_facet(panel = "SNP", data = snp_data, geom = geom_point, 
               mapping=aes(x = pos, color = location), shape = '|') +
  geom_facet(panel = "Trait", data = df_bar_data, geom = geom_col, 
             aes(x = dummy_bar_value, color = location, 
                 fill = location), orientation = 'y', width = .6) +
  theme_tree2(legend.position=c(.05, .85))



### good
ggsim <- ggtree(hrsip_tree, layout= 'rectangular',ladderize = TRUE) + 
  geom_tiplab(size=2.5, align = TRUE, hjust = -0.1) +     xlim_tree(0.5)
ggsim
l=gheatmap(p=ggsim, data=data, offset=0.3,
         font.size = 3, legend_title = "Log2 Fold Change",
         width=0.7,
         color='white', colnames=TRUE, colnames_angle = 55, 
         colnames_offset_y = -0.2,
         colnames_position = 'bottom') + scale_fill_continuous(type='viridis') +
  theme(legend.position="bottom")

#data1=data.frame(read[,13:18], row.names=read$ASVs)
max2=data.frame(read[,6], row.names=read$ASVs)
g=ggplot(max2,aes(x=read...6., y=1:21) ) + 
  geom_col(orientation='y') + theme_bw() + 
  labs(x="Log2 Fold Change", y = " ") + 
  theme(axis.text = element_text(color ='black',size=12, hjust=1,vjust=1),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(color='white'),
        axis.title = element_text(color='black',face='bold',size=14),
        legend.position = "none",
        panel.grid=element_blank(),
        strip.text.x = element_text(
          size = 10, color = "black", face = "bold"),
        strip.background = element_rect(
          color="white", fill="white", size=1, linetype="solid"),
        panel.spacing = unit(0.05, "lines"),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.border = element_blank())

g
library(cowplot)
library(ggplot2)
plot_grid(l,g, labels=c("a", "b"), rel_widths = c(2,1),
          align='hv')
library(aplot)
plot_list(l,g, ncol=2, tag_levels='A')
maxxx=max2$read...6.

x= l + geom_facet(panel='blah', data = max2, geom= geom_col, 
                  mapping= aes(x=maxxx, y=1:21), orientation='y',
                  width=0.8) +
  xlim_tree(5)
x

x + xlim_tree(6) + xlim_expand(c(-10, 10), 'Dot')




ggsim <- ggtree(hrsip_tree, layout= 'rectangular',ladderize = TRUE) + 
  geom_tiplab(size=2.5, align = TRUE, hjust = -0.1) +     xlim_tree(0.5)
ggsim
l= gheatmap(p=ggsim, data=data1, offset=0.3,
           font.size = 3, legend_title = "Log2 Fold Change",
           width=0.7,
           color='white', colnames=TRUE, colnames_angle = 55, 
           colnames_offset_y = -0.2,
           colnames_position = 'bottom') + 
  theme(legend.position="bottom")
l

f=ggsim +  geom_facet(panel='blah', data = max2, geom= geom_col, 
  mapping= aes(x=maxxx, y=1:21), orientation='y',width=0.8) 
s=facet_widths(f, widths = c(1, 2))

f
plot_list( l, f)


df_l2ft=read.csv('/Users/oliviaahern/Documents/GitHub/Exp3_MC2/data_greaterthan0.1.csv',header=T,row.names=1)
sub=subset(df_l2ft, Good=="Yes")
d=dcast(sub, OTU ~ timepoint, value.var = "log2FoldChange", fun.aggregate = sum)
read=read.csv(file="/Users/oliviaahern/Documents/GitHub/Exp3_MC2/asvs_subset_l2fchange.csv",header=T)
data1=data.frame(read[,13:18], row.names=read$ASVs)
max2=data.frame(read[,6], row.names=read$ASVs)


library(reshape2)
#data1=data.frame(read[,9:13], row.names=read$ASVs)
heats=d[,2:6]
rownames(heats)=d[,1]
vec=rownames(max2)
data_new1 <- heats[match(rev(vec), rownames(heats)), ]  
data_new1

#rownames(data_new1)=21:1
colnames(data_new1)=c("9MC2", "11MC2", "12MC2", "14MC2", "15MC2")
melt=melt(as.matrix(data_new1))
colnames(melt)=c("id", "Treat", "Max2")

heatnum=data_new1
rownames(heatnum)=21:1
melt2=melt(as.matrix(heatnum))
colnames(melt2)=c("id", "Treat", "Max2")

viv=viridis::viridis(50)
viv=viridis::cividis(50)
l= gheatmap(p=ggsim, data=data_new1, offset=0.3,
            font.size = 3, legend_title = "Log2 Fold Change",
            width=0.4, low = viv[1], high= viv[50],
            color='white', colnames=TRUE, colnames_angle = 55, 
            colnames_offset_y = -0.2,
            colnames_position = 'bottom') + 
  theme(legend.position="bottom") 
l

colnames(max2)=c("id", "maxxx")

sh=facet_plot(p=ggsim,panel='Maximum Log2 Fold Change', data=max2, geom=geom_col,
                  mapping= aes(x=maxxx),orientation='y',
                  width=1, color='white', fill='gray80') + 
  labs(x="Max Log2 Fold Change")
sh
?facet_grid



sh2= facet_plot(p=ggsim,panel='asdf', data=max2, geom=geom_col,
              mapping= aes(x=maxxx),orientation='y',
              width=1, color='black', fill='gray80') 
sh2
### OK
sh3=gheatmap(p=ggsim, data=data_new1, offset=0.3,
           font.size = 3, legend_title = "Log2 Fold Change",
           width=0.4, low = viv[1], high= viv[50],
           color='white', colnames=TRUE, colnames_angle = 55, 
           colnames_offset_y = -0.2,
           colnames_position = 'bottom') + 
  theme(legend.position="bottom")  +
  facet_plot(p=ggsim,panel='Log2', data=max2, geom=geom_col,
             mapping= aes(x=maxxx),orientation='y',
             width=0.7, color='white', fill=viv[1])
sh4 = sh3 + xlim_tree(6)  +
 # geom_text(aes(label=lab), data=max2) + 
  theme(plot.margin=margin(6, 6, 40, 6)) + 
 scale_x_continuous(name=c("Max Log2 Fold Change"))
sh4

g=ggplot(max2,aes(x=maxxx, y=1:21) ) + 
  geom_col(orientation='y') + theme_bw() + 
  labs(x="Log2 Fold Change", y = " ") + 
  theme(axis.text = element_text(color ='black',size=12, hjust=1,vjust=1),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(color='white'),
        axis.title = element_text(color='black',face='bold',size=14),
        legend.position = "none",
        panel.grid=element_blank(),
        strip.text.x = element_text(
          size = 10, color = "black", face = "bold"),
        strip.background = element_rect(
          color="white", fill="white", size=1, linetype="solid"),
        panel.spacing = unit(0.05, "lines"),
        plot.margin = margin(1,1,1,1, "cm"))
g

plot_grid(l, g, labels=c("a", "b"), rel_widths = c(2,1),
          axis='t',
          align='hv')
library(aplot)
plot_list(l,g, ncol=2, tag_levels='A')

## try w heatmap and them off to sides 
sh3=gheatmap(p=ggsim, data=data_new1, offset=0.3,
             font.size = 3, legend_title = "Log2 Fold Change",
             width=0.4, low = viv[1], high= viv[50],
             color='white', colnames=TRUE, colnames_angle = 55, 
             colnames_offset_y = -0.2,
             colnames_position = 'bottom') + 
  theme(legend.position="bottom")
sh3
asvs=read$ASVs
maxsx=read$MaxL2FChange
peak=log(read$peak)
try1=data.frame(read[,2:7], row.names=asvs)
p <- ggtree::ggtree(tree) 
p[["data"]]$label

ggsim <- ggtree(hrsip_tree) + geom_tiplab(size=2)
s1=s[['data']]
try=s1$angle
try2=s1$label
orderme=data.frame(try[1:18], row.names=try2[1:18])
setdiff(asvs, try2)




---
title: "Trophic Cascades - MC2"
output: 
  html_document:
    keep_md: true
    toc: true
    toc_float: true
    toc_depth: 5
    code_folding: hide
    number_sections: true
    theme: cosmo

knit: (function(input_file, encoding) {
  out_dir <- 'docs';
  rmarkdown::render(input_file,
 encoding=encoding,
 output_file=file.path(dirname(input_file), out_dir, 'index.html'))})
---


![Experiment 3 - MC2 Timeline](/Users/ahern/downloads/Exp3_MC2_Timeline.png)


These ASV tables were generated using Dada2/Qiime2. Full code available here:  https://drive.google.com/file/d/1xo94K0Stehiay85RZSoHhuamBjvUQDDB/view?usp=sharing.


The goal of this document is to look at the 16S rRNA data from Experiment 3 replicate 2. First I will try to:

1. Look at the state data from the chemostat 
2. Determine whether the 16S community stays stable over time in the chemostat
  + alpha diversity 
    + general sequencing statistics
    + relative abundance over time
  + beta diversity indices (compositional)
    + PhILR phylogenetic isometric log ratio transformation
    + CLR centered log ratio transformation 
3. Find certain ASVs that are present in heavier SIP fractions
  + run with both ASVs and Genus level taxonomic ranks 
4. Bring it all together by highlighting the SIP ASVs

# State Data

## Ph/DO Probe

[pH/DO Probe: VisiFerm 325 interfaced with Hamilton software](https://www.hamiltoncompany.com/process-analytics/sensors)

```r
read=read.csv(file="/Users/ahern/R/trophic_cascades/mc2/m2_do_02.csv",header=T)
read=subset(read, T<275)
time=read$T
ph= read$pH.M2
do=read$DO.M2..uM.
{
  par(mar=c(5,5,1,5))
  plot(time, do,type='l',col='blue',xlab="Time (hours)", ylab = "")
  rect(xleft=0, xright=240, ybottom=111, ytop=281, col= rgb(0.7,0.7,0.7,alpha=0.2),
    border=NA)
  mtext(side=2,"Dissolved Oxygen (uM)", col = 'blue', padj=-5)
  par(new = TRUE)
  plot(time, ph,col='red',type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
  mtext(side=4, "pH", col ='red', line=3)
  axis(side=4, at = c(6.5, 6.7, 7, 7.2, 7.4))

  }
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

pH and DO over time with the 13C addition highlighted. 

## Carbon Concentrations 

Generated using some sort of coloration with a plate reader? Have to find kit info. 


```r
read=read.csv(file="/Users/ahern/R/trophic_cascades/mc2/states.csv",header=T)

{par(mfrow=c(2,3),mar=c(2,6,3,1))
plot(read$Time,read$ethanol_um,type='o',ylim=c(0,22),ylab="Concentration\n (um)", col='gray30',bg='#ba4f45',pch=21,cex.axis=1.3,cex.lab=1.2,cex=1.3, main="Ethanol")
plot(read$Time,read$methanol_um,type='o',ylim=c(0,20),ylab="Concentration\n (um)", col='gray30',bg='#ad993c',pch=21,cex.axis=1.3,cex.lab=1.2,cex=1.3, main="Methanol")
plot(read$Time,read$Acetate_um,type='o',ylim=c(0,12),xlab= "Time (hours)", ylab = "Concentration\n (um)",col='gray30',bg='#56ae6c',pch=21,cex.axis=1.3,cex.lab=1.2,cex=1.3, main="Acetate")
plot(read$Time,read$Xylose_um,type='o',ylim=c(0,12),xlab= "Time (hours)", ylab = "Concentration\n (um)",col='gray30',bg='#7065bb',pch=21,cex.axis=1.3,cex.lab=1.2,cex=1.3, main="Xylose")
plot(read$Time,read$Glucose_um,type='o',ylim=c(0,12),xlab= "Time (hours)", ylab = "Concentration\n (um)",col='gray30',bg='#b65090',pch=21,cex.axis=1.3,cex.lab=1.2,cex=1.3,main="Glucose")}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-2-1.png)<!-- -->


```r
read=read.csv(file="/Users/ahern/R/trophic_cascades/mc2/states.csv",header=T)

{par(mfrow=c(5,1),mar=c(2,7,3,1),xpd=T)
plot(read$Time,read$ethanol_um,type='o',ylim=c(0,22),ylab="Chemostat\n Concentration\n (um)", col='gray30',bg='#ba4f45',pch=21,cex.axis=1.3,cex.lab=1.2,cex=1.3)
text(-200,25,"Ethanol Feed (um)")
segments(-144,0,-144,23,col='gray50',lty=2)
segments(-60,0,-60,23,col='gray50',lty=2)
text(-62,25,"49.9")
segments(-44,0,-44,23,col='gray50',lty=2)
text(-44,25,"18.6")
segments(-0.5,0,-0.5,23,col='gray50',lty=2)
text(-0.5,25,"31.7")
segments(120,0,120,23,col='gray50',lty=2)
text(120,25,"16")
segments(240,0,240,23,col='gray50',lty=2)
text(240,25,"53.2")

plot(read$Time,read$methanol_um,type='o',ylim=c(0,22),ylab="Chemostat\n Concentration\n (um)", col='gray30',bg='#ad993c',pch=21,cex.axis=1.3,cex.lab=1.2,cex=1.3)
text(-200,25,"Methanol Feed (um)")
segments(-144,0,-144,23,col='gray50',lty=2)
segments(-60,0,-60,23,col='gray50',lty=2)
text(-65,25,"45.01")
segments(-44,0,-44,23,col='gray50',lty=2)
text(-44,25,"18.71")
segments(-0.5,0,-0.5,23,col='gray50',lty=2)
text(-0.5,25,"25.93")
segments(120,0,120,23,col='gray50',lty=2)
text(120,25,"20.05")
segments(240,0,240,23,col='gray50',lty=2)
text(240,25,"48.09")

plot(read$Time,read$Acetate_um,type='o',ylim=c(0,12),xlab= "Time (hours)", ylab = "Chemostat\n Concentration\n (um)",col='gray30',bg='#56ae6c',pch=21,cex.axis=1.3,cex.lab=1.2,cex=1.3)
text(-200,14,"Acetate Feed (um)")
segments(-144,0,-144,12.5,col='gray50',lty=2)
text(-144,14,"154.36")
segments(-60,0,-60,12.5,col='gray50',lty=2)
text(-65,14,"102.14" )
segments(-44,0,-44,12.5,col='gray50',lty=2)
text(-40,14, "101.58" )
segments(-0.5,0,-0.5,12.5,col='gray50',lty=2)
text(-0.5,14,"167.42")
segments(120,0,120,12.5,col='gray50',lty=2)
text(120,14,"85.19")
segments(240,0,240,12.5,col='gray50',lty=2)
text(240,14,"123.53")

plot(read$Time,read$Xylose_um,type='o',ylim=c(0,12),xlab= "Time (hours)", ylab = "Chemostat\n Concentration\n (um)",col='gray30',bg='#7065bb',pch=21,cex.axis=1.3,cex.lab=1.2,cex=1.3)
text(-200,14,"Xylose Feed (um)")
segments(-144,0,-144,12.5,col='gray50',lty=2)
text(-144,14,"58.86")
segments(-60,0,-60,12.5,col='gray50',lty=2)
text(-62,14,"42" )
segments(-44,0,-44,12.5,col='gray50',lty=2)
text(-44,14, "37.43" )
segments(-0.5,0,-0.5,12.5,col='gray50',lty=2)
text(-0.5,14,"57.71")
segments(120,0,120,12.5,col='gray50',lty=2)
text(120,14,"24.29")
segments(240,0,240,12.5,col='gray50',lty=2)
text(240,14,"35.71")
plot(read$Time,read$Glucose_um,type='o',ylim=c(0,12),xlab= "Time (hours)", ylab = "Chemostat\n Concentration\n (um)",col='gray30',bg='#b65090',pch=21,cex.axis=1.3,cex.lab=1.2,cex=1.3)
text(-200,14,"Glucose Feed (um)")
segments(-144,0,-144,12.5,col='gray50',lty=2)
segments(-60,0,-60,12.5,col='gray50',lty=2)
segments(-44,0,-44,12.5,col='gray50',lty=2)
segments(-0.5,0,-0.5,12.5,col='gray50',lty=2)
segments(120,0,120,12.5,col='gray50',lty=2)
segments(240,0,240,12.5,col='gray50',lty=2)
text(240,14,"35.39")}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

Concentrations of different carbon substrates in the inflow. Lines represent sampling times for feed concentrations. Feed concentrations on top. Missing values represent no data. Unless there is the same feed for both MC1 and MC2? I'm not sure...

## Delta 13C

```r
read=read.csv(file="/Users/ahern/R/trophic_cascades/mc2/states.csv",header=T)
glucose=read$Glucose_um
time_glucose=read$Time

read=read.csv(file="/Users/ahern/R/trophic_cascades/mc2/c13_glucose.csv",
              header=T)


{par(mfcol=c(2,1),mar=c(4,5.5,1,1),xpd=T)
plot(time_glucose,glucose,type='o',ylim=c(0,12),xlab= "Time (hours)", ylab = "Glucose\n (um)",col='gray30',bg='#b65090',pch=21,cex.axis=1.3,cex.lab=1.2,cex=1.3)
  plot(read$T,read$Atom..13C,type='o',ylim=c(0,24),xlab="Time (hours)",
     ylab="Atom C13\n (%)", col='gray20',bg='gray50',
     pch=21,cex.axis=1.3,cex.lab=1.2,cex=1.3)
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

20% enrichment of 13C occured around 24 hours (24h = 19%)


## Buoyant Density 
Buoyant density data vs. r max. RT-qPCR primers were the universal bacterial 341F and 806R. 

```r
data<-read.csv("/Users/ahern/R/trophic_cascades/mc2/26oct21/ASV_run3/bd_data.csv",header=T,row.names=1)

summary(data$TP)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   -0.45   24.00   48.00   79.12  120.00  240.00
```

```r
MC5=data[1:10,]
MC9=data[11:23,]
MC11=data[24:32,]
MC12=data[33:41,]
MC14=data[42:51,]
MC15=data[52:60,]
library(viridis)
col=viridis(6)

{plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=21,
      bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum")
lines(MC9$BD,MC9$R_MAX,type='o',pch=23,
      bg=col[2])
lines(MC11$BD,MC11$R_MAX,type='o',pch=21,
      bg=col[3])
lines(MC12$BD,MC12$R_MAX,type='o',pch=21,
      bg=col[4])
lines(MC14$BD,MC14$R_MAX,type='o',pch=21,
      bg=col[5])
lines(MC15$BD,MC15$R_MAX,type='o',pch=21,
      bg=col[6])
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

Buoyant density vs ration of maxiumum broken up by timepoint. 

```r
## focus on two tiempoints first



{par(mfrow=c(1,5), mar=c(2,2,1,1))
  plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
      bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum")
lines(MC9$BD,MC9$R_MAX,type='o',pch=21,
      bg=col[2])
legend('topright', legend=c('5MC2','9MC2'), pt.bg=c(col[1:2]),
       pch=c(23,21),bty='n')

plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
      bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum")
lines(MC11$BD,MC11$R_MAX,type='o',pch=21,
      bg=col[3])
legend('topright', legend=c('5MC2','11MC2'), pt.bg=c(col[1], col[3]),
       pch=c(23,21),bty='n')

plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
      bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum")
lines(MC12$BD,MC12$R_MAX,type='o',pch=21,
      bg=col[4])
legend('topright', legend=c('5MC2','12MC2'), pt.bg=c(col[1], col[4]),
       pch=c(23,21),bty='n')

plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
      bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum")
lines(MC14$BD,MC14$R_MAX,type='o',pch=21,
      bg=col[5])
legend('topright', legend=c('5MC2','14MC2'), pt.bg=c(col[1], col[5]),
       pch=c(23,21),bty='n')

plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
      bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum")
lines(MC15$BD,MC15$R_MAX,type='o',pch=21,
      bg=col[6])
legend('topright', legend=c('5MC2','15MC2'), pt.bg=c(col[1], col[6]),
       pch=c(23,21),bty='n')
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-6-1.png)<!-- -->



# Bacterial Community


```r
#setwd('/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/')
library(HTSSIP)
library(phyloseq)
library(ggplot2)
library(ape)
library(vegan)
library(tidyverse)
library(HTSSIP)
# load data in 
x<-read.csv(file='/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/asv-table.csv',header=TRUE,row.names=1)
OTU = otu_table(x, taxa_are_rows=T)
taxa<-read.csv(file="/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/taxonomy_4nov21.csv",header=TRUE,row.names=1)
t<-as.matrix(taxa)
tax2<-tax_table(t)
map<-import_qiime_sample_data("/Users/ahern/R/trophic_cascades/mc2/26oct21/ASV_run3/map.txt")
tree=read.tree("/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/tree.nwk")

phyo = phyloseq(OTU, tax2,map,tree)
phyo
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1026 taxa and 70 samples ]
## sample_data() Sample Data:       [ 70 samples by 16 sample variables ]
## tax_table()   Taxonomy Table:    [ 1026 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 1026 tips and 1025 internal nodes ]
```

```r
phyo1 = subset_taxa(phyo, !Order=="Chloroplast")
phyo1 # 91 chloroplasts
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 935 taxa and 70 samples ]
## sample_data() Sample Data:       [ 70 samples by 16 sample variables ]
## tax_table()   Taxonomy Table:    [ 935 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 935 tips and 934 internal nodes ]
```

```r
phyo2 = subset_taxa(phyo1, !Family=="Mitochondria")
phyo2 # 12 mitochondria
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 923 taxa and 70 samples ]
## sample_data() Sample Data:       [ 70 samples by 16 sample variables ]
## tax_table()   Taxonomy Table:    [ 923 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 923 tips and 922 internal nodes ]
```
## Alpha Diversity 

### No of Sequences per Sample 

```r
{par(mar=c(5,5,1,1))
barplot(colSums(OTU),las=2,horiz=T,cex.names=0.5)
box(which='plot')}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/no-asvs-1.png)<!-- -->

```r
# samples with lowest number of sequences 
head(sort(colSums(OTU)),decreasing=F)
```

```
##  MCP_9MC205  MCP_9MC203 MCP_11MC206 MCP_12MC206 MCP_15MC214 MCP_14MC206 
##       10541       12485       18938       19525       23009       25324
```

### Basic Alpha Diversity Measurements
Assessing alpha diversity using species number and Pielou's Evenness Index (ranges from 0 - uneven to 1 totally even).

```r
# calculate no. asvs
s=specnumber(t(OTU))

# calculate Pielou's evenness index 
H <- diversity(t(OTU))
J <- H/log(s)

{par(mar=c(4,4,1,1),mfrow=c(1,2))
b=barplot(s,las=2,horiz=T, col='gray50',xlab='No. ASVs', yaxt='n')
box(which='plot')
mtext(side=2, at=b, text=colnames(OTU),las=2,cex=0.5,adj=1.1)


barplot(J, las=2,horiz=T,yaxt='n', col='gray30',
        xlim=c(0,1),xlab="Pielous Index")
box(which='plot')
mtext(side=2, at=b, text=colnames(OTU),las=2,cex=0.5,adj=1.1)
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/alpha-div-plots-1.png)<!-- -->

### Class Level Relative Abundance   

```r
# transform counts into relative abudance 
phyo_abund=transform_sample_counts(phyo2, function(x) x / sum(x) )

## I'm not really good at ggplot... So I use a different style for plotting 
class=tax_glom(phyo_abund, taxrank="Class")
classy=as.matrix((otu_table(class)))
tax_class=data.frame(tax_table(class))
colors=read.csv('/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/colors_923.csv',row.names=1,header=T)
#colors=randomcoloR::randomColor(50, hue="bright") # used a random color generator for colors 

{par(mar=c(4,4,1,7),xpd=T)
barplot(classy, las=2, cex.names = 0.5, col=colors$x,
        ylab="Relative Abundance Bacterial Classes", space=0)
box(which='plot')
legend(73,1, legend=tax_class$Class,cex=0.6, pt.bg=colors$x, pch =22,col=colors$x,
       bty='n', ncol=1)
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/class-abund-all-asvs-1.png)<!-- -->

```r
c1=tax_glom(phyo2, taxrank="Class")
f1=sort( ((rowSums(otu_table(c1))/sum(otu_table(c1)))*100),decreasing=T)
# head(f1)

# my_subset <- subset(tax_table(c1), rownames(tax_table(c1)) %in% c('cf3381318906e1e08f58bca905636dc0', 'efbe1f58b1e2984ddc53a64f047d94ff', '63413adbe8227f80412e824dffac1f5a'))
```


Top Bacterial Classes are: 

1. Alphaproteobacteria 52.39%
2. Gammaproteobacteria 35.99%
3. Bacteroidia 8.60%


### ASV Level Relative Abundance All    

```r
asvy=as.matrix((otu_table(phyo_abund)))

{par(mar=c(5,4,1,1),xpd=T)
barplot(asvy, las=2, cex.names = 0.3, col=colors$x,
        ylab="Relative Abundance Bacterial ASVs",space=0)
box(which='plot')
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/asv-abund-all-asvs-1.png)<!-- -->

```r
# top ASVs
f1=sort( ((rowSums(otu_table(phyo2))/sum(otu_table(phyo2)))*100),decreasing=T)
head(f1)
```

```
## cf3381318906e1e08f58bca905636dc0 efbe1f58b1e2984ddc53a64f047d94ff 
##                        44.467291                        18.480045 
## 1e104b993416995498346089f57fd3bf 63413adbe8227f80412e824dffac1f5a 
##                         8.327097                         8.078652 
## 1583d8af9fd9cb078f88c0e1a8b8a888 83c246c4ed456dea2a15d054f16009b3 
##                         3.060769                         2.808620
```

```r
my_subset <- subset(tax_table(phyo2), rownames(tax_table(phyo2)) %in% c('cf3381318906e1e08f58bca905636dc0', 'efbe1f58b1e2984ddc53a64f047d94ff', '1e104b993416995498346089f57fd3bf',"63413adbe8227f80412e824dffac1f5a","1583d8af9fd9cb078f88c0e1a8b8a888","83c246c4ed456dea2a15d054f16009b3"))
```

Top 6 ASVs:

1. Alphaproteobacteria / Rhizobiales / Rhizobiaceae / Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium ASV cf3381 (44.47%)
2. Gammaproteobacteria / Betaproteobacteriales / Burkholderiaceae / Janthinobacterium ASV efbe1f (18.48%)
3. Gammaproteobacteria / Betaproteobacteriales/ Methylophilaceae / Methylophilus ASV 1e104b9 (8.33%)
4. Bacteroidia / Sphingobacteriales / NS11-12 marine group / ASV 63413ad (8.08%)
5. Alphaproteobacteria / Rhizobiales / Xanthobacteraceae / Ancylobacter ASV 1583d8a (3.06%)
6. Gammaproteobacteria / Betaproteobacteriales / Rhodocyclaceae / Zoogloea ASV 83c246


### ASV Level Relative Abundance Fraction  
Higher buoyant density in the lower fractions. This is all of the ASVs in all fractions. 

```r
phyo_abundf=subset_samples(phyo_abund, Comm_Frac=="Frac")
asvy=as.matrix((otu_table(phyo_abundf)))

space=c(rep("0",10),1,rep("0",12),1,rep("0",8),1,rep("0",8),1,rep("0",9),1,(rep("0",8)))
{par(mar=c(5,4,2,1),xpd=T)
barplot(asvy, las=2, cex.names = 0.6, col=colors$x,
        ylab="Relative Abundance Bacterial ASVs",space=as.numeric(space))
box(which='plot')
text(5,1.05,"5MC2 -0.45h",font=1)
text(17,1.05, "9MC2 24h",font=1)
text(29,1.05,"11MC2 48h",font=1)
text(40, 1.05,'12MC2 72h',font=1)
text(50,1.05, "14MC2 120h",font=1)
text(60,1.05,"15MC2 240h",font=1)
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/asv-abund-all-asvs-frac-1.png)<!-- -->

### ASV Level Relative Abundance Community   
Just community samples  

```r
phyo_abundf=subset_samples(phyo_abund, Comm_Frac=="Community")
asvy=as.matrix((otu_table(phyo_abundf)))

{par(mar=c(7,4,2,1),xpd=T)
barplot(asvy, las=2, col=colors$x, space=0,ylab="Relative Abundance Bacterial ASVs")  
box(which='plot')
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/asv-abund-all-asvs-comm-1.png)<!-- -->



### Class Community Top Asvs 
Just community samples  

```r
phyo_abundf=subset_samples(phyo2, Comm_Frac=="Community")

# Compile taxa by Order (filtering out low abundance taxa)
physeq2 = filter_taxa(phyo_abundf, function(x) mean(x) > 0.1, TRUE)
physeq2
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 462 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 16 sample variables ]
## tax_table()   Taxonomy Table:    [ 462 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 462 tips and 461 internal nodes ]
```

```r
physeq3 = transform_sample_counts(physeq2, function(x) x / sum(x) )
physeq3
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 462 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 16 sample variables ]
## tax_table()   Taxonomy Table:    [ 462 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 462 tips and 461 internal nodes ]
```

```r
stom <- psmelt(physeq3)
glom <- tax_glom(physeq3, taxrank = 'Class')

data <- psmelt(glom) # create dataframe from phyloseq object
data$Class <- as.character(data$Class) #convert to character

#simple way to rename phyla with < 1% abundance
data$Class[data$Abundance < 0.01] <- "< 1% abund."

library(plyr)
medians <- ddply(data, ~Class, function(x) c(median=median(x$Abundance)))
remainder <- medians[medians$median <= 0.01,]$Class

data[data$Class %in% remainder,]$Class <- "Phyla < 1% abund."
#rename phyla with < 1% relative abundance
data$Class[data$Abundance < 0.01] <- "Phyla < 1% abund."

is.numeric(data$Abundance)
```

```
## [1] TRUE
```

```r
library(reshape2)
d=dcast(data, Class ~ ID, value.var = "Abundance", fun.aggregate = sum)

{par(mar=c(7,4,2,10),xpd=T)
barplot(as.matrix(d[,2:11]), las=2, space=0,ylab="Relative Abundance Bacterial Classes",
        col=c('#d86624','#693dc6','#be6bf9','#6af76a','#e83e7c','gray70','#0eb55e'))  
box(which='plot')
legend(10.5,0.6,pt.bg=c('#d86624','#693dc6','#be6bf9','#6af76a','#e83e7c','#0eb55e','gray70'),pch=22,legend=c('Actinobacteria','Alphaproteobacteria','Bacteroidia','Deltaproteobacteria','Gammaproteobacteria','Verrucomicrobiae','Less than 1%'),bty=
         'n')
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/class-abund-top-ones-comm-1.png)<!-- -->

## Beta Diversity 
16S data is inherently compositional. Therefore I am re-analyzing the data in terms of a compositional approach. 

[Compositional data analysis paper](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full)

### PhILR
Phylogenetic approach using algorithms appropriate for compositional datasets. [PhILR: Phylogenetic Isometric Log Ratio Transform; Silverman et al 2017](https://elifesciences.org/articles/21887)

#### All samples 

```r
library(philr)
GP <- transform_sample_counts(phyo2, function(x) x+1)
#
GP
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 923 taxa and 70 samples ]
## sample_data() Sample Data:       [ 70 samples by 16 sample variables ]
## tax_table()   Taxonomy Table:    [ 923 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 923 tips and 922 internal nodes ]
```

```r
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "Species_/Kingdom_Bacteria"
```

```r
otu.table <- t(otu_table(GP))
tree <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')


gp.dist <- dist(gp.philr, method="euclidean")
gp.pcoa <- ordinate(GP, 'PCoA', distance=gp.dist)
{par(mar=c(5,5,1,8),xpd=T)
  plot(gp.pcoa$vectors[,1],gp.pcoa$vectors[,2],pch=sample_data(GP)$pch,
      bg=as.character(sample_data(GP)$col),cex=2, 
      xlab = "PCoA1 61.62%", ylab= "PCoA2 21.31%")
  ordiellipse(gp.pcoa$vectors,
              group=as.factor(sample_data(GP)$Treatment),
              label=T,xpd=F,
              kind ='sd',conf=0.8)
  legend(27,15,legend=c('Community',"Fraction","5MC2","6MC2","7MC2","8MC2",
                          "9MC2","11MC2","12MC2","14MC2","15MC2"),
       pt.bg=c('white','white','#b94c3f','#d08f36','#a67e3b','#99af3e','#62a352','#45c097','#697cd4','#9350a1','#b94a73'),
       pc=c(22,21,22,22,22,22,22,22,22,22,22),bty='n')
  }
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/philr-all-1.png)<!-- -->

```r
adonis(gp.dist~ sample_data(GP)$Treatment)
```

```
## 
## Call:
## adonis(formula = gp.dist ~ sample_data(GP)$Treatment) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##                           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## sample_data(GP)$Treatment  9      6054  672.66  12.282 0.64818  0.001 ***
## Residuals                 60      3286   54.77         0.35182           
## Total                     69      9340                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anosim(gp.dist,sample_data(GP)$Treatment)
```

```
## 
## Call:
## anosim(x = gp.dist, grouping = sample_data(GP)$Treatment) 
## Dissimilarity: euclidean 
## 
## ANOSIM statistic R: 0.5885 
##       Significance: 0.001 
## 
## Permutation: free
## Number of permutations: 999
```

#### Fractions



```r
library(philr)
phyo_frac=subset_samples(phyo2, Comm_Frac=="Frac")
GP <- transform_sample_counts(phyo_frac, function(x) x+1)
#
GP
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 923 taxa and 60 samples ]
## sample_data() Sample Data:       [ 60 samples by 16 sample variables ]
## tax_table()   Taxonomy Table:    [ 923 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 923 tips and 922 internal nodes ]
```

```r
phy_tree(GP) <- makeNodeLabel(phy_tree(GP), method="number", prefix='n')
name.balance(phy_tree(GP), tax_table(GP), 'n1')
```

```
## [1] "Species_/Kingdom_Bacteria"
```

```r
otu.table <- t(otu_table(GP))
tree <- phy_tree(GP)
metadata <- sample_data(GP)
tax <- tax_table(GP)

gp.philr <- philr(otu.table, tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')


gp.dist <- dist(gp.philr, method="euclidean")
gp.pcoa <- ordinate(GP, 'PCoA', distance=gp.dist)
bd3=(sample_data(phyo_frac)$X)
cex_legs=c(4,3.5,3,2.5,2,1.5,1,0.5)


{par(mar=c(5,5,1,8),xpd=T)
  plot(gp.pcoa$vectors[,1],gp.pcoa$vectors[,2],pch=sample_data(GP)$pch,
      bg=as.character(sample_data(GP)$col),cex=bd3,
      xlab = "PCoA1 65.51%", ylab= "PCoA2 21.78%")
  ordiellipse(gp.pcoa$vectors,
              group=as.factor(sample_data(GP)$Treatment),
              label=T,xpd=F,
              kind ='sd',conf=0.8)
  legend(28,15,legend=c("5MC2","6MC2","7MC2","8MC2","9MC2","11MC2","12MC2","14MC2","15MC2","1.819-1.81","1.809-1.8","1.799-17.9","1.789-1.78","1.779-1.77","1.769-1.76","1.759-1.75","1.749-1.73"), pt.bg=c('#b94c3f','#d08f36','#a67e3b','#99af3e','#62a352','#45c097','#697cd4','#9350a1','#b94a73', rep('gray70',8)), pch=21,bty='n', 
         pt.cex = c(1,1,1,1,1,1,1,1,1, 4,3.5,3,2.5,2,1.5,1,0.5))
  }
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

Phylogenetic distance of Fractions with points scaled to the buoyant density. 



```r
adonis(gp.dist~sample_data(phyo_frac)$Treatment)
```

```
## 
## Call:
## adonis(formula = gp.dist ~ sample_data(phyo_frac)$Treatment) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##                                  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
## sample_data(phyo_frac)$Treatment  5    5261.0 1052.21   15.72 0.59276  0.001
## Residuals                        54    3614.4   66.93         0.40724       
## Total                            59    8875.4                 1.00000       
##                                     
## sample_data(phyo_frac)$Treatment ***
## Residuals                           
## Total                               
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anosim(gp.dist,sample_data(phyo_frac)$Treatment)
```

```
## 
## Call:
## anosim(x = gp.dist, grouping = sample_data(phyo_frac)$Treatment) 
## Dissimilarity: euclidean 
## 
## ANOSIM statistic R: 0.5815 
##       Significance: 0.001 
## 
## Permutation: free
## Number of permutations: 999
```

  
### CLR 
Centered log ratio transformation according to [Gloor et al 2017](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full).

I will be removing ASVs that don't have a minimum of 2 reads and represent less than 0.1% of the dataset. I will start with 923 ASVs and widdle that down to 250. 

#### Community PCA


```r
library(CoDaSeq)
library(compositions)
input=t(data.frame(otu_table(phyo2)))
dim(input)
```

```
## [1]  70 923
```

```r
# filter by most abundant taxa 
# must have minimum of 2 reads, and represent 0.1 % of the dataset 
d.subset <- codaSeq.filter(input, samples.by.row=T,min.reads=2,min.prop	=0.001)

dim(d.subset) 
```

```
## [1] 250  70
```

```r
log_rats <- (compositions::clr(t(d.subset)))
OTU=otu_table(t(log_rats),taxa_are_rows = T)
comp_clr=phyloseq(OTU,map,tax2,tree)

# pca on variances 
otu=t(otu_table(comp_clr))
p=prcomp(otu)

{par(mar=c(5,5,1,7), xpd=T)
  plot(p$x[,1],p$x[,2],pch=sample_data(comp_clr)$pch,
     bg=as.character(sample_data(comp_clr)$col),cex=2,
     xlab = "PCoA1 24.17%", ylab= "PCA2 14.21%")
ordiellipse(p$x,
           group=as.factor(sample_data(comp_clr)$Treatment),
           kind ='sd',conf=0.8,xpd=F)
legend(6.5,5,legend=c('Community',"Fraction","5MC2","6MC2","7MC2","8MC2",
                          "9MC2","11MC2","12MC2","14MC2","15MC2", "SidersPond"),
       pt.bg=c('white','white','#b94c3f','#d08f36','#a67e3b','#99af3e','#62a352','#45c097','#697cd4','#9350a1','#b94a73', 'black'),
       col=c('black','black','#b94c3f','#d08f36','#a67e3b','#99af3e','#62a352','#45c097','#697cd4','#9350a1','#b94a73','black'),
       pc=c(22,21,22,22,22,22,22,22,22,22,22,22),bty='n')
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/clr-pcoa-all-1.png)<!-- -->

```r
euc=dist(otu, 'euc')
adonis(euc~sample_data(comp_clr)$Treatment)
```

```
## 
## Call:
## adonis(formula = euc ~ sample_data(comp_clr)$Treatment) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##                                 Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
## sample_data(comp_clr)$Treatment  9    1981.5 220.167  4.8623 0.42174  0.001 ***
## Residuals                       60    2716.8  45.281         0.57826           
## Total                           69    4698.4                 1.00000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anosim(euc,sample_data(comp_clr)$Treatment)
```

```
## 
## Call:
## anosim(x = euc, grouping = sample_data(comp_clr)$Treatment) 
## Dissimilarity: euclidean 
## 
## ANOSIM statistic R: 0.5721 
##       Significance: 0.001 
## 
## Permutation: free
## Number of permutations: 999
```


#### Community w/o Siders


```r
# PCA Aitchinson's without pond 
comp_clr_nopond=subset_samples(comp_clr, Treatment!="SidersPond")
otu=t(otu_table(comp_clr_nopond))

p=prcomp(otu)

{par(mar=c(5,5,1,7),xpd=T)
  plot(p$x[,1],p$x[,2],pch=sample_data(comp_clr_nopond)$pch,
     bg=as.character(sample_data(comp_clr_nopond)$col),cex=2,
     xlab = "PCA1 25.17%", ylab= "PCA2 14.05%")
ordiellipse(p$x,
           group=as.factor(sample_data(comp_clr_nopond)$Treatment),
           kind ='sd',conf=0.8, label=T,xpd=F)
legend(9,5,legend=c('Community',"Fraction","5MC2","6MC2","7MC2","8MC2",
                          "9MC2","11MC2","12MC2","14MC2","15MC2"),
       pt.bg=c('white','white','#b94c3f','#d08f36','#a67e3b','#99af3e','#62a352','#45c097','#697cd4','#9350a1','#b94a73'),
       col=c('black','black','#b94c3f','#d08f36','#a67e3b','#99af3e','#62a352','#45c097','#697cd4','#9350a1','#b94a73'),
       pc=c(22,21,22,22,22,22,22,22,22,22,22),bty='n')
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/clr-nopond-1.png)<!-- -->

```r
euc=dist(otu,'euc')
adonis(euc~sample_data(comp_clr_nopond)$Treatment)
```

```
## 
## Call:
## adonis(formula = euc ~ sample_data(comp_clr_nopond)$Treatment) 
## 
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##                                        Df SumsOfSqs MeanSqs F.Model      R2
## sample_data(comp_clr_nopond)$Treatment  8    1773.4 221.670  4.8954 0.39494
## Residuals                              60    2716.8  45.281         0.60506
## Total                                  68    4490.2                 1.00000
##                                        Pr(>F)    
## sample_data(comp_clr_nopond)$Treatment  0.001 ***
## Residuals                                        
## Total                                            
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
anosim(euc,sample_data(comp_clr_nopond)$Treatment)
```

```
## 
## Call:
## anosim(x = euc, grouping = sample_data(comp_clr_nopond)$Treatment) 
## Dissimilarity: euclidean 
## 
## ANOSIM statistic R: 0.5575 
##       Significance: 0.001 
## 
## Permutation: free
## Number of permutations: 999
```

#### Fractions



```r
comp_clr_frac=subset_samples(comp_clr, Comm_Frac=="Frac")
otu=t(otu_table(comp_clr_frac))
p=prcomp(otu)

bd3=sample_data(comp_clr_frac)$X



{par(mar=c(5,5,1,7),xpd=T)
  plot(p$x[,1],p$x[,2],pch=sample_data(comp_clr_frac)$pch,
     bg=as.character(sample_data(comp_clr_frac)$col),cex=bd3,
     xlab = "PCA1 23.83%", ylab= "PCA2 15.06%")
ordiellipse(p$x,
           group=as.factor(sample_data(comp_clr_frac)$Treatment),
           kind ='sd',conf=0.8, label=T,xpd=F)
  legend(9,6,legend=c("5MC2","6MC2","7MC2","8MC2","9MC2","11MC2","12MC2","14MC2","15MC2","1.819-1.81","1.809-1.8","1.799-17.9","1.789-1.78","1.779-1.77","1.769-1.76","1.759-1.75","1.749-1.73"), pt.bg=c('#b94c3f','#d08f36','#a67e3b','#99af3e','#62a352','#45c097','#697cd4','#9350a1','#b94a73', rep('gray70',8)), pch=21,bty='n', 
         pt.cex = c(1,1,1,1,1,1,1,1,1, 4,3.5,3,2.5,2,1.5,1,0.5))
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

PCA of centered log ratio transformed ASV abundance of fractionated samples. Points are scaled to the buoyant density of the library. 

# Can I find heavy ASVs?

## HR-SIP: ASVs 
HR-SIP = high resolution stable isotope probing


HR-SIP is good because... 

 + Removes sparse OTUs unders the assumption they are sequencing errors 
 + Uses DESeq log2 fold changes to detect enrichment 
 + More sensitive and less true positives than qSIP 

[For more info, check out the tutorial](https://mran.microsoft.com/snapshot/2017-02-04/web/packages/HTSSIP/vignettes/MW_HR_SIP.html)
[and the Paper](https://www.nature.com/articles/ismej2015106)

But first input parameters and get treatments to test against for DESeq.
#### Inputting data
Ok so I had some issues with this. With the l2fc function from DESeq2 you have to make sure that your control is the one being tested.

for example if you run this code:

```r
test=subset_samples(physeq_S2D2_l[[5]], Buoyant_density > 1.76)
d=phyloseq_to_deseq2(test, design=~Treatment)

res = results(diagdds, cooksCutoff = FALSE)
res
alpha = 0.1
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_S2D2_l[[1]])[rownames(sigtab), ], "matrix"))
head(sigtab)
```

The head of the signal tab should look like:

```r
> res
log2 fold change (MLE): Treatment 15MC2 vs 5MC2 
Wald test p-value: Treatment 15MC2 vs 5MC2 
```

Where 15MC2 is has the 13C substrate and 5MC2 has the 12C substrate (aka the control). The key line of code in this is:

```r
sample_data(phyo_frac)$Treatment <- relevel(as.factor(sample_data(phyo_frac)$Treatment) , ref = "5MC2")
```

That is where you set up the reference sample. 


```r
# get just the Fractionated samples 
x<-read.csv(file='/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/asv-table_2.csv',header=TRUE,row.names=1)
OTU = otu_table(x, taxa_are_rows=T)
taxa<-read.csv(file="/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/taxonomy_4nov21.csv",header=TRUE,row.names=1)
t<-as.matrix(taxa)
tax2<-tax_table(t)
map<-import_qiime_sample_data("/Users/ahern/R/trophic_cascades/mc2/26oct21/ASV_run3/map_2.txt")
tree=read.tree("/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/tree.nwk")
phyo = phyloseq(OTU, tax2,map,tree)
phyo1 = subset_taxa(phyo, !Order=="Chloroplast")
phyo2 = subset_taxa(phyo1, !Family=="Mitochondria")

phyo_frac=subset_samples(phyo2, Comm_Frac=="Frac")
input=otu_table(phyo_frac)
d.subset <- codaSeq.filter(input, samples.by.row=F,min.reads=2,min.prop	=0.001)
OTU=otu_table(d.subset,taxa_are_rows = T)
phyo_frac=phyloseq(OTU,map,tax2,tree)

params = get_treatment_params(phyo_frac, c('Treatment', 'TP'))

params = dplyr::filter(params, Treatment!="5MC2")

ex = "(Treatment=='${Treatment}' & TP=='${TP}') | (Treatment=='5MC2' & TP == '-0.45')"
sample_data(phyo_frac)$Treatment <- relevel(as.factor(sample_data(phyo_frac)$Treatment) , ref = "5MC2")
#sample_data(all)$Treatment <- relevel(as.factor(get_variable(all, "Treatment")), ref="5MC2")

physeq_S2D2_l = phyloseq_subset(phyo_frac, params, ex)
physeq_S2D2_l
```

```
## $`(Treatment=='9MC2' & TP=='24.00') | (Treatment=='5MC2' & TP == '-0.45')`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 130 taxa and 23 samples ]
## sample_data() Sample Data:       [ 23 samples by 15 sample variables ]
## tax_table()   Taxonomy Table:    [ 130 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
## 
## $`(Treatment=='11MC2' & TP=='48.00') | (Treatment=='5MC2' & TP == '-0.45')`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 130 taxa and 19 samples ]
## sample_data() Sample Data:       [ 19 samples by 15 sample variables ]
## tax_table()   Taxonomy Table:    [ 130 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
## 
## $`(Treatment=='12MC2' & TP=='72.00') | (Treatment=='5MC2' & TP == '-0.45')`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 130 taxa and 19 samples ]
## sample_data() Sample Data:       [ 19 samples by 15 sample variables ]
## tax_table()   Taxonomy Table:    [ 130 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
## 
## $`(Treatment=='14MC2' & TP=='120.00') | (Treatment=='5MC2' & TP == '-0.45')`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 130 taxa and 20 samples ]
## sample_data() Sample Data:       [ 20 samples by 15 sample variables ]
## tax_table()   Taxonomy Table:    [ 130 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
## 
## $`(Treatment=='15MC2' & TP=='240.00') | (Treatment=='5MC2' & TP == '-0.45')`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 130 taxa and 19 samples ]
## sample_data() Sample Data:       [ 19 samples by 15 sample variables ]
## tax_table()   Taxonomy Table:    [ 130 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 130 tips and 129 internal nodes ]
```
I also filtered the data to only look at the most abundant ASVs to remove potential sequencing errors, rare ASVs, or singletons. So they had to have at least 2 reads/library and represent at least 0.1% of the total dataset. This resulted in 130 ASVs. 

#### Beta Diversity Ordinations
Weighted UniFrac of timepoints vs. control (5MC2)

```r
physeq_S2D2_l_df = SIP_betaDiv_ord(physeq_S2D2_l, parallel=TRUE, method="unifrac", weighted=TRUE,
                                   normalized=FALSE)
```

```
## Run 0 stress 0.03922615 
## Run 1 stress 0.128354 
## Run 2 stress 0.1216414 
## Run 3 stress 0.03922615 
## ... Procrustes: rmse 4.584072e-06  max resid 1.53514e-05 
## ... Similar to previous best
## Run 4 stress 0.03922615 
## ... Procrustes: rmse 9.80909e-06  max resid 4.042597e-05 
## ... Similar to previous best
## Run 5 stress 0.03922615 
## ... New best solution
## ... Procrustes: rmse 2.822336e-06  max resid 9.037165e-06 
## ... Similar to previous best
## Run 6 stress 0.09551262 
## Run 7 stress 0.09551262 
## Run 8 stress 0.1122066 
## Run 9 stress 0.1124714 
## Run 10 stress 0.1278181 
## Run 11 stress 0.1240119 
## Run 12 stress 0.09551262 
## Run 13 stress 0.1122065 
## Run 14 stress 0.1118232 
## Run 15 stress 0.03922615 
## ... Procrustes: rmse 1.85189e-06  max resid 7.222178e-06 
## ... Similar to previous best
## Run 16 stress 0.1122067 
## Run 17 stress 0.03922615 
## ... Procrustes: rmse 1.13369e-06  max resid 2.778324e-06 
## ... Similar to previous best
## Run 18 stress 0.03922615 
## ... Procrustes: rmse 3.745663e-06  max resid 1.32177e-05 
## ... Similar to previous best
## Run 19 stress 0.03922615 
## ... Procrustes: rmse 7.042487e-06  max resid 2.249849e-05 
## ... Similar to previous best
## Run 20 stress 0.03922615 
## ... Procrustes: rmse 5.853971e-06  max resid 2.304946e-05 
## ... Similar to previous best
## *** Solution reached
## Run 0 stress 0.04587765 
## Run 1 stress 0.04587765 
## ... Procrustes: rmse 1.966849e-06  max resid 6.612184e-06 
## ... Similar to previous best
## Run 2 stress 0.04587765 
## ... New best solution
## ... Procrustes: rmse 5.72649e-06  max resid 1.99611e-05 
## ... Similar to previous best
## Run 3 stress 0.04587765 
## ... Procrustes: rmse 2.456937e-06  max resid 6.664494e-06 
## ... Similar to previous best
## Run 4 stress 0.04587765 
## ... New best solution
## ... Procrustes: rmse 2.191324e-06  max resid 7.533483e-06 
## ... Similar to previous best
## Run 5 stress 0.04587765 
## ... Procrustes: rmse 4.344692e-06  max resid 1.505954e-05 
## ... Similar to previous best
## Run 6 stress 0.04587765 
## ... Procrustes: rmse 8.491254e-07  max resid 1.78157e-06 
## ... Similar to previous best
## Run 7 stress 0.04587765 
## ... Procrustes: rmse 3.013982e-06  max resid 8.19701e-06 
## ... Similar to previous best
## Run 8 stress 0.04587765 
## ... Procrustes: rmse 3.042719e-06  max resid 1.041467e-05 
## ... Similar to previous best
## Run 9 stress 0.04587765 
## ... Procrustes: rmse 3.521038e-06  max resid 1.217504e-05 
## ... Similar to previous best
## Run 10 stress 0.04587765 
## ... Procrustes: rmse 5.670067e-06  max resid 1.967128e-05 
## ... Similar to previous best
## Run 11 stress 0.2353787 
## Run 12 stress 0.04587765 
## ... Procrustes: rmse 4.564387e-06  max resid 1.584811e-05 
## ... Similar to previous best
## Run 13 stress 0.04587765 
## ... Procrustes: rmse 2.278527e-06  max resid 7.816363e-06 
## ... Similar to previous best
## Run 14 stress 0.04587765 
## ... Procrustes: rmse 6.190628e-06  max resid 2.152752e-05 
## ... Similar to previous best
## Run 15 stress 0.04587765 
## ... Procrustes: rmse 3.185882e-06  max resid 1.103838e-05 
## ... Similar to previous best
## Run 16 stress 0.04587765 
## ... Procrustes: rmse 6.47797e-07  max resid 1.529028e-06 
## ... Similar to previous best
## Run 17 stress 0.04587765 
## ... Procrustes: rmse 4.674339e-06  max resid 1.628022e-05 
## ... Similar to previous best
## Run 18 stress 0.04587765 
## ... Procrustes: rmse 6.555662e-06  max resid 2.283536e-05 
## ... Similar to previous best
## Run 19 stress 0.04587765 
## ... Procrustes: rmse 8.66799e-07  max resid 2.546354e-06 
## ... Similar to previous best
## Run 20 stress 0.04587765 
## ... Procrustes: rmse 8.537304e-07  max resid 2.441352e-06 
## ... Similar to previous best
## *** Solution reached
## Run 0 stress 0.03459215 
## Run 1 stress 0.1114156 
## Run 2 stress 0.03459207 
## ... New best solution
## ... Procrustes: rmse 9.701884e-05  max resid 0.0003750711 
## ... Similar to previous best
## Run 3 stress 0.06024651 
## Run 4 stress 0.03459208 
## ... Procrustes: rmse 1.387748e-05  max resid 5.404455e-05 
## ... Similar to previous best
## Run 5 stress 0.1111054 
## Run 6 stress 0.03459209 
## ... Procrustes: rmse 4.951614e-05  max resid 0.0001924361 
## ... Similar to previous best
## Run 7 stress 0.03459211 
## ... Procrustes: rmse 5.964867e-05  max resid 0.0002322846 
## ... Similar to previous best
## Run 8 stress 0.1020202 
## Run 9 stress 0.03459212 
## ... Procrustes: rmse 9.021983e-05  max resid 0.0003513389 
## ... Similar to previous best
## Run 10 stress 0.03459214 
## ... Procrustes: rmse 0.0001016413  max resid 0.0003954546 
## ... Similar to previous best
## Run 11 stress 0.03459208 
## ... Procrustes: rmse 2.825838e-05  max resid 0.000109882 
## ... Similar to previous best
## Run 12 stress 0.1084131 
## Run 13 stress 0.03459208 
## ... Procrustes: rmse 4.986914e-06  max resid 1.275368e-05 
## ... Similar to previous best
## Run 14 stress 0.3098689 
## Run 15 stress 0.1020201 
## Run 16 stress 0.03459207 
## ... Procrustes: rmse 3.46868e-06  max resid 1.348562e-05 
## ... Similar to previous best
## Run 17 stress 0.03459208 
## ... Procrustes: rmse 2.331745e-05  max resid 9.081427e-05 
## ... Similar to previous best
## Run 18 stress 0.1063514 
## Run 19 stress 0.06024689 
## Run 20 stress 0.06024725 
## *** Solution reached
## Run 0 stress 0.03673642 
## Run 1 stress 0.1133692 
## Run 2 stress 0.03673642 
## ... New best solution
## ... Procrustes: rmse 2.991359e-06  max resid 9.829231e-06 
## ... Similar to previous best
## Run 3 stress 0.03673642 
## ... New best solution
## ... Procrustes: rmse 1.085277e-05  max resid 3.347209e-05 
## ... Similar to previous best
## Run 4 stress 0.03673642 
## ... New best solution
## ... Procrustes: rmse 1.56139e-06  max resid 4.740351e-06 
## ... Similar to previous best
## Run 5 stress 0.03673642 
## ... New best solution
## ... Procrustes: rmse 4.354168e-06  max resid 1.097473e-05 
## ... Similar to previous best
## Run 6 stress 0.2152178 
## Run 7 stress 0.03673642 
## ... Procrustes: rmse 5.460026e-06  max resid 1.549449e-05 
## ... Similar to previous best
## Run 8 stress 0.03673642 
## ... Procrustes: rmse 1.142662e-05  max resid 3.548929e-05 
## ... Similar to previous best
## Run 9 stress 0.03673642 
## ... New best solution
## ... Procrustes: rmse 3.132212e-06  max resid 5.810786e-06 
## ... Similar to previous best
## Run 10 stress 0.03673642 
## ... Procrustes: rmse 6.785589e-06  max resid 2.106319e-05 
## ... Similar to previous best
## Run 11 stress 0.376911 
## Run 12 stress 0.03673642 
## ... Procrustes: rmse 1.031848e-05  max resid 3.194369e-05 
## ... Similar to previous best
## Run 13 stress 0.03673642 
## ... New best solution
## ... Procrustes: rmse 3.673003e-06  max resid 9.546399e-06 
## ... Similar to previous best
## Run 14 stress 0.03673642 
## ... Procrustes: rmse 8.265277e-06  max resid 2.683435e-05 
## ... Similar to previous best
## Run 15 stress 0.03673642 
## ... Procrustes: rmse 8.343531e-06  max resid 2.847194e-05 
## ... Similar to previous best
## Run 16 stress 0.03673642 
## ... Procrustes: rmse 4.595913e-06  max resid 1.194085e-05 
## ... Similar to previous best
## Run 17 stress 0.03673642 
## ... Procrustes: rmse 5.444032e-06  max resid 1.908349e-05 
## ... Similar to previous best
## Run 18 stress 0.03673642 
## ... New best solution
## ... Procrustes: rmse 1.043443e-06  max resid 2.414741e-06 
## ... Similar to previous best
## Run 19 stress 0.03673642 
## ... Procrustes: rmse 2.860029e-06  max resid 9.089341e-06 
## ... Similar to previous best
## Run 20 stress 0.03673642 
## ... Procrustes: rmse 5.194677e-06  max resid 1.500825e-05 
## ... Similar to previous best
## *** Solution reached
## Run 0 stress 0.01746445 
## Run 1 stress 0.01746446 
## ... Procrustes: rmse 2.03515e-05  max resid 5.019518e-05 
## ... Similar to previous best
## Run 2 stress 0.01746446 
## ... Procrustes: rmse 0.0001147042  max resid 0.0002589563 
## ... Similar to previous best
## Run 3 stress 0.03804692 
## Run 4 stress 0.01746445 
## ... Procrustes: rmse 0.0001046187  max resid 0.000236036 
## ... Similar to previous best
## Run 5 stress 0.03804662 
## Run 6 stress 0.0174645 
## ... Procrustes: rmse 0.0001407141  max resid 0.0002995706 
## ... Similar to previous best
## Run 7 stress 0.03804677 
## Run 8 stress 0.01746453 
## ... Procrustes: rmse 0.0001608986  max resid 0.0003509907 
## ... Similar to previous best
## Run 9 stress 0.01746445 
## ... New best solution
## ... Procrustes: rmse 0.0001069948  max resid 0.0002451488 
## ... Similar to previous best
## Run 10 stress 0.01746445 
## ... Procrustes: rmse 4.043705e-06  max resid 1.028168e-05 
## ... Similar to previous best
## Run 11 stress 0.01746454 
## ... Procrustes: rmse 5.960252e-05  max resid 0.0001365624 
## ... Similar to previous best
## Run 12 stress 0.01746447 
## ... Procrustes: rmse 1.564455e-05  max resid 3.980453e-05 
## ... Similar to previous best
## Run 13 stress 0.01746448 
## ... Procrustes: rmse 2.203186e-05  max resid 5.302009e-05 
## ... Similar to previous best
## Run 14 stress 0.01746449 
## ... Procrustes: rmse 2.827064e-05  max resid 6.287197e-05 
## ... Similar to previous best
## Run 15 stress 0.01746449 
## ... Procrustes: rmse 2.742101e-05  max resid 5.467426e-05 
## ... Similar to previous best
## Run 16 stress 0.01746445 
## ... New best solution
## ... Procrustes: rmse 0.0001004719  max resid 0.0002259958 
## ... Similar to previous best
## Run 17 stress 0.01746443 
## ... New best solution
## ... Procrustes: rmse 7.365844e-05  max resid 0.0001677117 
## ... Similar to previous best
## Run 18 stress 0.01746449 
## ... Procrustes: rmse 5.724408e-05  max resid 0.0001334746 
## ... Similar to previous best
## Run 19 stress 0.01746448 
## ... Procrustes: rmse 5.064653e-05  max resid 0.0001106308 
## ... Similar to previous best
## Run 20 stress 0.0174645 
## ... Procrustes: rmse 6.694249e-05  max resid 0.0001407224 
## ... Similar to previous best
## *** Solution reached
```

```r
physeq_S2D2_l_df %>% .$phyloseq_subset %>% unique
```

```
## [1] (Treatment=='9MC2' & TP=='24.00') | (Treatment=='5MC2' & TP == '-0.45')  
## [2] (Treatment=='11MC2' & TP=='48.00') | (Treatment=='5MC2' & TP == '-0.45') 
## [3] (Treatment=='12MC2' & TP=='72.00') | (Treatment=='5MC2' & TP == '-0.45') 
## [4] (Treatment=='14MC2' & TP=='120.00') | (Treatment=='5MC2' & TP == '-0.45')
## [5] (Treatment=='15MC2' & TP=='240.00') | (Treatment=='5MC2' & TP == '-0.45')
## 5 Levels: (Treatment=='9MC2' & TP=='24.00') | (Treatment=='5MC2' & TP == '-0.45') ...
```

```r
physeq_S2D2_l_df = physeq_S2D2_l_df %>%
  dplyr::mutate(phyloseq_subset = gsub(' \\| ', '\n', phyloseq_subset))
physeq_S2D2_l_df %>% .$phyloseq_subset %>% unique
```

```
## [1] "(Treatment=='9MC2' & TP=='24.00')\n(Treatment=='5MC2' & TP == '-0.45')"  
## [2] "(Treatment=='11MC2' & TP=='48.00')\n(Treatment=='5MC2' & TP == '-0.45')" 
## [3] "(Treatment=='12MC2' & TP=='72.00')\n(Treatment=='5MC2' & TP == '-0.45')" 
## [4] "(Treatment=='14MC2' & TP=='120.00')\n(Treatment=='5MC2' & TP == '-0.45')"
## [5] "(Treatment=='15MC2' & TP=='240.00')\n(Treatment=='5MC2' & TP == '-0.45')"
```

```r
phyloseq_ord_plot(physeq_S2D2_l_df, point_fill = "Treatment")
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


#### Buoyant Density Shifts
Identifying the magnitude of the 13C incoporation. We do this by calculating the percent overlap of beta diversity between treatments (weighted UniFrac). 



```r
wmean = plyr::ldply(physeq_S2D2_l, BD_shift, nperm=10, ex="Treatment=='5MC2'", perm_method='control')

wmean %>% head(n=3)
```

```
##                                                                       .id
## 1 (Treatment=='9MC2' & TP=='24.00') | (Treatment=='5MC2' & TP == '-0.45')
## 2 (Treatment=='9MC2' & TP=='24.00') | (Treatment=='5MC2' & TP == '-0.45')
## 3 (Treatment=='9MC2' & TP=='24.00') | (Treatment=='5MC2' & TP == '-0.45')
##   perm_id   sample.x   distance Replicate.x IS__CONTROL.x BD_min.x BD_max.x
## 1       0 MCP_5MC206 0.07641905           1          TRUE    1.813    1.820
## 2       0 MCP_5MC207 0.11436009           1          TRUE    1.807    1.813
## 3       0 MCP_5MC208 0.11954049           1          TRUE    1.801    1.807
##   BD_range.x perc_overlap n_over_fracs wmean_dist wmean_dist_CI_low_global
## 1      0.007     14.28571            2  0.0856853                0.0106623
## 2      0.006    100.00000            1  0.1143601                0.0106623
## 3      0.006    100.00000            1  0.1195405                0.0106623
##   wmean_dist_CI_high_global wmean_dist_CI_low wmean_dist_CI_high
## 1                 0.0515278       0.046755962         0.07143507
## 2                 0.0515278       0.021058408         0.03270502
## 3                 0.0515278       0.009750047         0.03126197
```

```r
wmean = wmean %>%
  mutate(Treatment = gsub('.+(13C-[A-z]+).+', '\\1', .id))

wmean = wmean %>%
  mutate(.id = gsub(' \\| ', '\n', .id))
wmean %>% .$.id %>% unique
```

```
## [1] "(Treatment=='9MC2' & TP=='24.00')\n(Treatment=='5MC2' & TP == '-0.45')"  
## [2] "(Treatment=='11MC2' & TP=='48.00')\n(Treatment=='5MC2' & TP == '-0.45')" 
## [3] "(Treatment=='12MC2' & TP=='72.00')\n(Treatment=='5MC2' & TP == '-0.45')" 
## [4] "(Treatment=='14MC2' & TP=='120.00')\n(Treatment=='5MC2' & TP == '-0.45')"
## [5] "(Treatment=='15MC2' & TP=='240.00')\n(Treatment=='5MC2' & TP == '-0.45')"
```

```r
# calculating BD shift windows
wmean = wmean %>%
  mutate(BD_shift = wmean_dist > wmean_dist_CI_high) %>%
  arrange(Treatment, BD_min.x) %>%
  group_by(Treatment) %>%
  mutate(window = (BD_shift == TRUE & lag(BD_shift) == TRUE & lag(BD_shift, 2) == TRUE) |
                  (BD_shift == TRUE & lag(BD_shift) == TRUE & lead(BD_shift) == TRUE) |
                  (BD_shift == TRUE & lead(BD_shift) == TRUE & lead(BD_shift, 2) == TRUE),
         BD_shift = BD_shift == TRUE & window == TRUE,
         BD_shift = ifelse(is.na(BD_shift), FALSE, BD_shift)) %>%
  ungroup()
```

Ok lets see how the beta diversity looks in the Treatments vs. the control (5MC2). We can use this to pick buoyant density shift windows. 

```r
y_lab="Beta Diversity"
x_lab="Buoyant Density"
# plotting, with facetting by 13C-treatment
ggplot(wmean, aes(BD_min.x, wmean_dist)) +
  geom_line(alpha=0.7) +
  geom_linerange(aes(ymin=wmean_dist_CI_low,
                     ymax=wmean_dist_CI_high),
                 alpha=0.3) +
  geom_point(aes(color=BD_shift)) +
  labs(x=x_lab, y=y_lab, 
       title='Beta diversity of 5MC2 vs 13-C Treatments') +
  facet_grid( ~ .id) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

Ok since they are like all significant I will just choose the higher fractions: 

 * 1.76-1.81
 * 1.77-1.815
 * 1.78-1.82

Potential significant changes at all timepoints could be due to only one 12-C control. We should have multiple 12-C controls. Potentially have bottles to sample with natural 12-C chemostat inputs for controls to take at all timepoints. 


```r
## focus on two tiempoints first
data<-read.csv("/Users/ahern/R/trophic_cascades/mc2/26oct21/ASV_run3/bd_data.csv",header=T,row.names=1)
MC5=data[1:10,]
MC9=data[11:23,]
MC11=data[24:32,]
MC12=data[33:41,]
MC14=data[42:51,]
MC15=data[52:60,]
library(viridis)
col=viridis(6)

density_windows = data.frame(density_min = c(1.76,1.77,1.78), density_max = c(1.81,1.815,1.82))



{par(mfrow=c(5,1), mar=c(2,2,1,1))
  plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
       bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum")
  rect(xleft=1.76, xright=1.81, ybottom=0, ytop=1, col= rgb(0,.7,1.0,alpha=0.2),
       lty=2, lwd=0.5)
  rect(xleft=1.77, xright=1.815, ybottom=0, ytop=1, col= rgb(0,.8,1.0,alpha=0.2), ,
       lty=2, lwd=0.5)
  rect(xleft=1.78, xright=1.82, ybottom=0, ytop=1, col= rgb(0,1,1.0,alpha=0.2),,
       lty=2, lwd=0.5)
    lines(MC5$BD,MC5$R_MAX,type='o',pch=23,
         bg=col[1])
  lines(MC9$BD,MC9$R_MAX,type='o',pch=21,
        bg=col[2])
  legend('topright', legend=c('5MC2','9MC2'), pt.bg=c(col[1:2]),
         pch=c(23,21),bty='n')
  
  plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
       bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum")
  rect(xleft=1.76, xright=1.81, ybottom=0, ytop=1, col= rgb(0,.7,1.0,alpha=0.2),
       lty=2, lwd=0.5)
  rect(xleft=1.77, xright=1.815, ybottom=0, ytop=1, col= rgb(0,0.8,1.0,alpha=0.2), ,
       lty=2, lwd=0.5)
  rect(xleft=1.78, xright=1.82, ybottom=0, ytop=1, col= rgb(0,1,1.0,alpha=0.2),,
       lty=2, lwd=0.5)
  lines(MC5$BD,MC5$R_MAX,type='o',pch=23,
        bg=col[1])
  lines(MC11$BD,MC11$R_MAX,type='o',pch=21,
        bg=col[3])
  legend('topright', legend=c('5MC2','11MC2'), pt.bg=c(col[1], col[3]),
         pch=c(23,21),bty='n')
  
  plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
       bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum")
  rect(xleft=1.76, xright=1.81, ybottom=0, ytop=1, col= rgb(0,.7,1.0,alpha=0.2),
       lty=2, lwd=0.5)
  rect(xleft=1.77, xright=1.815, ybottom=0, ytop=1, col= rgb(0,0.8,1.0,alpha=0.2), ,
       lty=2, lwd=0.5)
  rect(xleft=1.78, xright=1.82, ybottom=0, ytop=1, col= rgb(0,1,1.0,alpha=0.2),,
       lty=2, lwd=0.5)
  lines(MC5$BD,MC5$R_MAX,type='o',pch=23,
        bg=col[1])
  lines(MC12$BD,MC12$R_MAX,type='o',pch=21,
        bg=col[4])
  legend('topright', legend=c('5MC2','12MC2'), pt.bg=c(col[1], col[4]),
         pch=c(23,21),bty='n')
  
  plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
       bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum")
  rect(xleft=1.76, xright=1.81, ybottom=0, ytop=1, col= rgb(0,.7,1.0,alpha=0.2),
       lty=2, lwd=0.5)
  rect(xleft=1.77, xright=1.815, ybottom=0, ytop=1, col= rgb(0,0.8,1.0,alpha=0.2), ,
       lty=2, lwd=0.5)
  rect(xleft=1.78, xright=1.82, ybottom=0, ytop=1, col= rgb(0,1,1.0,alpha=0.2),,
       lty=2, lwd=0.5)
  lines(MC5$BD,MC5$R_MAX,type='o',pch=23,
        bg=col[1])
  lines(MC14$BD,MC14$R_MAX,type='o',pch=21,
        bg=col[5])
  legend('topright', legend=c('5MC2','14MC2'), pt.bg=c(col[1], col[5]),
         pch=c(23,21),bty='n')
  
  plot(MC5$BD,MC5$R_MAX,type='o', xlim=c(1.731,1.819),pch=23,
       bg=col[1], xlab = "Buoyant Density", ylab = "Ratio of Maximum")
  rect(xleft=1.76, xright=1.81, ybottom=0, ytop=1, col= rgb(0,.7,1.0,alpha=0.2),
       lty=2, lwd=0.5)
  rect(xleft=1.77, xright=1.815, ybottom=0, ytop=1, col= rgb(0,0.8,1.0,alpha=0.2), ,
       lty=2, lwd=0.5)
  rect(xleft=1.78, xright=1.82, ybottom=0, ytop=1, col= rgb(0,1,1.0,alpha=0.2),,
       lty=2, lwd=0.5)
  lines(MC5$BD,MC5$R_MAX,type='o',pch=23,
        bg=col[1])
  lines(MC15$BD,MC15$R_MAX,type='o',pch=21,
        bg=col[6])
  legend('topright', legend=c('5MC2','15MC2'), pt.bg=c(col[1], col[6]),
         pch=c(23,21),bty='n')
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

Illustration of BD windows relative to RNA ratio of maximum


#### The code
HR-SIP code runs each timepoint against the control (5MC2). Used density windows based off of beta diversity measurements (look up one section). This is how the authors describe sparsity thresholds

>The sparsity threshold was set to the value that maximized the number of P-values under a false discovery rate threshold. The specific sparsity threshold was 0.3 meaning that an OTU not found in at least 30% of heavy fractions (control and labeled gradients) in a given day were removed as statistically uninformative.
>
> --- Pepe-Ranney et al 2015

I chose sparsity windows of 0.2, 0.25, and 0.3. Therefore only looking at those ASVs that are at least in 20% of the libraries. 


```r
ncores=1
density_windows = data.frame(density_min = c(1.76,1.77,1.78), density_max = c(1.81,1.815,1.82))
doParallel::registerDoParallel(ncores)
df_l2fc = plyr::ldply(physeq_S2D2_l, 
                      HRSIP, 
                      design = ~Treatment, 
                      density_windows=density_windows,
                      padj_cutoff = 0.1,sparsity_apply = "all",
                      sparsity_threshold = c(0.2,0.25,0.3),
                      .parallel=FALSE)
# save it 
write.csv(df_l2fc, '/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/HRSIP_asvs_19nov21.csv')
```


####  The results

```r
df_l2fc = read.csv('/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/HRSIP_asvs_19nov21.csv')
df_l2fc %>% head(n=3)
```

```
##     X                                                                      .id
## 1  33  (Treatment=='9MC2' & TP=='24.00') | (Treatment=='5MC2' & TP == '-0.45')
## 2  59 (Treatment=='11MC2' & TP=='48.00') | (Treatment=='5MC2' & TP == '-0.45')
## 3 111 (Treatment=='12MC2' & TP=='72.00') | (Treatment=='5MC2' & TP == '-0.45')
##                                OTU log2FoldChange         p padj  Kingdom
## 1 07eb166d091e1aebe641105881f4f2d3     0.08421855 0.7838914    1 Bacteria
## 2 07eb166d091e1aebe641105881f4f2d3    -0.29756958 0.8199338    1 Bacteria
## 3 07eb166d091e1aebe641105881f4f2d3    -2.57469345 0.9994998    1 Bacteria
##            Phyla               Class                 Order         Family
## 1 Proteobacteria Gammaproteobacteria Betaproteobacteriales Rhodocyclaceae
## 2 Proteobacteria Gammaproteobacteria Betaproteobacteriales Rhodocyclaceae
## 3 Proteobacteria Gammaproteobacteria Betaproteobacteriales Rhodocyclaceae
##      Genus Species                           Strain density_min density_max
## 1 Zoogloea         07eb166d091e1aebe641105881f4f2d3        1.77       1.815
## 2 Zoogloea         07eb166d091e1aebe641105881f4f2d3        1.76       1.810
## 3 Zoogloea         07eb166d091e1aebe641105881f4f2d3        1.76       1.810
##   sparsity_threshold sparsity_apply l2fc_threshold    ID timepoint colors
## 1                0.3            all           0.25  9MC2        24 gray70
## 2                0.2            all           0.25 11MC2        48 gray70
## 3                0.2            all           0.25 12MC2        72 gray70
##   colors1 incorp
## 1  gray80  FALSE
## 2  gray80  FALSE
## 3  gray80  FALSE
```

Adjust the ID column 


```r
df_l2fc %>% .$.id %>% unique
```

```
## [1] "(Treatment=='9MC2' & TP=='24.00') | (Treatment=='5MC2' & TP == '-0.45')"  
## [2] "(Treatment=='11MC2' & TP=='48.00') | (Treatment=='5MC2' & TP == '-0.45')" 
## [3] "(Treatment=='12MC2' & TP=='72.00') | (Treatment=='5MC2' & TP == '-0.45')" 
## [4] "(Treatment=='14MC2' & TP=='120.00') | (Treatment=='5MC2' & TP == '-0.45')"
## [5] "(Treatment=='15MC2' & TP=='240.00') | (Treatment=='5MC2' & TP == '-0.45')"
```

```r
df_l2fc = df_l2fc %>%
  mutate(.id = gsub(' \\| ', '\n', .id))
df_l2fc %>% .$.id %>% unique
```

```
## [1] "(Treatment=='9MC2' & TP=='24.00')\n(Treatment=='5MC2' & TP == '-0.45')"  
## [2] "(Treatment=='11MC2' & TP=='48.00')\n(Treatment=='5MC2' & TP == '-0.45')" 
## [3] "(Treatment=='12MC2' & TP=='72.00')\n(Treatment=='5MC2' & TP == '-0.45')" 
## [4] "(Treatment=='14MC2' & TP=='120.00')\n(Treatment=='5MC2' & TP == '-0.45')"
## [5] "(Treatment=='15MC2' & TP=='240.00')\n(Treatment=='5MC2' & TP == '-0.45')"
```

How many incorporators?

```r
# greater than 0.1
padj_cutoff = 0.1
df_l2fc %>% 
  filter(padj < padj_cutoff) %>%
  group_by(.id, sparsity_threshold) %>%
  summarize(n_incorp_OTUs = OTU %>% unique %>% length) %>%
  as.data.frame
```

```
##   n_incorp_OTUs
## 1            26
```

No of ASVs that incorporated 13-C 

| Treatment     | No ASV        |
|---------------|:-------------:|
| 9MC2 24h      | 10            |
| 11MC2 48h     | 11            |
| 12MC2 72h     | 13            |
| 14MC2 120h    | 18            |
| 15MC2 240h    | 18            |

####  HR-SIP ASVs Plots

##### ASVs log2 fold change
Plot log2 fold change (relative to the control 5MC2) for each time timepoint.


```r
t24=subset(df_l2fc, timepoint=="24")
t48=subset(df_l2fc, timepoint=="48")
t72=subset(df_l2fc, timepoint=="72")
t120=subset(df_l2fc, timepoint=="120")
t240=subset(df_l2fc, timepoint=="240")


{
par(mfrow=c(1,5),mar=c(5,1,1,1))
plot(t24$log2FoldChange,1:47,xlab="Log2 Fold Change", ylab=" ",pch=23, col='white',bg='white',main="9MC2 24h",
     yaxt='n')
abline(v=0, col ='black')
segments(0,1:47,t24$log2FoldChange,1:47, col = t24$colors1)
points(t24$log2FoldChange,1:47, pch = 23, bg=t24$colors,cex=1)

plot(t48$log2FoldChange,1:50,xlab="Log2 Fold Change", ylab=" ",pch=23, col='white',bg='white',main="11MC2 48h",
     yaxt='n')
abline(v=0, col ='black')
segments(0,1:50,t48$log2FoldChange,1:50, col = t48$colors1)
points(t48$log2FoldChange,1:50, pch = 23, bg=t48$colors,cex=1)

plot(t72$log2FoldChange,1:50,xlab="Log2 Fold Change", ylab=" ",pch=23, col='white',bg='white',main="12MC2 72h",
     yaxt='n')
abline(v=0, col ='black')
segments(0,1:50,t72$log2FoldChange,1:50, col = t72$colors1)
points(t72$log2FoldChange,1:50, pch = 23, bg=t72$colors,cex=1)

plot(t120$log2FoldChange,1:51,xlab="Log2 Fold Change", ylab=" ",pch=23, col='white',bg='white',main="14MC2 120h",
     yaxt='n')
abline(v=0, col ='black')
segments(0,1:51,t120$log2FoldChange,1:51, col = t120$colors1)
points(t120$log2FoldChange,1:51, pch = 23, bg=t120$colors,cex=1)

plot(t240$log2FoldChange,1:54,xlab="Log2 Fold Change", ylab=" ",pch=23, col='white',bg='white',main="14MC2 120h",
     yaxt='n')
abline(v=0, col ='black')
segments(0,1:54,t240$log2FoldChange,1:54, col = t240$colors1)
points(t240$log2FoldChange,1:54, pch = 23, bg=t240$colors,cex=1)
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/hrsip-allasvs-allresults-1.png)<!-- -->

ASVs that took up 13-C are in black with the red lines from the vertical line. 


##### Relative Abundance ASVs
Get just the ASVs that incorporated 13C

```r
sub=subset(df_l2fc, incorp=="TRUE")
otu=sub$OTU
length(unique(otu))
```

```
## [1] 28
```

```r
# 28 ASVs incorporated 13C

my_subset <- subset(tax_table(phyo_frac), rownames(tax_table(phyo_frac)) %in% otu)
x<-read.csv(file='/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/asv-table.csv',header=TRUE,row.names=1)
OTU = otu_table(x, taxa_are_rows=T)
sums=colSums(OTU)
phyo_hrsip=phyloseq(OTU, map, my_subset, tree)
subset2=otu_table(phyo_hrsip)
sw=sweep(subset2, 2, sums, "/")
sw2=otu_table(sw, taxa_are_rows = T)
phyo_hrsip_trans=phyloseq(sw2, map, tax2, tree)
```

###### Rel Abund Community 
QSip community fraction. Which HR-SIP ASVs are present in the community samples?

```r
comm_hrsip=subset_samples(phyo_hrsip_trans, Comm_Frac=="Community")
comm_tab=as.matrix((otu_table(comm_hrsip)))
{
par(mar=c(7,5,1,1),xpd=T)
barplot(comm_tab, las=2, cex.names = 1, col=colors$x,
        ylab="Relative Abundance Bacterial ASVs",space=0)
box(which='plot')
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/asv-qsip-comm-1.png)<!-- -->

Ok turns out that the most abundant ASVs are in both teh community and fractionated samples. This makes sense because of the sparsity filtering. 

###### Rel Abund Fracs

```r
frac_hrsip=subset_samples(phyo_hrsip_trans, Comm_Frac=="Frac")
frac_tab=as.matrix((otu_table(frac_hrsip)))

space=c(rep("0",10),0.5,rep("0",12),0.5,rep("0",8),0.5,rep("0",8),0.5,rep("0",9),0.5,(rep("0",8)))
{par(mar=c(6,4,3,1),xpd=T)
barplot(frac_tab, las=2, cex.names = 0.7, col=colors$x,
        ylab="Relative Abundance Bacterial ASVs",space=as.numeric(space))
box(which='plot')
text(5,1.05,"5MC2 -0.45h",font=1)
text(17,1.05, "9MC2 24h",font=1)
text(29,1.05,"11MC2 48h",font=1)
text(40, 1.05,'12MC2 72h',font=1)
text(50,1.05, "14MC2 120h",font=1)
text(60,1.05,"15MC2 240h",font=1)
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/asv-qsip-barplot-1.png)<!-- -->

Was there a shift in community and a shift in the fractions? 

###### ASVs in most timepoints 
Five ASVs showed up in all five timepoints:

| ASV        | Class               | Order                 | Family           | Genus        | Species                      |
|------------|:-------------------:|----------------------:|:----------------:|:------------:|:----------------------------:|
| 5ed8f8f3    | Verrucomicrobiae    | Opitutales            | Opitutaceae      | Diplosphaera | uncultured Desulfococcus sp. |
| 801efb1    | Alphaproteobacteria | Rhodospirillales      | Terasakiellaceae |              |                              |
| a6db2d87   | Deltaproteobacteria | Myxococcales          | Nannocystaceae   | Nannocystis  | Nannocystis sp.              |
| acc9603a   | Gammaproteobacteria | Betaproteobacteriales | Rhodocyclaceae   |              |                              |
| cf33813   | Alphaproteobacteria | Rhizobiales           | Rhizobiaceae     | *below*      |                              |

cf33813 Genus = Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium

Two ASVs showed up in four timepoints: 

| ASV        | Class               | Order                 | Family            | Genus              |
|------------|:-------------------:|----------------------:|:-----------------:|:------------------:|
| 3fa2974b   | Alphaproteobacteria | Rhizobiales           | Xanthobacteraceae | Pseudoxanthobacter |
| 83c246     | Gammaproteobacteria | Betaproteobacteriales | Rhodocyclaceae    | Zoogloea           |

ASV 3fa2974b was statistically significant in timepoints 24, 48, 72, 120
ASV 83c246 was statistically significant in timepoints 24, 48, 72, 240

Lets look at the relative abudance of those five ASVs in the dense and light fractions 


```r
phyo_frac2=subset_samples(phyo2, Comm_Frac=="Frac")
phyo_frac_abund=transform_sample_counts(phyo_frac2, function(x) (x/sum(x)))

phyo_frac_abund=subset_samples(phyo_frac_abund, Buoyant_density > 1.759)

five_points=subset_taxa(phyo_frac_abund, Strain=="5ed8f8f3a98e669d81b30c56825e0870" | Strain=="801efb14d6457f799fd2becee6d07b72" | Strain=="a6db2d87f1e7a32216e5041a1920f12e" | Strain == "acc9603a819dc23a606ac6c1b9105b6d" | Strain=="cf3381318906e1e08f58bca905636dc0")
otu=otu_table(five_points)

{
par(mar=c(5,5,1,1),mfcol=c(3,2))
plot(sample_data(five_points)$TP, otu[1,], bg=sample_data(five_points)$col, 
     pch =sample_data(five_points)$pch_sip,cex=2,log='y', xlab="Time (hours)", ylab="Relative Abundance ASV cf33813",
     xaxt='n')
axis(side=1, at =c(0,24,48,72,120,240))
plot(sample_data(five_points)$TP, otu[2,], bg=sample_data(five_points)$col, 
     pch =sample_data(five_points)$pch_sip,cex=2,log='y',xlab="Time (hours)", ylab="Relative Abundance ASV 801efb1",xaxt='n')
axis(side=1, at =c(0,24,48,72,120,240))
plot(sample_data(five_points)$TP, otu[3,], bg=sample_data(five_points)$col, 
     pch =sample_data(five_points)$pch_sip,cex=2,log='y', xlab = "Time (hours)", ylab="Relative Abundance ASV a6db2d87",xaxt='n')
axis(side=1, at =c(0,24,48,72,120,240))
plot(sample_data(five_points)$TP, otu[4,], bg=sample_data(five_points)$col, 
     pch =sample_data(five_points)$pch_sip,cex=2,log='y', xlab = "Time (hours)", ylab="Relative Abundance ASV acc9603a",xaxt='n')
axis(side=1, at =c(0,24,48,72,120,240))
plot(sample_data(five_points)$TP, otu[5,], bg=sample_data(five_points)$col, 
     pch =sample_data(five_points)$pch_sip,cex=2,log='y', xlab = "Time (hours)", ylab="Relative Abundance ASV 5ed8f8f3",xaxt='n')
axis(side=1, at =c(0,24,48,72,120,240))
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/relative-abundance-dense-1.png)<!-- -->


Lets look at the relative abundance of the ASVs that showed up in the most timepoints. 

###### Top ASVs Rel Abund


```r
phyo_frac2=subset_samples(phyo2, Comm_Frac=="Frac")
phyo_frac_abund=transform_sample_counts(phyo_frac2, function(x) (x/sum(x)))
five_four_points=subset_taxa(phyo_frac_abund, Strain=="5ed8f8f3a98e669d81b30c56825e0870" | Strain=="801efb14d6457f799fd2becee6d07b72" | Strain=="a6db2d87f1e7a32216e5041a1920f12e" | Strain == "acc9603a819dc23a606ac6c1b9105b6d" | Strain=="cf3381318906e1e08f58bca905636dc0" | Strain=="3fa2974b1f70109275811beb1ae321c9" | Strain == "83c246c4ed456dea2a15d054f16009b3")
frac_tab=otu_table(five_four_points)

space=c(rep("0",10),0.5,rep("0",12),0.5,rep("0",8),0.5,rep("0",8),0.5,rep("0",9),0.5,(rep("0",8)))
col_asv=c("#ba4450","#b56437","#c99837","#85a444","#50b47b","#7066bc","#b54f90")

{par(mar=c(10,4,3,1),xpd=T)
barplot(frac_tab, las=2,  cex.names = 0.7, col=col_asv,
        ylab="Relative Abundance Bacterial ASVs",space=as.numeric(space))
box(which='plot')
text(5,1.05,"5MC2 -0.45h",font=1)
text(17,1.05, "9MC2 24h",font=1)
text(29,1.05,"11MC2 48h",font=1)
text(40, 1.05,'12MC2 72h',font=1)
text(50,1.05, "14MC2 120h",font=1)
text(60,1.05,"15MC2 240h",font=1)
legend(0,-.35, legend=c("3fa2974b Pseudoxanthobacter", "cf33813 Desulfococcus", "801efb1 Terasakiellaceae", "a6db2d87 Nannocystis", "83c246 Zoogloea", "acc9603a Rhodocyclaceae" , "5ed8f8f3 Rhizobiaceae"), pch =22, pt.bg= col_asv, bty='n',ncol=3, cex=0.8)
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/asv-top5-hrsip-1.png)<!-- -->


```r
read=read.csv(file="/Users/ahern/TC_Exp3/subset_asvs.csv",header=T,row.names = 1)
data=data.frame(read[,2:6], row.names=read$OTU)
phyo_hrsip_f=subset_samples(phyo_hrsip,Comm_Frac=="Frac" )
hrsip_tree=phy_tree(phyo_hrsip_f)
library(viridis)
library(phytools)
c2=c('gray96',(rev(magma(10, begin=0.4))))
{phylo.heatmap(hrsip_tree, data, colors=c2, mar=c(1,5,1,1),lwd=3,
              fsize=c(0.7,1,1),split=c(1,0.3),standardize = F)}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

Two ASVs only had significant values for 13C Glucose assimilation at 9MC2 T24


```r
# d=dcast(sub, OTU ~ timepoint, value.var = "log2FoldChange", fun.aggregate = sum)

phyo_hrsip_f=subset_samples(phyo_hrsip,Comm_Frac=="Frac" )
phyo_hrsip_c=subset_samples(phyo_hrsip,Comm_Frac=="Community" )
f=subset_taxa(phyo_hrsip_f, Strain=="f1e34308d952529287ada31976faec95")
c=subset_taxa(phyo_hrsip_c, Strain=="f1e34308d952529287ada31976faec95")
f1=subset_taxa(phyo_hrsip_f, Strain=="e416d0916760d2fc17b616e2ac3ad855")
c1=subset_taxa(phyo_hrsip_c, Strain=="e416d0916760d2fc17b616e2ac3ad855")
{
par(mar=c(5,5,1,1), mfrow=c(2,2))
plot(sample_data(c)$TP,otu_table(c), bg = sample_data(c)$col, pch =sample_data(c)$pch, cex=2,
     xlab="Time", ylab="Abund Community", main="Abund Comm ASV f1e34308")
plot(sample_data(f)$TP,otu_table(f), bg = sample_data(f)$col, pch =sample_data(f)$pch_sip, 
    cex=(sample_data(f)$X), main="Abund Frac ASV f1e34308", xlab="Time", ylab="Rel Abund Frac")
plot(sample_data(c1)$TP,otu_table(c1), bg = sample_data(c1)$col, pch =sample_data(c1)$pch, cex=2,
     xlab="Time", ylab="Abund Community", main="Abund Comm ASV e416d09")
plot(sample_data(f1)$TP,otu_table(f1), bg = sample_data(f1)$col, pch =sample_data(f1)$pch_sip, 
    cex=(sample_data(f)$X), main="Abund Frac ASV e416d09", xlab="Time", ylab="Abund Frac")
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

One ASV only had significant values at 9MC2 T24 and 11MC2 T48


```r
phyo_hrsip_f=subset_samples(phyo_hrsip,Comm_Frac=="Frac" )
phyo_hrsip_c=subset_samples(phyo_hrsip,Comm_Frac=="Community" )
f=subset_taxa(phyo_hrsip_f, Strain=="efbe1f58b1e2984ddc53a64f047d94ff")
c=subset_taxa(phyo_hrsip_c, Strain=="efbe1f58b1e2984ddc53a64f047d94ff")
{
par(mar=c(5,5,1,1), mfrow=c(1,2))
plot(sample_data(c)$TP,otu_table(c), bg = sample_data(c)$col, pch =sample_data(c)$pch, cex=2,
     xlab="Time", ylab="Abund Community", main="Abund Comm ASV f1e34308", log='y')
plot(sample_data(f)$TP,otu_table(f), bg = sample_data(f)$col, pch =sample_data(f)$pch_sip, 
    cex=(sample_data(f)$X), main="Abund Frac ASV f1e34308", xlab="Time", ylab="Rel Abund Frac",log='y')
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

Three ASVs had significant values from 


```r
phyo_hrsip_f=subset_samples(phyo_hrsip,Comm_Frac=="Frac" )
phyo_hrsip_c=subset_samples(phyo_hrsip,Comm_Frac=="Community" )
f=subset_taxa(phyo_hrsip_f, Strain=="f1e34308d952529287ada31976faec95")
c=subset_taxa(phyo_hrsip_c, Strain=="f1e34308d952529287ada31976faec95")
f1=subset_taxa(phyo_hrsip_f, Strain=="e416d0916760d2fc17b616e2ac3ad855")
c1=subset_taxa(phyo_hrsip_c, Strain=="e416d0916760d2fc17b616e2ac3ad855")
{
par(mar=c(5,5,1,1), mfrow=c(2,2))
plot(sample_data(c)$TP,otu_table(c), bg = sample_data(c)$col, pch =sample_data(c)$pch, cex=2,
     xlab="Time", ylab="Abund Community", main="Abund Comm ASV f1e34308")
plot(sample_data(f)$TP,otu_table(f), bg = sample_data(f)$col, pch =sample_data(f)$pch_sip, 
    cex=(sample_data(f)$X), main="Abund Frac ASV f1e34308", xlab="Time", ylab="Rel Abund Frac")
plot(sample_data(c1)$TP,otu_table(c1), bg = sample_data(c1)$col, pch =sample_data(c1)$pch, cex=2,
     xlab="Time", ylab="Abund Community", main="Abund Comm ASV e416d09")
plot(sample_data(f1)$TP,otu_table(f1), bg = sample_data(f1)$col, pch =sample_data(f1)$pch_sip, 
    cex=(sample_data(f)$X), main="Abund Frac ASV e416d09", xlab="Time", ylab="Abund Frac")
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-26-1.png)<!-- -->



## HR-SIP: Genus  
Ok we did the genus level and now I want to look at the genus level to look at changes. No ASV was found across all timepoints. Though, after thinking this through with Joe and Julie, there are either two explanations to this:

1. Microdiversity perhaps hampering the 13-C signal. 100% ASVs resulting in less incorporators? Microdiversity may over represent the actual diversity in the chemostat. So I want to see if there are any genera that are present at all timepoints. 

2. No ASVs at all timepoints due to the nature of the chemostat. At TP 0, Ashleigh input 13-C labelled glucose at a turnover rate of 1 day^-1. So, only at the end of the 24 hours did the chemostat have glucose that was 100% 13-C enriched. I'm not 100% on this, but a theory is that it will take a bit for the bacteria to see the 13-C glucose, take it up, and then incorporate into its RNA. Therefore, we only see the signal 48 hours later and any ASVs that show up at 24 hours are potentially noise. 
 + Other papers I've seen have only used a 48 hour timepoint for 13-C incorporation, not 24 hour (i.e [this paper](https://www.nature.com/articles/ismej2015106)). Need to do more research. 
 + if this is correct, we should gut out 24 hour RNA sampling from next round of experiments. only take cell counts and delta 13-C


```r
# get just the Fractionated samples 
phyo_frac=subset_samples(phyo2, Comm_Frac=="Frac")
phyo_frac_genus=tax_glom(phyo_frac,taxrank="Genus")
params = get_treatment_params(phyo_frac_genus, c('Treatment', 'TP'))
params = dplyr::filter(params, Treatment!="5MC2")
ex = "(Treatment=='5MC2' & TP=='-0.45') | (Treatment=='${Treatment}' & TP == '${TP}')"
sample_data(phyo_frac_genus)$Treatment <- relevel(as.factor(sample_data(phyo_frac_genus)$Treatment) , ref = "5MC2")
physeq_S2D2_g = phyloseq_subset(phyo_frac_genus, params, ex)
physeq_S2D2_g
```

```
## $`(Treatment=='5MC2' & TP=='-0.45') | (Treatment=='9MC2' & TP == '24.00')`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 361 taxa and 23 samples ]
## sample_data() Sample Data:       [ 23 samples by 15 sample variables ]
## tax_table()   Taxonomy Table:    [ 361 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 361 tips and 360 internal nodes ]
## 
## $`(Treatment=='5MC2' & TP=='-0.45') | (Treatment=='11MC2' & TP == '48.00')`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 361 taxa and 19 samples ]
## sample_data() Sample Data:       [ 19 samples by 15 sample variables ]
## tax_table()   Taxonomy Table:    [ 361 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 361 tips and 360 internal nodes ]
## 
## $`(Treatment=='5MC2' & TP=='-0.45') | (Treatment=='12MC2' & TP == '72.00')`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 361 taxa and 19 samples ]
## sample_data() Sample Data:       [ 19 samples by 15 sample variables ]
## tax_table()   Taxonomy Table:    [ 361 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 361 tips and 360 internal nodes ]
## 
## $`(Treatment=='5MC2' & TP=='-0.45') | (Treatment=='14MC2' & TP == '120.00')`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 361 taxa and 20 samples ]
## sample_data() Sample Data:       [ 20 samples by 15 sample variables ]
## tax_table()   Taxonomy Table:    [ 361 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 361 tips and 360 internal nodes ]
## 
## $`(Treatment=='5MC2' & TP=='-0.45') | (Treatment=='15MC2' & TP == '240.00')`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 361 taxa and 19 samples ]
## sample_data() Sample Data:       [ 19 samples by 15 sample variables ]
## tax_table()   Taxonomy Table:    [ 361 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 361 tips and 360 internal nodes ]
```




### The code

```r
ncores=1
density_windows = data.frame(density_min = c(1.76,1.77,1.78), density_max = c(1.81,1.815,1.82))
doParallel::registerDoParallel(ncores)
df_l2fc = plyr::ldply(physeq_S2D2_g, 
                      HRSIP, 
                      design = ~Treatment, 
                      density_windows=density_windows,
                      padj_cutoff = 0.1,sparsity_apply = "all",
                      sparsity_threshold = c(0.2,0.25,0.3),
                      .parallel=FALSE)
# save it 
write.csv(df_l2fc, '/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/HRSIP_genus.csv')
```

No of Genera that incorporated 13-C 

| Treatment     | No ASV        |
|---------------|:-------------:|
| 9MC2 24h      | 4             |
| 11MC2 48h     | 7             |
| 12MC2 72h     | 9             |
| 14MC2 120h    | 11            |
| 15MC2 240h    | 11            |


### The plots

#### Log2 Fold Change



```r
df_qsip_genus=read.csv(file='/Users/ahern/R/trophic_cascades/mc2/4Nov21/ASVs/HRSIP_genus.csv',header=T)

t24=subset(df_qsip_genus, tp=="24")
t48=subset(df_qsip_genus, tp=="48")
t72=subset(df_qsip_genus, tp=="72")
t120=subset(df_qsip_genus, tp=="120")
t240=subset(df_qsip_genus, tp=="240")

{
par(mfrow=c(1,5),mar=c(5,1,1,1))
plot(t24$log2FoldChange,1:67,xlab="Log2 Fold Change", ylab=" ",pch=23, col='white',bg='white',main="9MC2 24h",
     yaxt='n')
abline(v=0, col ='black')
segments(0,1:67,t24$log2FoldChange,1:67, col = t24$colors1)
points(t24$log2FoldChange,1:67, pch = 23, bg=t24$colors,cex=1)

plot(t48$log2FoldChange,1:65,xlab="Log2 Fold Change", ylab=" ",pch=23, col='white',bg='white',main="11MC2 48h",
     yaxt='n')
abline(v=0, col ='black')
segments(0,1:65,t48$log2FoldChange,1:65, col = t48$colors1)
points(t48$log2FoldChange,1:65, pch = 23, bg=t48$colors,cex=1)

plot(t72$log2FoldChange,1:54,xlab="Log2 Fold Change", ylab=" ",pch=23, col='white',bg='white',main="12MC2 72h",
     yaxt='n')
abline(v=0, col ='black')
segments(0,1:54,t72$log2FoldChange,1:54, col = t72$colors1)
points(t72$log2FoldChange,1:54, pch = 23, bg=t72$colors,cex=1)

plot(t120$log2FoldChange,1:61,xlab="Log2 Fold Change", ylab=" ",pch=23, col='white',bg='white',main="14MC2 120h",
     yaxt='n')
abline(v=0, col ='black')
segments(0,1:61,t120$log2FoldChange,1:61, col = t120$colors1)
points(t120$log2FoldChange,1:61, pch = 23, bg=t120$colors,cex=1)

plot(t240$log2FoldChange,1:51,xlab="Log2 Fold Change", ylab=" ",pch=23, col='white',bg='white',main="14MC2 120h",
     yaxt='n')
abline(v=0, col ='black')
segments(0,1:51,t240$log2FoldChange,1:51, col = t240$colors1)
points(t240$log2FoldChange,1:51, pch = 23, bg=t240$colors,cex=1)
}
```

![](/Users/ahern/Documents/GitHub/Exp3_MC2/docs/index_files/figure-html/unnamed-chunk-29-1.png)<!-- -->


Genera that incorporated 13C in more than 3 time points

| ID        | Class               | Order                 | Family           | Genus         | No TP Incorp 13C |
|------------|:-------------------:|----------------------:|:----------------:|:------------:|:----------------:|
| 5ed8f8f3    | Verrucomicrobiae    | Opitutales            | Opitutaceae      | Diplosphaera | 4                |
| 801efb1    | Alphaproteobacteria | Rhodospirillales      | Terasakiellaceae |              | 4                |
| a6db2d87   | Deltaproteobacteria | Myxococcales          | Nannocystaceae   | Nannocystis  | 5                |
| cf33813   | Alphaproteobacteria | Rhizobiales           | Rhizobiaceae     | *below*      | 4                |

cf33813 Genus = Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium


#### Top Genera Rel Abund


```r
phyo_frac_genus=tax_glom(phyo2,taxrank="Genus")
phyo_frac2=subset_samples(phyo_frac_genus, Comm_Frac=="Frac")
phyo_genus_abund=transform_sample_counts(phyo_frac2, function(x) (x/sum(x)))
five_four_points=subset_taxa(phyo_genus_abund, row.names(tax_table(phyo_genus_abund))=="a6db2d87f1e7a32216e5041a1920f12e" | row.names(tax_table(phyo_genus_abund))=="cf3381318906e1e08f58bca905636dc0" | row.names(tax_table(phyo_genus_abund))=="5ed8f8f3a98e669d81b30c56825e0870" | row.names(tax_table(phyo_genus_abund)) == "801efb14d6457f799fd2becee6d07b72")
frac_tab=otu_table(five_four_points)

space=c(rep("0",10),0.5,rep("0",12),0.5,rep("0",8),0.5,rep("0",8),0.5,rep("0",9),0.5,(rep("0",8)))
col_asv=c("#bb7438","#72ac5c","#7f64b9","#b94b75")

{par(mar=c(15,4,3,1),xpd=T)
barplot(frac_tab, las=2,  cex.names = 0.7, col=col_asv,
        ylab="Relative Abundance Bacterial ASVs",space=as.numeric(space))
box(which='plot')
text(5,1.05,"5MC2 -0.45h",font=1)
text(17,1.05, "9MC2 24h",font=1)
text(29,1.05,"11MC2 48h",font=1)
text(40, 1.05,'12MC2 72h',font=1)
text(50,1.05, "14MC2 120h",font=1)
text(60,1.05,"15MC2 240h",font=1)
legend(0,-.2, legend=c("Alphaproteobacteria Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","Alphaproteobacteria Terasakiellaceae", "Deltaproteobacteria Nannocystis","Verrucomicrobiae Diplosphaera"), pch =22, pt.bg= col_asv, bty='n',ncol=1, cex=0.8)
}
```



Relative abundance in the community 

```r
phyo_g=subset_samples(phyo2, Comm_Frac=="Community")
phyo_comm_g=tax_glom(phyo_g, taxrank = "Genus")
phyo_genus_abund=transform_sample_counts(phyo_comm_g, function(x) (x/sum(x)))
five_four_points=subset_taxa(phyo_genus_abund, row.names(tax_table(phyo_genus_abund))=="a6db2d87f1e7a32216e5041a1920f12e" | row.names(tax_table(phyo_genus_abund))=="cf3381318906e1e08f58bca905636dc0" | row.names(tax_table(phyo_genus_abund))=="5ed8f8f3a98e669d81b30c56825e0870" | row.names(tax_table(phyo_genus_abund)) == "801efb14d6457f799fd2becee6d07b72")
frac_tab=otu_table(five_four_points)
comm_tab=as.matrix((otu_table(five_four_points)))
col_asv=c("#bb7438","#72ac5c","#7f64b9","#b94b75")

{
par(mar=c(20,5,1,1),xpd=T)
barplot(comm_tab, las=2, cex.names = 1, col=col_asv,
        ylab="Relative Abundance Bacterial Genera",space=0)
box(which='plot')
legend(0,-.2, legend=c("Alphaproteobacteria Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","Alphaproteobacteria Terasakiellaceae", "Deltaproteobacteria Nannocystis","Verrucomicrobiae Diplosphaera"), pch =22, pt.bg= col_asv, bty='n',ncol=1, cex=0.8)
}
```


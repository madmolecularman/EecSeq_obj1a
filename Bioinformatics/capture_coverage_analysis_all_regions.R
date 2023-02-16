getwd()
setwd("/home/jgreen/EAGER_OBJ1a/DNA/ddocent")
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)

print(files <- list.files(pattern="C1.hist."))
(files <- c("C1.hist.AllCDS.all.split.txt","C1.hist.AllExon.all.split.txt","C1.hist.AllGene.all.split.txt"))
print(f)

labs <- c("CDS","Exon","Gene","UTR")


cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")

cov <- list()
for (i in 1:length(files)) {
  cov[[i]] <- read.table(files[i])[,c(2,5)]
  cov_cumul=1-cumsum(cov[[i]][,2])
  cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
  cov[[i]]$sample=labs[i]
}

cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7")

p <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(0,100)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.75,0.75))

png(filename="Figure2C1.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p
dev.off()



# Capture 2
print(files <- list.files(pattern="C1.hist."))
(files <- c("C2.hist.AllCDS.all.split.txt","C2.hist.AllExon.all.split.txt","C2.hist.AllGene.all.split.txt","C2.hist.allUTR.all.split.txt"))
print(f)

labs <- c("CDS","Exon","Gene","UTR")


cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")

cov <- list()
for (i in 1:length(files)) {
  cov[[i]] <- read.table(files[i])[,c(2,5)]
  cov_cumul=1-cumsum(cov[[i]][,2])
  cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
  cov[[i]]$sample=labs[i]
}

cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#0072B2")

p <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(0,100)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.75,0.75))

png(filename="Figure2c2.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p
dev.off()

# Capture 3
print(files <- list.files(pattern="C1.hist."))
(files <- c("C3.hist.AllCDS.all.split.txt","C3.hist.AllExon.all.split.txt","C3.hist.AllGene.all.split.txt","C3.hist.allUTR.all.split.txt"))
print(f)

labs <- c("CDS","Exon","Gene","UTR")


cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")

cov <- list()
for (i in 1:length(files)) {
  cov[[i]] <- read.table(files[i])[,c(2,5)]
  cov_cumul=1-cumsum(cov[[i]][,2])
  cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
  cov[[i]]$sample=labs[i]
}

cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#0072B2")

p <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(0,100)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.75,0.75))

png(filename="Figure2c3.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p
dev.off()


# Capture 3
print(files <- list.files(pattern="C1.hist."))
(files <- c("C4.hist.AllCDS.all.split.txt","C4.hist.AllExon.all.split.txt","C4.hist.AllGene.all.split.txt","C4.hist.allUTR.all.split.txt"))
print(f)

labs <- c("CDS","Exon","Gene","UTR")


cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")

cov <- list()
for (i in 1:length(files)) {
  cov[[i]] <- read.table(files[i])[,c(2,5)]
  cov_cumul=1-cumsum(cov[[i]][,2])
  cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
  cov[[i]]$sample=labs[i]
}

cov_df=do.call("rbind",cov)
names(cov_df)[1:2]=c("depth","fraction")

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#0072B2")

p <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(0,100)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.75,0.75))

png(filename="Figure2c4.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p
dev.off()

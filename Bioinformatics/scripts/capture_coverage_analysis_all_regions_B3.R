getwd()
setwd("/home/jgreen/EAGER_OBJ1a/04_coverage_analysis/01_genome_region")
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)

print(files <- list.files(pattern="Capture1_B3.hist"))
(files <- c("Capture1_B3.hist.AllCDS.all.split.txt", "Capture1_B3.hist.AllExon.all.split.txt", "Capture1_B3.hist.AllGene.all.split.txt", "Capture1_B3.hist.AllUTR.all.split.txt"))
print(files)

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

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#56B4E9")

p1 <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(0,25)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ggtitle("Capture1")+
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.50,0.75))

png(filename="Figure2_Capture1_B3.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p1
dev.off()

# Capture 2 B3
print(files <- list.files(pattern="Capture2_B3.hist"))
(files <- c("Capture2_B3.hist.AllCDS.all.split.txt", "Capture2_B3.hist.AllExon.all.split.txt", "Capture2_B3.hist.AllGene.all.split.txt", "Capture2_B3.hist.AllUTR.all.split.txt"))
print(files)

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

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#56B4E9")

p2 <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(0,25)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ggtitle("Capture2")+
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.50,0.75))

png(filename="Figure2_Capture2_B3.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p2
dev.off()

# Capture 3 B3
print(files <- list.files(pattern="Capture3_B3.hist"))
(files <- c("Capture3_B3.hist.AllCDS.all.split.txt", "Capture3_B3.hist.AllExon.all.split.txt", "Capture3_B3.hist.AllGene.all.split.txt", "Capture3_B3.hist.AllUTR.all.split.txt"))
print(files)

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

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#56B4E9")

p3 <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(0,25)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ggtitle("Capture3")+
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.50,0.75))

png(filename="Figure2_Capture3_B3.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p3
dev.off()


# Capture 4 B3
print(files <- list.files(pattern="Capture4_B3.hist"))
(files <- c("Capture4_B3.hist.AllCDS.all.split.txt", "Capture4_B3.hist.AllExon.all.split.txt", "Capture4_B3.hist.AllGene.all.split.txt", "Capture4_B3.hist.AllUTR.all.split.txt"))
print(files)

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

pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#56B4E9")

p4 <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(0,25)+
  scale_alpha(guide = 'none') +
  geom_line(size=1.5)+ 
  #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
  scale_color_manual(values=pcbPalette) +
  scale_fill_manual(values=pcbPalette) +
  ggtitle("Capture4")+
  ylab("% of Bases > Depth")+
  xlab("Depth")+
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme(legend.title = element_blank()) + 
  theme(legend.position=c(0.50,0.75))

png(filename="Figure2_Capture4_B3.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
p4
dev.off()

library(gridExtra)
library(grid)
p5 <- grid.arrange(p1,p2,p3,p4, ncol = 2, nrow = 2,top=textGrob("Capture[1-4]_B3 Reads across genome"))

p5

#png(filename="Figure2_allcapture_B3.png", type="cairo",units="px", width=5600, height=3000, res=600, bg="transparent")
ggsave(file="Figure2_allcapture_B3.png", p5)

dev.off()

## R code for Exon stats

library(MASS)
library(fields)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)

TotalExon <- read.table("./C1.ExonMeanCoverage.txt", header = TRUE)
TotalExon <-as.data.frame(TotalExon)

TotalExon$Exon_Size_Class <-factor(TotalExon$Exon_Size_Class, levels=c("Lower","Middle","Upper"))

TotalExon$Exon_Size_Class <- revalue(TotalExon$Exon_Size_Class, c("Lower"="Lower 10%", "Upper"="Upper 10%", "Middle"="Middle 80%"))

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
TotalExon$density <- get_density(TotalExon$DNA_Coverage)

cbPalette <- c("#009E73","#D55E00","#56B4E9", "#0072B2","#E69F00", "#F0E442" , "#999999","#CC79A7")
t <- cbPalette[1]
cbPalette[1] <- cbPalette[2]
cbPalette[2] <- t

TotalExon$DNA_Coverage <- TotalExon$DNA_Coverage/6
TotalExon$RNA_Coverage <- TotalExon$RNA_Coverage/4


b <- bb<- ggplot(TotalExon, aes(x=RNA_Coverage+1,y=DNA_Coverage+1,alpha = 1/density),fill=Exon_Size_Class,color=Exon_Size_Class) + 
  geom_point(aes(color=TotalExon$Exon_Size_Class,fill=TotalExon$Exon_Size_Class,shape=TotalExon$Exon_Size_Class)) +
  geom_smooth(method="auto", alpha = 0.5, size = 0, se=TRUE)+
  stat_smooth(geom="line", alpha=0.75, size=0.5, linetype="dashed") +
  scale_alpha_continuous(guide = "none",range = c(.2, .95)) + 
  scale_shape_manual(values=c(15,16,17), name="Exon Size Percentile") +
  scale_fill_manual(values=cbPalette , name="Exon Size Percentile")+
  scale_color_manual(values=cbPalette, name="Exon Size Percentile") +
  scale_x_log10(limits=c(1,1300),expand=c(0.02,0), breaks = c(0,1,11,101,1001),
                labels = c("0","0","10","100","1,000"))+
  scale_y_log10(limits=c(1,1800),expand=c(0.02,0), breaks = c(0,1,11,101,1001),
                labels = c("0","0","10","100","1,000"))+
  xlab("Mean RNA Reads per Exon Base Pair")+
  ylab("Mean DNA Reads per Exon Base Pair") +
  theme_bw() +
  theme(legend.position = c(0.9,0.15))  

png(filename="Figure3.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
b
dev.off()    
```
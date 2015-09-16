library(ggplot2)
library(grid)

###################################################################################
# Part I: Mutation Frequency Plot
################################################################################### 
dat1<-read.csv(file="Mutation_Freq.csv",sep='\t',header=FALSE)
colnames(dat1)<-c("BN","MID","Freq")

dat1$BN <- dat1$BN  * 5 

c1 <- ggplot(dat1,aes(x=BN,y=Freq,group=MID))
c1 <- c1 + geom_line()
c1 <- c1 + xlab("Number of mitotic divisions") + ylab("Frequency of mutations")
c1 <- c1 + ggtitle("Frequency of mutations")

###################################################################################
# Part II: Mean mutations per cell plot
###################################################################################
dat1<-read.csv(file="Mean_Mutations.csv",sep='\t',header=FALSE)
colnames(dat1)<-c("BN","Mean_Mut")

dat1$BN <- dat1$BN  * 5 

c2 <- ggplot(dat1,aes(x=BN,y=Mean_Mut))
c2 <- c2 + geom_point()
c2 <- c2 + xlab("Number of mitotic divisions") + ylab("Mean number of mutations per cell")
c2 <- c2 + ggtitle("Mean number of mutations/cell")

###################################################################################
# Part II: Mean mutations per cell plot
###################################################################################
dat1<-read.csv(file="Simulator_Results.csv",sep='\t',header=FALSE)
colnames(dat1)<-c("BN","Mean_CDT","LSC")

dat1$BN <- dat1$BN  * 5 

c3 <- ggplot(dat1,aes(x=BN,y=LSC))
c3 <- c3 + geom_point()
c3 <- c3 + xlab("Number of mitotic divisions") + ylab("Growth rate relative to wild type")
c3 <- c3 + ylim(0,1.0)
c3 <- c3 + ggtitle("Growth rate")

pdf("Simulator_Results.pdf")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,2)))
vplayout<-function(x,y)
viewport(layout.pos.row =x,layout.pos.col = y)
print(c1, vp = vplayout(1,1))
print(c2, vp = vplayout(1,2))
print(c3, vp = vplayout(2,1:2))
dev.off()
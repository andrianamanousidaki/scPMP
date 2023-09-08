library(ggplot2)
library(plyr) 


# Colorblind-friendly palette
cbfly <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Plot Balls:
mydata <- read.csv(file = '../Count_data/Balls_n400_r53.csv',header = FALSE)
true_labels <- read.csv('../Count_data/Balls_n400_r53_Labels.csv',header = FALSE)
mydata$true_labels <- as.factor(true_labels$V1)
p1<-ggplot(mydata,aes(x=V1, y=V2)) 
p1<-p1+geom_point(aes( color=true_labels))+theme(legend.position = "none",
                                                  axis.title.x=element_blank(),
                                                  axis.title.y=element_blank(),
                                                  aspect.ratio=1)
p1<-p1+scale_colour_manual(values=cbfly)

pdf('../Figures/Fig_1_Balls.pdf')
p1
dev.off()

# Plot Elongated with Bridge:
mydata <- read.csv(file = '../Count_data/ElongatedGaussiansWithBridge3.csv', header = FALSE)
true_labels <- read.csv('../Count_data/ElongatedGaussiansWithBridge3_Labels.csv', header = FALSE)
mydata$true_labels <- as.factor(true_labels$V1)
p2<-ggplot(mydata,aes(x=V1, y=V2))
p2<-p2+geom_point(aes( color=true_labels))+theme(legend.position = "none",
                                                  axis.title.x=element_blank(),
                                                  axis.title.y=element_blank(),
                                                  aspect.ratio=1)
p2<- p2 + scale_colour_manual(values=cbfly)

pdf('../Figures/Fig_1_Elongated_with_bridge.pdf')
p2
dev.off()

# Plot Swiss Roll:

mydata <- read.csv(file = '../Count_data/SwissRoll1.csv',header = FALSE)
true_labels <- read.csv('../Count_data/SwissRoll1_Labels.csv',header = FALSE)
mydata$true_labels <- as.factor(true_labels$V1)

p3 <- ggplot(mydata,aes(x=V1, y=V3))

p3 <- p3 + geom_point(aes( color=true_labels))+theme(legend.position = "none",
                                                 axis.title.x=element_blank(),
                                                 axis.title.y=element_blank(),
                                                 aspect.ratio=1)
p3 <- p3 + scale_colour_manual(values=cbfly)


pdf('../Figures/Fig_1_Swiss_roll.pdf')
p3
dev.off()


# Plot GOL manifold:

mydata <- read.csv(file = '../Count_data/GLmanifold_d9_N3000_k3_sig0075.csv',header = FALSE)
true_labels <- read.csv('../Count_data/GLmanifold_d9_N3000_k3_sig0075_Labels.csv',header = FALSE)

pca_res <- prcomp(mydata, center = TRUE,scale. =FALSE, rank.=2)
pcadata <- as.data.frame(pca_res$x)
pcadata$true_labels <- as.factor(true_labels$V1)

p4<-ggplot(pcadata,aes(x=PC1, y=PC2)) 
p4<-p4+geom_point(aes( color=true_labels))+theme(legend.position = "none",
                                                 axis.title.x=element_blank(),
                                                 axis.title.y=element_blank(),
                                                 aspect.ratio=1)
p4<-p4+scale_colour_manual(values=cbfly)


pdf('../Figures/Fig_1_GOL_manifold.pdf')
p4
dev.off()

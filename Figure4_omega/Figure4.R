setwd("/path/too/wd")

#calculate tree-wide bootstrap (BS) averages (each row is a node)
allBS<-read.csv("BSvaltotal.csv")
meanBSs<-colMeans(allBS,na.rm=TRUE)

omeganonroids<-read.csv("omegadonevals.csv")

meanBSswithomegas<-meanBSs[names(meanBSs) %in% omeganonroids$gene]


meanBSswithomegasdf<-as.data.frame(meanBSswithomegas)
rownames(omeganonroids)<-omeganonroids$gene

megedBSandOMEGA<-merge(meanBSswithomegasdf, omeganonroids, by = 0) 


#keep only the genes with mean BS over 70
megedBSandOMEGAover70<-megedBSandOMEGA[which(megedBSandOMEGA$meanBSswithomegas>70),]


#omega values of steroid genes from Codeml output
#K00222incassettes.cml.model0.out:omega (dN/dS) =  0.13044
#K01852incassettes.cml.model0.out:omega (dN/dS) =  0.17216
#K05917incassettes.cml.model0.out:omega (dN/dS) =  0.26342


#write.csv(megedBSandOMEGAover70,'megedBSandOMEGAover70.csv')

#took this file and added the omega values for the steroid genes in CSGs (above) and a column describing whether the gene is one of the threee core CSG genes or not
#called the file megedBSandOMEGAover70brokendown.csv

over70<-read.csv('megedBSandOMEGAover70brokendown.csv')


library(ggplot2)

over70justnonsteroid<-over70[1:344,]
library(grid)
omegaplot<-ggplot(over70, aes(omega)) +
  geom_histogram(data=subset(over70,group=='TM7SF2'), fill="#3C5488FF",binwidth = 0.01) +
  geom_histogram(data=subset(over70,group=='CYP51'), fill="#E64B35FF",binwidth = 0.01) +
  geom_histogram(data=subset(over70,group=='LSS'), fill="#91D1C2FF",binwidth = 0.01) +
    geom_histogram(data=subset(over70,group=='non-steroid'), fill="light gray",binwidth = 0.01) +
  scale_x_continuous(lim = c(0.025, 0.815), breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0),lim = c(0, 40.18))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_density(data=over70justnonsteroid,aes(color = 'density'), color="dark gray",size = 0.75) +
  xlab("Omega") + ylab("Number of genes")+
  theme(text = element_text(size=10))+
  theme(axis.line = element_line(size = 0.25, color = 'black'))+
  theme(axis.ticks = element_line(colour = "black", size = 0.5))+ 
  annotate(geom="text", x=0.12044, y=2, label="ERG24",fontface = 'italic',
           color="black",angle = 45,hjust=0, size=3.5)+
  annotate(geom="text", x=0.16216, y=2, label="ERG7",fontface = 'italic',
           color="black",angle = 45,hjust=0, size=3.5)+
  annotate(geom="text", x=0.25342, y=2, label="CYP51",fontface = 'italic',
           color="black",angle = 45,hjust=0, size=3.5)
omegaplot

ggsave(filename='Fig_4.svg',plot = omegaplot,dpi = 300,width=10, height=10, units="cm", device="svg")






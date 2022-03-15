setwd("/path/to/wd")

library(phytools)
library(ape)
library(geiger)
library(tidyr)
library(dplyr)
library(dbplyr)
library(ggtree)
library(ggplot2)
library(treeio)
library(phangorn)


K00222treehalfcompat <- read.tree("K00222halfcompatfinal.nwk")
K01852treehalfcompat<- read.tree("K01852halfcompatfinal.nwk")
K05917treehalfcompat<- read.tree("K05917halfcompatfinal.nwk")

K00222treehalfcompatmid<-midpoint.root(K00222treehalfcompat)
K01852treehalfcompatmid<-midpoint.root(K01852treehalfcompat)
K05917treehalfcompatmid<-midpoint.root(K05917treehalfcompat)


radialtreeK00222<- ggtree(K00222treehalfcompatmid, layout = "fan")
radialtreeK01852<- ggtree(K01852treehalfcompatmid, layout = "fan")
radialtreeK05917<- ggtree(K05917treehalfcompatmid, layout = "fan")


taxonomyannotations<-read.csv("taxonomyforcasettetrees.csv", row.names = 1, header=F)

sortedtaxonomyfortree<-as.data.frame(taxonomyannotations$V4, drop=F)
rownames(sortedtaxonomyfortree)<-rownames(taxonomyannotations)

tiplabelsfortree<-as.data.frame(taxonomyannotations$V2, drop=F)
rownames(tiplabelsfortree)<-rownames(taxonomyannotations)
colnames(tiplabelsfortree)<-'MAGs'


d = data.frame(label=row.names(taxonomyannotations), label2 = taxonomyannotations$V2)
row.names(d)<-row.names(taxonomyannotations)



p<-radialtreeK00222 %<+% d + geom_tiplab(aes(label=label2), size=1.5, align=T, linetype="dashed", linesize = 0.1, offset=0.25)
p
K00222annotated<-gheatmap(p, sortedtaxonomyfortree, offset=0.45, width=0.14,
         colnames=F) +
  scale_fill_manual(values=c("yellow2",
                             "burlywood4",
                             "lightsalmon2", 
                           "forestgreen",
                             "gray70", 
                             "black",
                             "firebrick3",
                             "lightgreen", 
                             "blue",
                             "orange",
                             "gray80"),
                           name="taxonomy")

svg("./S4_C_K00222annotated.svg",height=10, width=10)
K00222annotated
dev.off()

  
p<-radialtreeK01852 %<+% d + geom_tiplab(aes(label=label2), size=1.5, align=T, linetype="dashed", linesize = 0.1, offset=0.2)
p
K01852annotated<-gheatmap(p, sortedtaxonomyfortree, offset=0.45, width=0.12,
                          colnames=F) +
  scale_fill_manual(values=c("yellow2",
                             "burlywood4",
                             "lightsalmon2", 
                             "forestgreen",
                             "gray70", 
                             "black",
                             "firebrick3",
                             "lightgreen", 
                             "blue",
                             "orange",
                             "gray80"),
                    name="taxonomy")

svg("./S4_B_K01852annotated.svg",height=10, width=10)
K01852annotated
dev.off()





p<-radialtreeK05917 %<+% d + geom_tiplab(aes(label=label2), size=1.5, align=T, linetype="dashed", linesize = 0.1, offset=0.02)
p
K05917annotated<-gheatmap(p, sortedtaxonomyfortree, offset=0.45, width=0.1,
                          colnames=F) +
  scale_fill_manual(values=c("yellow2",
                             "burlywood4",
                             "lightsalmon2", 
                             "forestgreen",
                             "gray70", 
                             "black",
                             "firebrick3",
                             "lightgreen", 
                             "blue",
                             "orange",
                             "gray80"),
                    name="taxonomy")

svg("./S4_A_K05917annotated.svg",height=10, width=10)
K05917annotated
dev.off()





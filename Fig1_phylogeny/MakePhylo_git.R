setwd("/path/to/wd")


library(phytools)
library(ggtree)
library(ggplot2)
library(viridis)
library(ggnewscale)

MAGbacTree<-read.newick("IrciniaAndTaraOceansBacteria.nwk")
TaxonomyFile<-read.csv("TaxonomyFile.csv",row.names = 1, header=T)
MetadataFile<-read.csv("Irciniastats2.csv",row.names = 1, header=T)
steroidannotations<-read.csv("SteroidAnnotationFile.csv",row.names = 1, header=T)

#prepare the metadata for annnotation

sortedtaxonomyfortreephylum<-as.data.frame(TaxonomyFile$phylum2, drop=F)
rownames(sortedtaxonomyfortreephylum)<-rownames(TaxonomyFile)

sortedtaxonomyfortreeclass<-as.data.frame(TaxonomyFile$class, drop=F)
rownames(sortedtaxonomyfortreeclass)<-rownames(TaxonomyFile)

hostdata<-as.data.frame(MetadataFile$Source, drop=F)
rownames(hostdata)<-rownames(MetadataFile)

steroiddata<-as.data.frame(steroidannotations$CSG, drop=F)
rownames(steroiddata)<-rownames(steroidannotations)




radialtree1<- ggtree(MAGbacTree, layout = "fan")+ geom_tiplab(size=0, align=T, linetype="dotted", linesize = 0.1)+theme(text=element_text(size=16,  family="Arial"))

p44<-gheatmap(radialtree1, hostdata, offset=0.0, width=0.05,
              colnames=F)  + scale_fill_manual(values=c("#440154FF", "#482576FF", "#414487FF", "#35608DFF", "#2A788EFF", "#21908CFF", "#22A884FF",
                                                        "#43BF71FF", "#7AD151FF", "#BBDF27FF", "#FDE725FF", "white","gray50"),
                                               name="Host")

p45 <- p44 + new_scale_fill()

p46<-gheatmap(p45, steroiddata, offset=0.1, width=0.05,
              colnames=F)  +scale_fill_manual(values=c("white","orangered2"),
                                              name="CSG")

p47 <- p46 + new_scale_fill()

p48<-gheatmap(p47, sortedtaxonomyfortreephylum, offset=0.2, width=0.05,
         colnames=F) + scale_fill_manual(values=c(rep(c("grey18","grey88"),10)),
                                         name="Taxonomy")


svg("./Figure1.svg",height=20, width=20)
p48
dev.off()



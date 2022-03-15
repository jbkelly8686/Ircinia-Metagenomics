setwd("/path/to/wd")


library(phytools)
library(ape)
library(geiger)


tr<-read.newick("TreeIrciniaTara_git.nwk")
sumstats<-read.table("combinedsumstatshostregions_git.tsv", row.names = 1, header=T)


length(which(tr$tip.label %in% row.names(sumstats))) #should be 913





#let's map the genes onto the trees
library(dplyr)
library(ggplot2)
library(gggenes)
library(ggtree)


Ko02221reordereddircol<-read.csv("genesforplottingontree_git.csv")


numberofgenes1<-Ko02221reordereddircol %>% count(gene)


get_genes <- function(data, genome) {
  filter(data, molecule == genome) %>% pull(gene)
}

g <- unique(Ko02221reordereddircol[,1])
n <- length(g)
d <- matrix(nrow = n, ncol = n)
rownames(d) <- colnames(d) <- g
genes <- lapply(g, get_genes, data = Ko02221reordereddircol)

for (i in 1:n) {
  for (j in 1:i) {
    jaccard_sim <- length(intersect(genes[[i]], genes[[j]])) / 
      length(union(genes[[i]], genes[[j]]))
    d[j, i] <- d[i, j] <- 1 - jaccard_sim
  }
}



tree<-tr

taxabac<-read.table("gtdbtk.bac120.summary_git.tsv", row.names = 1)

#this is required for the clade labels, but some manual editing was performed on Adobe Illustrator

proteoencrows<-c(grep("*Proteobacteria*",taxabac$V2))
proteoencmags<-as.vector(row.names(taxabac)[proteoencrows])
Proteoenc<-findMRCA(tree,proteoencmags, type="node")


Acidoencrows<-c(grep("*Acidobac*",taxabac$V2))
Acidoencmags<-as.vector(row.names(taxabac)[Acidoencrows])
Acidoenc<-findMRCA(tree,Acidoencmags, type="node")

Actinobacteriotarows<-c(grep("*Actinobacteriota*",taxabac$V2))
Actinobacteriotamags<-as.vector(row.names(taxabac)[Actinobacteriotarows])
Actinoenc<-findMRCA(tree,Actinobacteriotamags, type="node")

Bacteroidotaencrows<-c(grep("*Bacteroidota*",taxabac$V2))
Bacteroidotaencmags<-as.vector(row.names(taxabac)[Bacteroidotaencrows])
Bacenc<-findMRCA(tree,Bacteroidotaencmags, type="node")

Binatotaencrows<-c(grep("*Binatota*",taxabac$V2))
Binatotaaencmags<-as.vector(row.names(taxabac)[Binatotaencrows])
Binaenc<-findMRCA(tree,Binatotaaencmags, type="node")


chloroencrows<-c(grep("*Chloroflexota*",taxabac$V2))
chloroencmags<-as.vector(row.names(taxabac)[chloroencrows])
chloroenc<-findMRCA(tree,chloroencmags, type="node")


Cyanobacteriaencrows<-c(grep("*Cyanobacteria*",taxabac$V2))
Cyanobacteriaencmags<-as.vector(row.names(taxabac)[Cyanobacteriaencrows])
CyanoBenc<-findMRCA(tree,Cyanobacteriaencmags, type="node")

Dadaencrows<-c(grep("*Dadabacteria*",taxabac$V2))
Dadaencmags<-as.vector(row.names(taxabac)[Dadaencrows])
Dadaenc<-findMRCA(tree,Dadaencmags, type="node")

Gemmatimonadotacrows<-c(grep("*Gemmatimonadota*",taxabac$V2))
Gemmatimonadotamags<-as.vector(row.names(taxabac)[Gemmatimonadotacrows])
Gemmatimonadotaenc<-findMRCA(tree,Gemmatimonadotamags, type="node")

Latescibacterotarows<-c(grep("*Latescibacterota*",taxabac$V2))
Latescibacterotamags<-as.vector(row.names(taxabac)[Latescibacterotarows])
Latescibacterotaenc<-findMRCA(tree,Latescibacterotamags, type="node")

Nitrospirotarows<-c(grep("*Nitrospirota*",taxabac$V2))
Nitrospirotamags<-as.vector(row.names(taxabac)[Nitrospirotarows])
Nitrospirotaenc<-findMRCA(tree,Nitrospirotamags, type="node")

Poribacteriarows<-c(grep("*Poribacteria*",taxabac$V2))
Poribacteriamags<-as.vector(row.names(taxabac)[Poribacteriarows])
Poribacteriaenc<-findMRCA(tree,Poribacteriamags, type="node")

UBA8248rows<-c(grep("*UBA8248*",taxabac$V2))
UBA8248mags<-as.vector(row.names(taxabac)[UBA8248rows])
UBA8248magsenc<-findMRCA(tree,UBA8248mags, type="node")

Verrucomicrobiotarows<-c(grep("*Verrucomicrobiota*",taxabac$V2))
Verrucomicrobiotamags<-as.vector(row.names(taxabac)[Verrucomicrobiotarows])
Verrucomicrobiotaenc<-findMRCA(tree,Verrucomicrobiotamags, type="node")


Planctomycetotaencrows<-c(grep("*Planctomycetota*",taxabac$V2))
Planctomycetotaencmags<-as.vector(row.names(taxabac)[Planctomycetotaencrows])
Planctomycetotaenc<-findMRCA(tree,Planctomycetotaencmags, type="node")


Bdellovibrionotaencrows<-c(grep("*Bdellovibrionota*",taxabac$V2))
Bdellovibrionotaencmags<-as.vector(row.names(taxabac)[Bdellovibrionotaencrows])
Bdellovibrionotaenc<-findMRCA(tree,Bdellovibrionotaencmags, type="node")


Myxococcotaencrows<-c(grep("*Myxococcota*",taxabac$V2))
Myxococcotaencmags<-as.vector(row.names(taxabac)[Myxococcotaencrows])
Myxococcotaenc<-findMRCA(tree,Myxococcotaencmags, type="node")




UBA8248encrows<-c(grep("*UBA8248*",taxabac$V2))
UBA8248encmags<-as.vector(row.names(taxabac)[UBA8248encrows])
UBA8248enc<-findMRCA(tree,UBA8248encmags, type="node")


Marinisomatotaencrows<-c(grep("*Marinisomatota*",taxabac$V2))
Marinisomatotaencmags<-as.vector(row.names(taxabac)[Marinisomatotaencrows])
Marinisomatotaenc<-findMRCA(tree,Marinisomatotaencmags, type="node")



#now classes
alphaencrows<-c(grep("*Alphaproteo*",taxabac$V2))
alphaencmags<-as.vector(row.names(taxabac)[alphaencrows])
alphaenc<-findMRCA(tree,alphaencmags, type="node")

Gammaencrows<-c(grep("*Gammaproteo*",taxabac$V2))
Gammaencmags<-as.vector(row.names(taxabac)[Gammaencrows])
Gammaenc<-findMRCA(tree,Gammaencmags, type="node")

Cyanobacteriiaencrows<-c(grep("*Cyanobacteriia*",taxabac$V2))
Cyanobacteriiaencmags<-as.vector(row.names(taxabac)[Cyanobacteriiaencrows])
Cyanobacteriiaenc<-findMRCA(tree,Cyanobacteriiaencmags, type="node")



Actinobacteriaencrows<-c(grep("*Actinobacteria*",taxabac$V2))
Actinobacteriaencmags<-as.vector(row.names(taxabac)[Actinobacteriaencrows])
Actinobacteriaenc<-findMRCA(tree,Actinobacteriaencmags, type="node")


Dehalococcoidiaencrows<-c(grep("*Dehalococcoidia*",taxabac$V2))
Dehalococcoidiaencmags<-as.vector(row.names(taxabac)[Dehalococcoidiaencrows])
Dehalococcoidiaenc<-findMRCA(tree,Dehalococcoidiaencmags, type="node")


Anaerolineaeencrows<-c(grep("*Anaerolineae*",taxabac$V2))
Anaerolineaeencmags<-as.vector(row.names(taxabac)[Anaerolineaeencrows])
Anaerolineaeenc<-findMRCA(tree,Anaerolineaeencmags, type="node")


Verrucomicrobiaeencrows<-c(grep("*Verrucomicrobiae*",taxabac$V2))
Verrucomicrobiaeencmags<-as.vector(row.names(taxabac)[Verrucomicrobiaeencrows])
Verrucomicrobiaeenc<-findMRCA(tree,Verrucomicrobiaeencmags, type="node")


Kiritimatiellaeencrows<-c(grep("*Kiritimatiellae*",taxabac$V2))
Kiritimatiellaeencmags<-as.vector(row.names(taxabac)[Kiritimatiellaeencrows])
Kiritimatiellaeenc<-findMRCA(tree,Kiritimatiellaeencmags, type="node")


Planctomycetesencrows<-c(grep("*Planctomycetes*",taxabac$V2))
Planctomycetesencmags<-as.vector(row.names(taxabac)[Planctomycetesencrows])
Planctomycetesenc<-findMRCA(tree,Planctomycetesencmags, type="node")


Phycisphaeraeencrows<-c(grep("*Phycisphaerae*",taxabac$V2))
Phycisphaeraeencmags<-as.vector(row.names(taxabac)[Phycisphaeraeencrows])
Phycisphaeraeenc<-findMRCA(tree,Phycisphaeraeencmags, type="node")


UBA1135encrows<-c(grep("*UBA1135*",taxabac$V2))
UBA1135encmags<-as.vector(row.names(taxabac)[UBA1135encrows])
UBA1135enc<-findMRCA(tree,UBA1135encmags, type="node")


Bacteroidiaencrows<-c(grep("*Bacteroidia*",taxabac$V2))
Bacteroidiaencmags<-as.vector(row.names(taxabac)[Bacteroidiaencrows])
Bacteroidiaenc<-findMRCA(tree,Bacteroidiaencmags, type="node")


Rhodothermiaencrows<-c(grep("*Rhodothermia*",taxabac$V2))
Rhodothermiaencmags<-as.vector(row.names(taxabac)[Rhodothermiaencrows])
Rhodothermiaenc<-findMRCA(tree,Rhodothermiaencmags, type="node")


Marinisomatiaencrows<-c(grep("*Marinisomatia*",taxabac$V2))
Marinisomatiaencmags<-as.vector(row.names(taxabac)[Marinisomatiaencrows])
Marinisomatiaenc<-findMRCA(tree,Marinisomatiaencmags, type="node")


Gemmatimonadetesencrows<-c(grep("*Gemmatimonadetes*",taxabac$V2))
Gemmatimonadetesencmags<-as.vector(row.names(taxabac)[Gemmatimonadetesencrows])
Gemmatimonadetesenc<-findMRCA(tree,Gemmatimonadetesencmags, type="node")


UBA2968encrows<-c(grep("*UBA2968*",taxabac$V2))
UBA2968encmags<-as.vector(row.names(taxabac)[UBA2968encrows])
UBA2968enc<-findMRCA(tree,UBA2968encmags, type="node")


WGA4encrows<-c(grep("*WGA-4*",taxabac$V2))
WGA4encmags<-as.vector(row.names(taxabac)[WGA4encrows])
WGA4enc<-findMRCA(tree,WGA4encmags, type="node")


UBA1144encrows<-c(grep("*UBA1144*",taxabac$V2))
UBA1144encmags<-as.vector(row.names(taxabac)[UBA1144encrows])
UBA1144enc<-findMRCA(tree,UBA1144encmags, type="node")


Nitrospiriaencrows<-c(grep("*Nitrospiria*",taxabac$V2))
Nitrospiriaencmags<-as.vector(row.names(taxabac)[Nitrospiriaencrows])
Nitrospiriaenc<-findMRCA(tree,Nitrospiriaencmags, type="node")


UBA8248encrows<-c(grep("*UBA8248*",taxabac$V2))
UBA8248encmags<-as.vector(row.names(taxabac)[UBA8248encrows])
UBA8248enc<-findMRCA(tree,UBA8248encmags, type="node")


Binatiaencrows<-c(grep("*Binatia*",taxabac$V2))
Binatiaencmags<-as.vector(row.names(taxabac)[Binatiaencrows])
Binatiaenc<-findMRCA(tree,Binatiaencmags, type="node")


Polyangiaencrows<-c(grep("*Polyangia*",taxabac$V2))
Polyangiaencmags<-as.vector(row.names(taxabac)[Polyangiaencrows])
Polyangiaenc<-findMRCA(tree,Polyangiaencmags, type="node")


UBA9160encrows<-c(grep("*UBA9160*",taxabac$V2))
UBA9160encmags<-as.vector(row.names(taxabac)[UBA9160encrows])
UBA9160enc<-findMRCA(tree,UBA9160encmags, type="node")


Vicinamibacteriaencrows<-c(grep("*Vicinamibacteria*",taxabac$V2))
Vicinamibacteriaencmags<-as.vector(row.names(taxabac)[Vicinamibacteriaencrows])
Vicinamibacteriaenc<-findMRCA(tree,Vicinamibacteriaencmags, type="node")


Acidimicrobiiaencrows<-c(grep("*Acidimicrobiia*",taxabac$V2))
Acidimicrobiiaencmags<-as.vector(row.names(taxabac)[Acidimicrobiiaencrows])
Acidimicrobiiaenc<-findMRCA(tree,Acidimicrobiiaencmags, type="node")


UBA2235encrows<-c(grep("*UBA2235*",taxabac$V2))
UBA2235encmags<-as.vector(row.names(taxabac)[UBA2235encrows])
UBA2235enc<-findMRCA(tree,UBA2235encmags, type="node")


UBA11872encrows<-c(grep("*UBA11872*",taxabac$V2))
UBA11872encmags<-as.vector(row.names(taxabac)[UBA11872encrows])
UBA11872enc<-findMRCA(tree,UBA11872encmags, type="node")


Acidobacteriaeencrows<-c(grep("*Acidobacteriae*",taxabac$V2))
Acidobacteriaeencmags<-as.vector(row.names(taxabac)[Acidobacteriaeencrows])
Acidobacteriaeenc<-findMRCA(tree,Acidobacteriaeencmags, type="node")


Bin61encrows<-c(grep("*_c__Bin61_*",taxabac$V2))
Bin61encmags<-as.vector(row.names(taxabac)[Bin61encrows])
Bin61enc<-findMRCA(tree,Bin61encmags, type="node")


Thermoanaerobaculiaencrows<-c(grep("*Thermoanaerobaculia*",taxabac$V2))
Thermoanaerobaculiaencmags<-as.vector(row.names(taxabac)[Thermoanaerobaculiaencrows])
Thermoanaerobaculiaenc<-findMRCA(tree,Thermoanaerobaculiaencmags, type="node")




#Figure 3: focus on Alphaproteobacteria

taxabac<-cbind(taxabac,taxabac$V2)

alphaproteo<-c(grep("*Alphaproteo*",taxabac$V2))
alphataxa<-taxabac[alphaproteo,]

tipstokeepalpha<-c(which(tree$tip.label %in% row.names(alphataxa)==TRUE))
length(tipstokeepalpha)
length(row.names(alphataxa))
alphatree<-keep.tip(tree, tipstokeepalpha)
plot(alphatree)

#get metadata for host annotations
hostdata1<-as.data.frame(sumstats$Source, drop=F)
rownames(hostdata1)<-rownames(sumstats)



 

#this is figure S3


mycolors<-c("white",alpha("#E64B35FF", 0.7),alpha("grey70", 0.7),alpha("yellow", 0.7),alpha("#3C5488FF", 0.7),"blue",alpha("orange",0.7),alpha("#91D1C2FF", 0.7),"pink","white", "white")

svg("FigureS3.svg", height=135,width=20)

p <- ggtree(tree) + 
  geom_tiplab( align=TRUE, linesize=.5) +
  geom_label2(aes(subset=(node==alphaenc),x=branch, label="Alphaproteobacteria"), fill='grey80', size=4.5)+
  geom_label2(aes(subset=(node==Gammaenc),x=branch, label="Gammaproteobacteria"), fill='grey80', size=4.5)+
  geom_label2(aes(subset=(node==Proteoenc),x=branch, label="Proteobacteria"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Gemmatimonadotaenc),x=branch, label="Gemmatimonadota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Latescibacterotaenc),x=branch, label="Latescibacterota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Nitrospirotaenc),x=branch, label="Nitrospirota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Poribacteriaenc),x=branch, label="Poribacteria"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==UBA8248magsenc),x=branch, label="UBA8248"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Acidoenc),x=branch, label="Acidobacteriota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Actinoenc),x=branch, label="Actinobacteriota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Binaenc),x=branch, label="Binatota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==chloroenc),x=branch, label="Chloroflexota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==CyanoBenc),x=branch, label="Cyanobacteria"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Bacenc),x=branch, label="Bacteroidota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Dadaenc),x=branch, label="Dadabacteria"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Verrucomicrobiotaenc),x=branch, label="Verrucomicrobiota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Planctomycetotaenc),x=branch, label="Planctomycetota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Bdellovibrionotaenc),x=branch, label="Bdellovibrionota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Myxococcotaenc),x=branch, label="Myxococcota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==UBA8248enc),x=branch, label="UBA8248"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Marinisomatotaenc),x=branch, label="Marinisomatota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Cyanobacteriiaenc),x=branch, label="Cyanobacteriia"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Actinobacteriaenc),x=branch, label="Actinobacteria"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Dehalococcoidiaenc),x=branch, label="Dehalococcoidia"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Anaerolineaeenc),x=branch, label="Anaerolineae"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Verrucomicrobiaeenc),x=branch, label="Verrucomicrobiae"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Kiritimatiellaeenc),x=branch, label="Kiritimatiellae"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Planctomycetesenc),x=branch, label="Planctomycetes"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Phycisphaeraeenc),x=branch, label="Phycisphaerae"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==UBA1135enc),x=branch, label="UBA1135"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Bacteroidiaenc),x=branch, label="Bacteroidia"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Rhodothermiaenc),x=branch, label="Rhodothermia"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Marinisomatiaenc),x=branch, label="Marinisomatia"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Gemmatimonadetesenc),x=branch, label="Gemmatimonadetes"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==UBA2968enc),x=branch, label="UBA2968"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==WGA4enc),x=branch, label="WGA-4"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==UBA1144enc),x=branch, label="UBA1144"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Nitrospiriaenc),x=branch, label="Nitrospiria"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==UBA8248enc),x=branch, label="UBA8248"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Binatiaenc),x=branch, label="Binatia"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Polyangiaenc),x=branch, label="Polyangia"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==UBA9160enc),x=branch, label="UBA9160"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Acidimicrobiiaenc),x=branch, label="Acidimicrobiia"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Vicinamibacteriaenc),x=branch, label="Vicinamibacteria"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==UBA2235enc),x=branch, label="UBA2235"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==UBA11872enc),x=branch, label="UBA11872"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Acidobacteriaeenc),x=branch, label="Acidobacteriae"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Bin61enc),x=branch, label="Bin61"), fill='grey80', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Thermoanaerobaculiaenc),x=branch, label="Thermoanaerobaculia"), fill='grey80', size=4.5,nudge_x = 0.075) +
  geom_treescale(y=-1)+ 
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene, forward = direction),
             data = Ko02221reordereddircol, geom = geom_motif, panel = 'Alignment',
             on='arb', align = 'centre' ) +
  scale_fill_manual(values = mycolors)

pB <- p + new_scale_fill()

pA<-gheatmap(pB, hostdata1, offset=1, width=0.15,
             colnames=F)  + scale_fill_manual(values=c("#440154FF","#482576FF", "#414487FF", "#35608DFF", "#2A788EFF", "#21908CFF", "#22A884FF",
                                                       "#43BF71FF", "#7AD151FF", "#BBDF27FF","red", "#FDE725FF", "orange","indianred2", "gray80","white","gray50"),
                                              name="Host")
facet_widths(pA, widths=c(1,1)) 
dev.off()




#and this is figure 3

mycolorsalpha<-c("white",alpha("#E64B35FF", 0.7),alpha("grey70", 0.7),alpha("yellow", 0.7),alpha("#3C5488FF", 0.7),alpha("orange",0.7),alpha("#91D1C2FF", 0.7),alpha("white",0.7),"white")

svg("Figure3.svg", height=30,width=20)
p<-ggtree(alphatree) + geom_treescale(y=-1)+ geom_tiplab(align=T)+ #xlim_tree(5.5) +
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene, forward = direction),
             data = Ko02221reordereddircol, geom = geom_motif, panel = 'Alignment',
             on='arb', align = 'centre' ) +
  scale_fill_manual(values = mycolorsalpha)

pB <- p + new_scale_fill()

pA<-gheatmap(pB, hostdata1, offset=1, width=0.15,
             colnames=F)  + scale_fill_manual(values=c("#482576FF", "#414487FF", "#35608DFF", "#2A788EFF", "#21908CFF", "#22A884FF",
                                                       "#43BF71FF", "#7AD151FF", "#BBDF27FF","red", "#FDE725FF", "orange","indianred2", "gray80","white","gray50"),
                                              name="Host")

facet_widths(pA, widths=c(1,1)) 
dev.off()






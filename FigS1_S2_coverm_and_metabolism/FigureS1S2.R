setwd("/path/to/wd")


library(phytools)
library(ape)
library(geiger)
library(factoextra)
library(vegan)
library(ggtree)
library(tidyverse)
library(ggtree)
library(ggimage)
library(ggplot2)
library(evobiR)
library(picante)
library(tibble)
library(ggfortify)

KOdata <- read.table("ko_frequency_table_git12.tsv", header = TRUE, row.names = 1)
KOdatatnotpresabs<-as.data.frame(t(KOdata))
presabsKO<-as.matrix((KOdatatnotpresabs > 0) + 0)

#need to keep just the Ircinia MAGs for plotting on the tree
nonIrciniaMags<-row.names(presabsKO)[grep("*GCA*",row.names(presabsKO))]
todrop<-row.names(presabsKO)%in%nonIrciniaMags
todrop1<-which(as.numeric(todrop)==1)

presabsKOircinia<-presabsKO[-todrop1,]

WoodLjundahlpathway<-c(which(colnames(presabsKOircinia)=="K00198"),
                       which(colnames(presabsKOircinia)=="K05299"),
                       which(colnames(presabsKOircinia)=="K15022"),
                       which(colnames(presabsKOircinia)=="K01938"),
                       which(colnames(presabsKOircinia)=="K01491"),
                       which(colnames(presabsKOircinia)=="K00297"),
                       which(colnames(presabsKOircinia)=="K15023"),
                       which(colnames(presabsKOircinia)=="K14138"),
                       which(colnames(presabsKOircinia)=="K00197"),
                       which(colnames(presabsKOircinia)=="K00194"))

dicarboxylate4hydroxybutyratecycle<-c(which(colnames(presabsKOircinia)=="K00169"),
                                      which(colnames(presabsKOircinia)=="K00170"),
                                      which(colnames(presabsKOircinia)=="K00171"),
                                      which(colnames(presabsKOircinia)=="K00172"),
                                      which(colnames(presabsKOircinia)=="K01007"),
                                      which(colnames(presabsKOircinia)=="K01595"),
                                      which(colnames(presabsKOircinia)=="K00024"),
                                      which(colnames(presabsKOircinia)=="K01676"),
                                      which(colnames(presabsKOircinia)=="K01677"),
                                      which(colnames(presabsKOircinia)=="K01678"),
                                      which(colnames(presabsKOircinia)=="K00239"),
                                      which(colnames(presabsKOircinia)=="K00240"),
                                      which(colnames(presabsKOircinia)=="K00241"),
                                      which(colnames(presabsKOircinia)=="K18860"),
                                      which(colnames(presabsKOircinia)=="K01902"),
                                      which(colnames(presabsKOircinia)=="K01903"),
                                      which(colnames(presabsKOircinia)=="K15038"),
                                      which(colnames(presabsKOircinia)=="K15017"),
                                      which(colnames(presabsKOircinia)=="K14465"),
                                      which(colnames(presabsKOircinia)=="K14467"),
                                      which(colnames(presabsKOircinia)=="K18861"),
                                      which(colnames(presabsKOircinia)=="K14534"),
                                      which(colnames(presabsKOircinia)=="K15016"),
                                      which(colnames(presabsKOircinia)=="K00626"))


threehydroxypropionatebicycle<-c(which(colnames(presabsKOircinia)=="K01961"),
                             which(colnames(presabsKOircinia)=="K01962"),
                             which(colnames(presabsKOircinia)=="K01963"),
                             which(colnames(presabsKOircinia)=="K02160"),
                             which(colnames(presabsKOircinia)=="K14468"),
                             which(colnames(presabsKOircinia)=="K14469"),
                             which(colnames(presabsKOircinia)=="K15052"),
                             which(colnames(presabsKOircinia)=="K05606"),
                             which(colnames(presabsKOircinia)=="K01847"),
                             which(colnames(presabsKOircinia)=="K01848"),
                             which(colnames(presabsKOircinia)=="K01849"),
                             which(colnames(presabsKOircinia)=="K14471"),
                             which(colnames(presabsKOircinia)=="K14472"),
                             which(colnames(presabsKOircinia)=="K00239"),
                             which(colnames(presabsKOircinia)=="K00240"),
                             which(colnames(presabsKOircinia)=="K00241"),
                             which(colnames(presabsKOircinia)=="K01679"),
                             which(colnames(presabsKOircinia)=="K08691"),
                             which(colnames(presabsKOircinia)=="K14449"),
                             which(colnames(presabsKOircinia)=="K14470"),
                             which(colnames(presabsKOircinia)=="K09709"))

threehydroxypropionate4hydroxybutyrate<-c(which(colnames(presabsKOircinia)=="K01964"),
                                      which(colnames(presabsKOircinia)=="K15036"),
                                      which(colnames(presabsKOircinia)=="K15037"),
                                      which(colnames(presabsKOircinia)=="K15017"),
                                      which(colnames(presabsKOircinia)=="K15039"),
                                      which(colnames(presabsKOircinia)=="K15018"),
                                      which(colnames(presabsKOircinia)=="K15019"),
                                      which(colnames(presabsKOircinia)=="K15020"),
                                      which(colnames(presabsKOircinia)=="K05606"),
                                      which(colnames(presabsKOircinia)=="K01848"),
                                      which(colnames(presabsKOircinia)=="K01849"),
                                      which(colnames(presabsKOircinia)=="K15038"),
                                      which(colnames(presabsKOircinia)=="K14465"),
                                      which(colnames(presabsKOircinia)=="K14466"),
                                      which(colnames(presabsKOircinia)=="K18861"),
                                      which(colnames(presabsKOircinia)=="K14534"),
                                      which(colnames(presabsKOircinia)=="K15016"),
                                      which(colnames(presabsKOircinia)=="K00626"))

reductivecitricacidcycle<-c(which(colnames(presabsKOircinia)=="K00169"),
                            which(colnames(presabsKOircinia)=="K00170"),
                            which(colnames(presabsKOircinia)=="K00171"),
                            which(colnames(presabsKOircinia)=="K00172"),
                            which(colnames(presabsKOircinia)=="K03737"),
                            which(colnames(presabsKOircinia)=="K01007"),
                            which(colnames(presabsKOircinia)=="K01006"),
                            which(colnames(presabsKOircinia)=="K01595"),
                            which(colnames(presabsKOircinia)=="K01959"),
                            which(colnames(presabsKOircinia)=="K01960"),
                            which(colnames(presabsKOircinia)=="K01958"),
                            which(colnames(presabsKOircinia)=="K00024"),
                            which(colnames(presabsKOircinia)=="K01676"),
                            which(colnames(presabsKOircinia)=="K01679"),
                            which(colnames(presabsKOircinia)=="K01677"),
                            which(colnames(presabsKOircinia)=="K01678"),
                            which(colnames(presabsKOircinia)=="K00239"),
                            which(colnames(presabsKOircinia)=="K00240"),
                            which(colnames(presabsKOircinia)=="K00241"),
                            which(colnames(presabsKOircinia)=="K00242"),
                            which(colnames(presabsKOircinia)=="K00244"),
                            which(colnames(presabsKOircinia)=="K00245"),
                            which(colnames(presabsKOircinia)=="K00246"),
                            which(colnames(presabsKOircinia)=="K00247"),
                            which(colnames(presabsKOircinia)=="K18556"),
                            which(colnames(presabsKOircinia)=="K18557"),
                            which(colnames(presabsKOircinia)=="K18558"),
                            which(colnames(presabsKOircinia)=="K18559"),
                            which(colnames(presabsKOircinia)=="K18560"),
                            which(colnames(presabsKOircinia)=="K01902"),
                            which(colnames(presabsKOircinia)=="K01903"),
                            which(colnames(presabsKOircinia)=="K00174"),
                            which(colnames(presabsKOircinia)=="K00175"),
                            which(colnames(presabsKOircinia)=="K00177"),
                            which(colnames(presabsKOircinia)=="K00176"),
                            which(colnames(presabsKOircinia)=="K00031"),
                            which(colnames(presabsKOircinia)=="K01681"),
                            which(colnames(presabsKOircinia)=="K01682"),
                            which(colnames(presabsKOircinia)=="K15230"),
                            which(colnames(presabsKOircinia)=="K15231"),
                            which(colnames(presabsKOircinia)=="K15232"),
                            which(colnames(presabsKOircinia)=="K15233"),
                            which(colnames(presabsKOircinia)=="K15234"))

CalvinBensonBasshamcycle<-c(which(colnames(presabsKOircinia)=="K00855"),
                            which(colnames(presabsKOircinia)=="K01601"),
                            which(colnames(presabsKOircinia)=="K01602"),
                            which(colnames(presabsKOircinia)=="K00927"),
                            which(colnames(presabsKOircinia)=="K05298"),
                            which(colnames(presabsKOircinia)=="K00150"),
                            which(colnames(presabsKOircinia)=="K00134"),
                            which(colnames(presabsKOircinia)=="K01623"),
                            which(colnames(presabsKOircinia)=="K01624"),
                            which(colnames(presabsKOircinia)=="K03841"),
                            which(colnames(presabsKOircinia)=="K02446"),
                            which(colnames(presabsKOircinia)=="K11532"),
                            which(colnames(presabsKOircinia)=="K01086"),
                            which(colnames(presabsKOircinia)=="K04041"),
                            which(colnames(presabsKOircinia)=="K00615"),
                            which(colnames(presabsKOircinia)=="K01100"),
                            which(colnames(presabsKOircinia)=="K01807"),
                            which(colnames(presabsKOircinia)=="K01808"))


Nitrogenmetabolism<-c(which(colnames(presabsKOircinia)=="K01455"),
                      which(colnames(presabsKOircinia)=="K02575"),
                      which(colnames(presabsKOircinia)=="K15576"),
                      which(colnames(presabsKOircinia)=="K15577"),
                      which(colnames(presabsKOircinia)=="K15578"),
                      which(colnames(presabsKOircinia)=="K15579"),
                      which(colnames(presabsKOircinia)=="K00367"),
                      which(colnames(presabsKOircinia)=="K10534"),
                      which(colnames(presabsKOircinia)=="K00370"),
                      which(colnames(presabsKOircinia)=="K00371"),
                      which(colnames(presabsKOircinia)=="K00374"),
                      which(colnames(presabsKOircinia)=="K02567"),
                      which(colnames(presabsKOircinia)=="K02568"),
                      which(colnames(presabsKOircinia)=="K00372"),
                      which(colnames(presabsKOircinia)=="K00360"),
                      which(colnames(presabsKOircinia)=="K17877"),
                      which(colnames(presabsKOircinia)=="K00362"),
                      which(colnames(presabsKOircinia)=="K00363"),
                      which(colnames(presabsKOircinia)=="K00366"),
                      which(colnames(presabsKOircinia)=="K03385"),
                      which(colnames(presabsKOircinia)=="K15876"),
                      which(colnames(presabsKOircinia)=="K00368"),
                      which(colnames(presabsKOircinia)=="K15864"),
                      which(colnames(presabsKOircinia)=="K04561"),
                      which(colnames(presabsKOircinia)=="K02305"),
                      which(colnames(presabsKOircinia)=="K15877"),
                      which(colnames(presabsKOircinia)=="K00376"),
                      which(colnames(presabsKOircinia)=="K02586"),
                      which(colnames(presabsKOircinia)=="K02591"),
                      which(colnames(presabsKOircinia)=="K02588"),
                      which(colnames(presabsKOircinia)=="K00531"),
                      which(colnames(presabsKOircinia)=="K22896"),
                      which(colnames(presabsKOircinia)=="K22897"),
                      which(colnames(presabsKOircinia)=="K22898"),
                      which(colnames(presabsKOircinia)=="K22899"),
                      which(colnames(presabsKOircinia)=="K20932"),
                      which(colnames(presabsKOircinia)=="K20933"),
                      which(colnames(presabsKOircinia)=="K20934"),
                      which(colnames(presabsKOircinia)=="K20935"),
                      which(colnames(presabsKOircinia)=="K10944"),
                      which(colnames(presabsKOircinia)=="K10945"),
                      which(colnames(presabsKOircinia)=="K10946"),
                      which(colnames(presabsKOircinia)=="K05601"),
                      which(colnames(presabsKOircinia)=="K10535"),
                      which(colnames(presabsKOircinia)=="K00459"),
                      which(colnames(presabsKOircinia)=="K19823"),
                      which(colnames(presabsKOircinia)=="K01501"),
                      which(colnames(presabsKOircinia)=="K15371"),
                      which(colnames(presabsKOircinia)=="K00260"),
                      which(colnames(presabsKOircinia)=="K00261"),
                      which(colnames(presabsKOircinia)=="K00262"),
                      which(colnames(presabsKOircinia)=="K01915"),
                      which(colnames(presabsKOircinia)=="K00264"),
                      which(colnames(presabsKOircinia)=="K00265"),
                      which(colnames(presabsKOircinia)=="K00266"),
                      which(colnames(presabsKOircinia)=="K00284"),
                      which(colnames(presabsKOircinia)=="K01948"),
                      which(colnames(presabsKOircinia)=="K01725"),
                      which(colnames(presabsKOircinia)=="K00926"),
                      which(colnames(presabsKOircinia)=="K01672"),
                      which(colnames(presabsKOircinia)=="K18245"),
                      which(colnames(presabsKOircinia)=="K18246"),
                      which(colnames(presabsKOircinia)=="K01673"),
                      which(colnames(presabsKOircinia)=="K01674"),
                      which(colnames(presabsKOircinia)=="K01743"))


Sulfurmetabolism<-c(which(colnames(presabsKOircinia)=="K02048"),
                    which(colnames(presabsKOircinia)=="K23163"),
                    which(colnames(presabsKOircinia)=="K02046"),
                    which(colnames(presabsKOircinia)=="K02047"),
                    which(colnames(presabsKOircinia)=="K02045"),
                    which(colnames(presabsKOircinia)=="K15551"),
                    which(colnames(presabsKOircinia)=="K15552"),
                    which(colnames(presabsKOircinia)=="K10831"),
                    which(colnames(presabsKOircinia)=="K03119"),
                    which(colnames(presabsKOircinia)=="K15553"),
                    which(colnames(presabsKOircinia)=="K15554"),
                    which(colnames(presabsKOircinia)=="K15555"),
                    which(colnames(presabsKOircinia)=="K04091"),
                    which(colnames(presabsKOircinia)=="K00299"),
                    which(colnames(presabsKOircinia)=="K13811"),
                    which(colnames(presabsKOircinia)=="K00955"),
                    which(colnames(presabsKOircinia)=="K00956"),
                    which(colnames(presabsKOircinia)=="K00957"),
                    which(colnames(presabsKOircinia)=="K00958"),
                    which(colnames(presabsKOircinia)=="K00988"),
                    which(colnames(presabsKOircinia)=="K22966"),
                    which(colnames(presabsKOircinia)=="K00860"),
                    which(colnames(presabsKOircinia)=="K01082"),
                    which(colnames(presabsKOircinia)=="K15759"),
                    which(colnames(presabsKOircinia)=="K15422"),
                    which(colnames(presabsKOircinia)=="K06881"),
                    which(colnames(presabsKOircinia)=="K00394"),
                    which(colnames(presabsKOircinia)=="K00395"),
                    which(colnames(presabsKOircinia)=="K05907"),
                    which(colnames(presabsKOircinia)=="K00390"),
                    which(colnames(presabsKOircinia)=="K00387"),
                    which(colnames(presabsKOircinia)=="K05301"),
                    which(colnames(presabsKOircinia)=="K00386"),
                    which(colnames(presabsKOircinia)=="K21307"),
                    which(colnames(presabsKOircinia)=="K21308"),
                    which(colnames(presabsKOircinia)=="K21309"),
                    which(colnames(presabsKOircinia)=="K17222"),
                    which(colnames(presabsKOircinia)=="K17223"),
                    which(colnames(presabsKOircinia)=="K17226"),
                    which(colnames(presabsKOircinia)=="K17227"),
                    which(colnames(presabsKOircinia)=="K17224"),
                    which(colnames(presabsKOircinia)=="K17225"),
                    which(colnames(presabsKOircinia)=="K22622"),
                    which(colnames(presabsKOircinia)=="K11180"),
                    which(colnames(presabsKOircinia)=="K11181"),
                    which(colnames(presabsKOircinia)=="K00380"),
                    which(colnames(presabsKOircinia)=="K00381"),
                    which(colnames(presabsKOircinia)=="K16950"),
                    which(colnames(presabsKOircinia)=="K16951"),
                    which(colnames(presabsKOircinia)=="K00385"),
                    which(colnames(presabsKOircinia)=="K00392"),
                    which(colnames(presabsKOircinia)=="K17218"),
                    which(colnames(presabsKOircinia)=="K17995"),
                    which(colnames(presabsKOircinia)=="K17996"),
                    which(colnames(presabsKOircinia)=="K17993"),
                    which(colnames(presabsKOircinia)=="K17994"),
                    which(colnames(presabsKOircinia)=="K17229"),
                    which(colnames(presabsKOircinia)=="K17230"),
                    which(colnames(presabsKOircinia)=="K16952"),
                    which(colnames(presabsKOircinia)=="K17219"),
                    which(colnames(presabsKOircinia)=="K17220"),
                    which(colnames(presabsKOircinia)=="K17221"),
                    which(colnames(presabsKOircinia)=="K22470"),
                    which(colnames(presabsKOircinia)=="K17725"),
                    which(colnames(presabsKOircinia)=="K16936"),
                    which(colnames(presabsKOircinia)=="K16937"),
                    which(colnames(presabsKOircinia)=="K05908"),
                    which(colnames(presabsKOircinia)=="K08357"),
                    which(colnames(presabsKOircinia)=="K08358"),
                    which(colnames(presabsKOircinia)=="K08359"),
                    which(colnames(presabsKOircinia)=="K08352"),
                    which(colnames(presabsKOircinia)=="K08353"),
                    which(colnames(presabsKOircinia)=="K08354"),
                    which(colnames(presabsKOircinia)=="K01011"),
                    which(colnames(presabsKOircinia)=="K02439"),
                    which(colnames(presabsKOircinia)=="K00640"),
                    which(colnames(presabsKOircinia)=="K23304"),
                    which(colnames(presabsKOircinia)=="K10150"),
                    which(colnames(presabsKOircinia)=="K01738"),
                    which(colnames(presabsKOircinia)=="K13034"),
                    which(colnames(presabsKOircinia)=="K17069"),
                    which(colnames(presabsKOircinia)=="K00651"),
                    which(colnames(presabsKOircinia)=="K00641"),
                    which(colnames(presabsKOircinia)=="K01739"),
                    which(colnames(presabsKOircinia)=="K10764"),
                    which(colnames(presabsKOircinia)=="K17217"),
                    which(colnames(presabsKOircinia)=="K17285"),
                    which(colnames(presabsKOircinia)=="K17228"),
                    which(colnames(presabsKOircinia)=="K16968"),
                    which(colnames(presabsKOircinia)=="K16969"),
                    which(colnames(presabsKOircinia)=="K15762"),
                    which(colnames(presabsKOircinia)=="K15765"),
                    which(colnames(presabsKOircinia)=="K07306"),
                    which(colnames(presabsKOircinia)=="K07307"),
                    which(colnames(presabsKOircinia)=="K07308"),
                    which(colnames(presabsKOircinia)=="K00184"),
                    which(colnames(presabsKOircinia)=="K00185"),
                    which(colnames(presabsKOircinia)=="K16964"),
                    which(colnames(presabsKOircinia)=="K16965"),
                    which(colnames(presabsKOircinia)=="K16966"),
                    which(colnames(presabsKOircinia)=="K16953"),
                    which(colnames(presabsKOircinia)=="K17486"),
                    which(colnames(presabsKOircinia)=="K20034"),
                    which(colnames(presabsKOircinia)=="K20035"),
                    which(colnames(presabsKOircinia)=="K20036"),
                    which(colnames(presabsKOircinia)=="K16967"),
                    which(colnames(presabsKOircinia)=="K21310"),
                    which(colnames(presabsKOircinia)=="K16954"),
                    which(colnames(presabsKOircinia)=="K16955"))


Methanemetabolism<-c(which(colnames(presabsKOircinia)=="K16157"),
                     which(colnames(presabsKOircinia)=="K16158"),
                     which(colnames(presabsKOircinia)=="K16159"),
                     which(colnames(presabsKOircinia)=="K16160"),
                     which(colnames(presabsKOircinia)=="K16161"),
                     which(colnames(presabsKOircinia)=="K16162"),
                     which(colnames(presabsKOircinia)=="K10944"),
                     which(colnames(presabsKOircinia)=="K10945"),
                     which(colnames(presabsKOircinia)=="K10946"),
                     which(colnames(presabsKOircinia)=="K14028"),
                     which(colnames(presabsKOircinia)=="K16254"),
                     which(colnames(presabsKOircinia)=="K16255"),
                     which(colnames(presabsKOircinia)=="K14029"),
                     which(colnames(presabsKOircinia)=="K16256"),
                     which(colnames(presabsKOircinia)=="K16257"),
                     which(colnames(presabsKOircinia)=="K16258"),
                     which(colnames(presabsKOircinia)=="K16259"),
                     which(colnames(presabsKOircinia)=="K16260"),
                     which(colnames(presabsKOircinia)=="K23995"),
                     which(colnames(presabsKOircinia)=="K17066"),
                     which(colnames(presabsKOircinia)=="K00093"),
                     which(colnames(presabsKOircinia)=="K00148"),
                     which(colnames(presabsKOircinia)=="K17067"),
                     which(colnames(presabsKOircinia)=="K17068"),
                     which(colnames(presabsKOircinia)=="K03396"),
                     which(colnames(presabsKOircinia)=="K00121"),
                     which(colnames(presabsKOircinia)=="K01070"),
                     which(colnames(presabsKOircinia)=="K00122"),
                     which(colnames(presabsKOircinia)=="K00123"),
                     which(colnames(presabsKOircinia)=="K22515"),
                     which(colnames(presabsKOircinia)=="K00124"),
                     which(colnames(presabsKOircinia)=="K00127"),
                     which(colnames(presabsKOircinia)=="K00126"),
                     which(colnames(presabsKOircinia)=="K22516"),
                     which(colnames(presabsKOircinia)=="K00125"),
                     which(colnames(presabsKOircinia)=="K05299"),
                     which(colnames(presabsKOircinia)=="K15022"),
                     which(colnames(presabsKOircinia)=="K00192"),
                     which(colnames(presabsKOircinia)=="K00195"),
                     which(colnames(presabsKOircinia)=="K00193"),
                     which(colnames(presabsKOircinia)=="K00197"),
                     which(colnames(presabsKOircinia)=="K00194"),
                     which(colnames(presabsKOircinia)=="K00198"),
                     which(colnames(presabsKOircinia)=="K00196"),
                     which(colnames(presabsKOircinia)=="K00600"),
                     which(colnames(presabsKOircinia)=="K00830"),
                     which(colnames(presabsKOircinia)=="K00018"),
                     which(colnames(presabsKOircinia)=="K11529"),
                     which(colnames(presabsKOircinia)=="K01689"),
                     which(colnames(presabsKOircinia)=="K01595"),
                     which(colnames(presabsKOircinia)=="K00024"),
                     which(colnames(presabsKOircinia)=="K08692"),
                     which(colnames(presabsKOircinia)=="K14067"),
                     which(colnames(presabsKOircinia)=="K08691"),
                     which(colnames(presabsKOircinia)=="K17100"),
                     which(colnames(presabsKOircinia)=="K00863"),
                     which(colnames(presabsKOircinia)=="K01623"),
                     which(colnames(presabsKOircinia)=="K11645"),
                     which(colnames(presabsKOircinia)=="K01624"),
                     which(colnames(presabsKOircinia)=="K01622"),
                     which(colnames(presabsKOircinia)=="K16305"),
                     which(colnames(presabsKOircinia)=="K16306"),
                     which(colnames(presabsKOircinia)=="K03841"),
                     which(colnames(presabsKOircinia)=="K02446"),
                     which(colnames(presabsKOircinia)=="K11532"),
                     which(colnames(presabsKOircinia)=="K01086"),
                     which(colnames(presabsKOircinia)=="K04041"),
                     which(colnames(presabsKOircinia)=="K00850"),
                     which(colnames(presabsKOircinia)=="K16370"),
                     which(colnames(presabsKOircinia)=="K21071"),
                     which(colnames(presabsKOircinia)=="K00918"),
                     which(colnames(presabsKOircinia)=="K08094"),
                     which(colnames(presabsKOircinia)=="K08093"),
                     which(colnames(presabsKOircinia)=="K13812"),
                     which(colnames(presabsKOircinia)=="K13831"),
                     which(colnames(presabsKOircinia)=="K00317"),
                     which(colnames(presabsKOircinia)=="K18277"),
                     which(colnames(presabsKOircinia)=="K07811"),
                     which(colnames(presabsKOircinia)=="K03532"),
                     which(colnames(presabsKOircinia)=="K03533"),
                     which(colnames(presabsKOircinia)=="K07821"),
                     which(colnames(presabsKOircinia)=="K07812"),
                     which(colnames(presabsKOircinia)=="K15228"),
                     which(colnames(presabsKOircinia)=="K15229"),
                     which(colnames(presabsKOircinia)=="K08685"),
                     which(colnames(presabsKOircinia)=="K22081"),
                     which(colnames(presabsKOircinia)=="K22082"),
                     which(colnames(presabsKOircinia)=="K22083"),
                     which(colnames(presabsKOircinia)=="K22084"),
                     which(colnames(presabsKOircinia)=="K22085"),
                     which(colnames(presabsKOircinia)=="K22086"),
                     which(colnames(presabsKOircinia)=="K22087"),
                     which(colnames(presabsKOircinia)=="K00200"),
                     which(colnames(presabsKOircinia)=="K00201"),
                     which(colnames(presabsKOircinia)=="K00202"),
                     which(colnames(presabsKOircinia)=="K00203"),
                     which(colnames(presabsKOircinia)=="K00204"),
                     which(colnames(presabsKOircinia)=="K00205"),
                     which(colnames(presabsKOircinia)=="K11260"),
                     which(colnames(presabsKOircinia)=="K11261"),
                     which(colnames(presabsKOircinia)=="K00672"),
                     which(colnames(presabsKOircinia)=="K01499"),
                     which(colnames(presabsKOircinia)=="K00319"),
                     which(colnames(presabsKOircinia)=="K00440"),
                     which(colnames(presabsKOircinia)=="K00441"),
                     which(colnames(presabsKOircinia)=="K00442"),
                     which(colnames(presabsKOircinia)=="K00443"),
                     which(colnames(presabsKOircinia)=="K00300"),
                     which(colnames(presabsKOircinia)=="K10714"),
                     which(colnames(presabsKOircinia)=="K13942"),
                     which(colnames(presabsKOircinia)=="K10713"),
                     which(colnames(presabsKOircinia)=="K00320"),
                     which(colnames(presabsKOircinia)=="K00577"),
                     which(colnames(presabsKOircinia)=="K00578"),
                     which(colnames(presabsKOircinia)=="K00579"),
                     which(colnames(presabsKOircinia)=="K00580"),
                     which(colnames(presabsKOircinia)=="K00581"),
                     which(colnames(presabsKOircinia)=="K00582"),
                     which(colnames(presabsKOircinia)=="K00583"),
                     which(colnames(presabsKOircinia)=="K00584"),
                     which(colnames(presabsKOircinia)=="K00399"),
                     which(colnames(presabsKOircinia)=="K00400"),
                     which(colnames(presabsKOircinia)=="K00401"),
                     which(colnames(presabsKOircinia)=="K00402"),
                     which(colnames(presabsKOircinia)=="K03421"),
                     which(colnames(presabsKOircinia)=="K03422"),
                     which(colnames(presabsKOircinia)=="K22480"),
                     which(colnames(presabsKOircinia)=="K22481"),
                     which(colnames(presabsKOircinia)=="K22482"),
                     which(colnames(presabsKOircinia)=="K03388"),
                     which(colnames(presabsKOircinia)=="K03389"),
                     which(colnames(presabsKOircinia)=="K03390"),
                     which(colnames(presabsKOircinia)=="K08264"),
                     which(colnames(presabsKOircinia)=="K08265"),
                     which(colnames(presabsKOircinia)=="K14127"),
                     which(colnames(presabsKOircinia)=="K14126"),
                     which(colnames(presabsKOircinia)=="K14128"),
                     which(colnames(presabsKOircinia)=="K00925"),
                     which(colnames(presabsKOircinia)=="K00625"),
                     which(colnames(presabsKOircinia)=="K13788"),
                     which(colnames(presabsKOircinia)=="K01895"),
                     which(colnames(presabsKOircinia)=="K00169"),
                     which(colnames(presabsKOircinia)=="K00170"),
                     which(colnames(presabsKOircinia)=="K00172"),
                     which(colnames(presabsKOircinia)=="K00171"),
                     which(colnames(presabsKOircinia)=="K01007"),
                     which(colnames(presabsKOircinia)=="K01834"),
                     which(colnames(presabsKOircinia)=="K15633"),
                     which(colnames(presabsKOircinia)=="K15634"),
                     which(colnames(presabsKOircinia)=="K15635"),
                     which(colnames(presabsKOircinia)=="K00058"),
                     which(colnames(presabsKOircinia)=="K00831"),
                     which(colnames(presabsKOircinia)=="K01079"),
                     which(colnames(presabsKOircinia)=="K02203"),
                     which(colnames(presabsKOircinia)=="K22305"),
                     which(colnames(presabsKOircinia)=="K14080"),
                     which(colnames(presabsKOircinia)=="K04480"),
                     which(colnames(presabsKOircinia)=="K14081"),
                     which(colnames(presabsKOircinia)=="K14082"),
                     which(colnames(presabsKOircinia)=="K14083"),
                     which(colnames(presabsKOircinia)=="K14084"),
                     which(colnames(presabsKOircinia)=="K16176"),
                     which(colnames(presabsKOircinia)=="K16177"),
                     which(colnames(presabsKOircinia)=="K16178"),
                     which(colnames(presabsKOircinia)=="K16179"),
                     which(colnames(presabsKOircinia)=="K08097"),
                     which(colnames(presabsKOircinia)=="K05979"),
                     which(colnames(presabsKOircinia)=="K05884"),
                     which(colnames(presabsKOircinia)=="K06034"),
                     which(colnames(presabsKOircinia)=="K13039"),
                     which(colnames(presabsKOircinia)=="K11779"),
                     which(colnames(presabsKOircinia)=="K11781"),
                     which(colnames(presabsKOircinia)=="K11780"),
                     which(colnames(presabsKOircinia)=="K14941"),
                     which(colnames(presabsKOircinia)=="K11212"),
                     which(colnames(presabsKOircinia)=="K12234"),
                     which(colnames(presabsKOircinia)=="K14940"),
                     which(colnames(presabsKOircinia)=="K10977"),
                     which(colnames(presabsKOircinia)=="K16792"),
                     which(colnames(presabsKOircinia)=="K16793"),
                     which(colnames(presabsKOircinia)=="K10978"),
                     which(colnames(presabsKOircinia)=="K18933"),
                     which(colnames(presabsKOircinia)=="K06914"),
                     which(colnames(presabsKOircinia)=="K09733"),
                     which(colnames(presabsKOircinia)=="K19793"),
                     which(colnames(presabsKOircinia)=="K07144"),
                     which(colnames(presabsKOircinia)=="K07072"))

orderofgenes<-c(WoodLjundahlpathway,dicarboxylate4hydroxybutyratecycle,threehydroxypropionatebicycle,threehydroxypropionate4hydroxybutyrate,reductivecitricacidcycle,CalvinBensonBasshamcycle,Nitrogenmetabolism,Sulfurmetabolism,Methanemetabolism)


presabsKOirciniareordered<-presabsKOircinia[,c(WoodLjundahlpathway,dicarboxylate4hydroxybutyratecycle,threehydroxypropionatebicycle,threehydroxypropionate4hydroxybutyrate,reductivecitricacidcycle,CalvinBensonBasshamcycle,Nitrogenmetabolism,Sulfurmetabolism,Methanemetabolism)]
#write.csv(presabsKOirciniareordered,"presabsKOirciniareordered_git.csv")


colnames(presabsKOirciniareordered)<-c("K00198",
                                                 "K05299",
                                                 "K15022",
                                                 "K01938",
                                                 "K01491",
                                                 "K00297",
                                                 "K15023",
                                                 "K14138",
                                                 "K00197",
                                                 "K00194",
                                                 "K00169",
                                                 "K00170",
                                                 "K00171",
                                                 "K00172",
                                                 "K01007",
                                                 "K01595",
                                                 "K00024",
                                                 "K01676",
                                                 "K01677",
                                                 "K01678",
                                                 "K00239",
                                                 "K00240",
                                                 "K00241",
                                                 "K18860",
                                                 "K01902",
                                                 "K01903",
                                                 " K15038",
                                                 " K15017",
                                                 " K14465",
                                                 "K14467",
                                                 " K18861",
                                                 " K14534",
                                                 " K15016",
                                                 " K00626",
                                                 "K01961",
                                                 "K01962",
                                                 "K01963",
                                                 "K02160",
                                                 "K14468",
                                                 "K14469",
                                                 "K15052",
                                                 " K05606",
                                                 "K01847",
                                                 " K01848",
                                                 " K01849",
                                                 "K14471",
                                                 "K14472",
                                                 " K00239",
                                                 " K00240",
                                                 " K00241",
                                                 " K01679",
                                                 " K08691",
                                                 "K14449",
                                                 "K14470",
                                                 "K09709",
                                                 "K01964",
                                                 "K15036",
                                                 "K15037",
                                                 "K15017",
                                                 "K15039",
                                                 "K15018",
                                                 "K15019",
                                                 "K15020",
                                                 "K05606",
                                                 "K01848",
                                                 "K01849",
                                                 "K15038",
                                                 "K14465",
                                                 "K14466",
                                                 "K18861",
                                                 "K14534",
                                                 "K15016",
                                                 "K00626",
                                                 " K00169",
                                                 " K00170",
                                                 " K00171",
                                                 " K00172",
                                                 "K03737",
                                                 " K01007",
                                                 "K01006",
                                                 " K01595",
                                                 "K01959",
                                                 "K01960",
                                                 "K01958",
                                                 " K00024",
                                                 " K01676",
                                                 "K01679",
                                                 " K01677",
                                                 " K01678",
                                                 "  K00239",
                                                 "  K00240",
                                                 "  K00241",
                                                 "K00242",
                                                 "K00244",
                                                 "K00245",
                                                 "K00246",
                                                 "K00247",
                                                 "K18556",
                                                 "K18557",
                                                 "K18558",
                                                 "K18559",
                                                 "K18560",
                                                 " K01902",
                                                 " K01903",
                                                 "K00174",
                                                 "K00175",
                                                 "K00177",
                                                 "K00176",
                                                 "K00031",
                                                 "K01681",
                                                 "K01682",
                                                 "K15230",
                                                 "K15231",
                                                 "K15232",
                                                 "K15233",
                                                 "K15234",
                                                 "K00855",
                                                 "K01601",
                                                 "K01602",
                                                 "K00927",
                                                 "K05298",
                                                 "K00150",
                                                 "K00134",
                                                 " K01623",
                                                 " K01624",
                                                 " K03841",
                                                 " K02446",
                                                 " K11532",
                                                 " K01086",
                                                 " K04041",
                                                 "K00615",
                                                 "K01100",
                                                 "K01807",
                                                 "K01808",
                                                 "K01455",
                                                 "K02575",
                                                 "K15576",
                                                 "K15577",
                                                 "K15578",
                                                 "K15579",
                                                 "K00367",
                                                 "K10534",
                                                 "K00370",
                                                 "K00371",
                                                 "K00374",
                                                 "K02567",
                                                 "K02568",
                                                 "K00372",
                                                 "K00360",
                                                 "K17877",
                                                 "K00362",
                                                 "K00363",
                                                 "K00366",
                                                 "K03385",
                                                 "K15876",
                                                 "K00368",
                                                 "K15864",
                                                 "K04561",
                                                 "K02305",
                                                 "K15877",
                                                 "K00376",
                                                 "K02586",
                                                 "K02591",
                                                 "K02588",
                                                 "K00531",
                                                 "K20932",
                                                 "K20933",
                                                 "K20934",
                                                 "K20935",
                                                 " K10944",
                                                 " K10945",
                                                 " K10946",
                                                 "K05601",
                                                 "K10535",
                                                 "K00459",
                                                 "K19823",
                                                 "K01501",
                                                 "K15371",
                                                 "K00260",
                                                 "K00261",
                                                 "K00262",
                                                 "K01915",
                                                 "K00264",
                                                 "K00265",
                                                 "K00266",
                                                 "K00284",
                                                 "K01948",
                                                 "K01725",
                                                 "K00926",
                                                 "K01672",
                                                 "K18245",
                                                 "K18246",
                                                 "K01673",
                                                 "K01674",
                                                 "K02048",
                                                 "K02046",
                                                 "K02047",
                                                 "K02045",
                                                 "K15551",
                                                 "K15552",
                                                 "K10831",
                                                 "K03119",
                                                 "K15553",
                                                 "K15554",
                                                 "K15555",
                                                 "K04091",
                                                 "K00299",
                                                 "K13811",
                                                 "K00955",
                                                 "K00956",
                                                 "K00957",
                                                 "K00958",
                                                 "K00988",
                                                 "K00860",
                                                 "K01082",
                                                 "K15759",
                                                 "K15422",
                                                 "K06881",
                                                 "K00394",
                                                 "K00395",
                                                 "K05907",
                                                 "K00390",
                                                 "K00387",
                                                 "K05301",
                                                 "K17222",
                                                 "K17223",
                                                 "K17226",
                                                 "K17227",
                                                 "K17224",
                                                 "K17225",
                                                 "K11180",
                                                 "K11181",
                                                 "K00380",
                                                 "K00381",
                                                 "K16950",
                                                 "K16951",
                                                 "K00385",
                                                 "K00392",
                                                 "K17218",
                                                 "K17995",
                                                 "K17996",
                                                 "K17993",
                                                 "K17994",
                                                 "K17229",
                                                 "K17230",
                                                 "K16952",
                                                 "K17219",
                                                 "K17220",
                                                 "K17221",
                                                 "K17725",
                                                 "K16936",
                                                 "K16937",
                                                 "K08357",
                                                 "K08358",
                                                 "K08359",
                                                 "K08352",
                                                 "K08353",
                                                 "K08354",
                                                 "K01011",
                                                 "K02439",
                                                 "K00640",
                                                 "K10150",
                                                 "K01738",
                                                 "K13034",
                                                 "K17069",
                                                 "K00651",
                                                 "K00641",
                                                 "K01739",
                                                 "K10764",
                                                 "K17217",
                                                 "K17285",
                                                 "K17228",
                                                 "K16968",
                                                 "K16969",
                                                 "K15762",
                                                 "K15765",
                                                 "K07306",
                                                 "K07307",
                                                 "K07308",
                                                 "K00184",
                                                 "K00185",
                                                 "K16964",
                                                 "K16965",
                                                 "K16966",
                                                 "K16953",
                                                 "K17486",
                                                 "K20034",
                                                 "K20035",
                                                 "K20036",
                                                 "K16967",
                                                 "K16954",
                                                 "K16955",
                                                 "K16157",
                                                 "K16158",
                                                 "K16159",
                                                 "K16160",
                                                 "K16161",
                                                 "K16162",
                                                 "K10944",
                                                 "K10945",
                                                 "K10946",
                                                 "K14028",
                                                 "K16254",
                                                 "K16255",
                                                 "K14029",
                                                 "K16256",
                                                 "K16257",
                                                 "K16258",
                                                 "K16259",
                                                 "K16260",
                                                 "K17066",
                                                 "K00093",
                                                 "K00148",
                                                 "K17067",
                                                 "K17068",
                                                 "K03396",
                                                 "K00121",
                                                 "K01070",
                                                 "K00122",
                                                 "K00123",
                                                 "K00124",
                                                 "K00127",
                                                 "K00126",
                                                 "K00125",
                                                 " K05299",
                                                 " K15022",
                                                 "K00192",
                                                 "K00195",
                                                 "K00193",
                                                 " K00197",
                                                 " K00194",
                                                 " K00198",
                                                 "K00196",
                                                 "K00600",
                                                 "K00830",
                                                 "K00018",
                                                 "K11529",
                                                 "K01689",
                                                 "  K01595",
                                                 "  K00024",
                                                 "K08692",
                                                 "K14067",
                                                 "K08691",
                                                 "K17100",
                                                 "K00863",
                                                 "K01623",
                                                 "K11645",
                                                 "K01624",
                                                 "K01622",
                                                 "K16305",
                                                 "K16306",
                                                 "K03841",
                                                 "K02446",
                                                 "K11532",
                                                 "K01086",
                                                 "K04041",
                                                 "K00850",
                                                 "K16370",
                                                 "K00918",
                                                 "K08094",
                                                 "K08093",
                                                 "K13812",
                                                 "K13831",
                                                 "K00317",
                                                 "K18277",
                                                 "K07811",
                                                 "K03532",
                                                 "K03533",
                                                 "K07821",
                                                 "K07812",
                                                 "K15228",
                                                 "K15229",
                                                 "K08685",
                                                 "K00200",
                                                 "K00201",
                                                 "K00202",
                                                 "K00203",
                                                 "K00204",
                                                 "K00205",
                                                 "K11260",
                                                 "K11261",
                                                 "K00672",
                                                 "K01499",
                                                 "K00319",
                                                 "K00440",
                                                 "K00441",
                                                 "K00442",
                                                 "K00443",
                                                 "K10714",
                                                 "K13942",
                                                 "K10713",
                                                 "K00320",
                                                 "K00577",
                                                 "K00578",
                                                 "K00579",
                                                 "K00580",
                                                 "K00581",
                                                 "K00582",
                                                 "K00583",
                                                 "K00584",
                                                 "K00399",
                                                 "K00400",
                                                 "K00401",
                                                 "K00402",
                                                 "K03421",
                                                 "K03422",
                                                 "K03388",
                                                 "K03389",
                                                 "K03390",
                                                 "K08264",
                                                 "K08265",
                                                 "K14127",
                                                 "K14126",
                                                 "K14128",
                                                 "K00925",
                                                 "K00625",
                                                 "K13788",
                                                 "K01895",
                                                 "  K00169",
                                                 "  K00170",
                                                 "  K00172",
                                                 "  K00171",
                                                 "  K01007",
                                                 "K01834",
                                                 "K15633",
                                                 "K15634",
                                                 "K15635",
                                                 "K00058",
                                                 "K00831",
                                                 "K01079",
                                                 "K02203",
                                                 "K14080",
                                                 "K04480",
                                                 "K14081",
                                                 "K14082",
                                                 "K14083",
                                                 "K14084",
                                                 "K16176",
                                                 "K16177",
                                                 "K16178",
                                                 "K16179",
                                                 "K08097",
                                                 "K05979",
                                                 "K05884",
                                                 "K06034",
                                                 "K13039",
                                                 "K11779",
                                                 "K11781",
                                                 "K11780",
                                                 "K14941",
                                                 "K11212",
                                                 "K12234",
                                                 "K14940",
                                                 "K10977",
                                                 "K16792",
                                                 "K16793",
                                                 "K10978",
                                                 "K18933",
                                                 "K06914",
                                                 "K09733",
                                                 "K19793",
                                                 "K07144",
                                                 "K07072")

#and now for the tree

tr<-read.newick("TreeIrciniaTara_git.nwk")

TaraMags<-tr$tip.label[grep("*GCA*",(tr$tip.label))]
Irciniatree<-drop.tip(tr,TaraMags)


length(row.names(presabsKOirciniareordered)%in%Irciniatree$tip.label)

tree<-geiger::treedata(Irciniatree,presabsKOirciniareordered,sort=F,warnings=T)$phy
data<-as.data.frame(geiger::treedata(tree,presabsKOirciniareordered,sort=F,warnings=T)$data)

#drop two archaea
reorderedredundnacny<-ReorderData(tree, data, taxa.names="row names")



taxabac<-read.table("gtdbtk.bac120.summary_git.tsv", row.names = 1)

taxabac<-cbind(taxabac,taxabac)
colnames(taxabac)<-c("V2","V3")
nonIRcinia<-c(grep("*GCA*",rownames(taxabac)))

taxabacIrcinia<-taxabac[-nonIRcinia,]


alphaencrows<-c(grep("*Alphaproteo*",taxabacIrcinia$V2))
alphaencmags<-as.vector(row.names(taxabacIrcinia)[alphaencrows])
alphaenc<-findMRCA(tree,alphaencmags, type="node")

Gammaencrows<-c(grep("*Gammaproteo*",taxabacIrcinia$V2))
Gammaencmags<-as.vector(row.names(taxabacIrcinia)[Gammaencrows])
Gammaenc<-findMRCA(tree,Gammaencmags, type="node")

proteoencrows<-c(grep("*Proteobacteria*",taxabacIrcinia$V2))
proteoencmags<-as.vector(row.names(taxabacIrcinia)[proteoencrows])
Proteoenc<-findMRCA(tree,proteoencmags, type="node")


Acidoencrows<-c(grep("*Acidobac*",taxabacIrcinia$V2))
Acidoencmags<-as.vector(row.names(taxabacIrcinia)[Acidoencrows])
Acidoenc<-findMRCA(tree,Acidoencmags, type="node")

Actinobacteriotarows<-c(grep("*Actinobacteriota*",taxabacIrcinia$V2))
Actinobacteriotamags<-as.vector(row.names(taxabacIrcinia)[Actinobacteriotarows])
Actinoenc<-findMRCA(tree,Actinobacteriotamags, type="node")

Bacteroidotaencrows<-c(grep("*Bacteroidota*",taxabacIrcinia$V2))
Bacteroidotaencmags<-as.vector(row.names(taxabacIrcinia)[Bacteroidotaencrows])
Bacenc<-findMRCA(tree,Bacteroidotaencmags, type="node")

Binatotaencrows<-c(grep("*Binatota*",taxabacIrcinia$V2))
Binatotaaencmags<-as.vector(row.names(taxabacIrcinia)[Binatotaencrows])
Binaenc<-findMRCA(tree,Binatotaaencmags, type="node")

chloroBencrows<-c(grep("*Chloroflexota_B*",taxabacIrcinia$V2))
chloroBencmags<-as.vector(row.names(taxabacIrcinia)[chloroBencrows])
chloroBenc<-findMRCA(tree,chloroBencmags, type="node")


chloroencrows<-c(grep("*Chloroflexota*",taxabacIrcinia$V2))
chloroencmags<-as.vector(row.names(taxabacIrcinia)[chloroencrows])
chloroenc<-findMRCA(tree,chloroencmags, type="node")

Cyanobacteriaencrows<-c(grep("*Cyanobacteria*",taxabacIrcinia$V2))
Cyanobacteriaencmags<-as.vector(row.names(taxabacIrcinia)[Cyanobacteriaencrows])
CyanoBenc<-findMRCA(tree,Cyanobacteriaencmags, type="node")

Dadaencrows<-c(grep("*Dadabacteria*",taxabacIrcinia$V2))
Dadaencmags<-as.vector(row.names(taxabacIrcinia)[Dadaencrows])
Dadaenc<-findMRCA(tree,Dadaencmags, type="node")

Gemmatimonadotacrows<-c(grep("*Gemmatimonadota*",taxabacIrcinia$V2))
Gemmatimonadotamags<-as.vector(row.names(taxabacIrcinia)[Gemmatimonadotacrows])
Gemmatimonadotaenc<-findMRCA(tree,Gemmatimonadotamags, type="node")

Latescibacterotarows<-c(grep("*Latescibacterota*",taxabacIrcinia$V2))
Latescibacterotamags<-as.vector(row.names(taxabacIrcinia)[Latescibacterotarows])
Latescibacterotaenc<-findMRCA(tree,Latescibacterotamags, type="node")

Nitrospirotarows<-c(grep("*Nitrospirota*",taxabacIrcinia$V2))
Nitrospirotamags<-as.vector(row.names(taxabacIrcinia)[Nitrospirotarows])
Nitrospirotaenc<-findMRCA(tree,Nitrospirotamags, type="node")

Poribacteriarows<-c(grep("*Poribacteria*",taxabacIrcinia$V2))
Poribacteriamags<-as.vector(row.names(taxabacIrcinia)[Poribacteriarows])
Poribacteriaenc<-findMRCA(tree,Poribacteriamags, type="node")

UBA8248rows<-c(grep("*UBA8248*",taxabacIrcinia$V2))
UBA8248mags<-as.vector(row.names(taxabacIrcinia)[UBA8248rows])
UBA8248magsenc<-findMRCA(tree,UBA8248mags, type="node")

Verrucomicrobiotarows<-c(grep("*Verrucomicrobiota*",taxabacIrcinia$V2))
Verrucomicrobiotamags<-as.vector(row.names(taxabacIrcinia)[Verrucomicrobiotarows])
Verrucomicrobiotaenc<-findMRCA(tree,Verrucomicrobiotamags, type="node")



p = ggtree(tree)
edge=data.frame(tree$edge, edge_num=1:nrow(tree$edge))
colnames(edge)=c("parent", "node", "edge_num")
p%<+% edge + geom_label(aes(x=branch, label=edge_num))

coverm <- read.csv("covermoutputFullMagSet_noapostrophe.csv", header = TRUE, row.names = 1)
 

#these are the percents of reads that mapped

covermsummary<-c(55.91,54.26,76.64,68.53,55.12,79.27,29.97,11.29,68.46,53.22,82.84,54.68,48.31,57.93,68.92,37.04,54.66,69.68,71.82,54.16)
sd(covermsummary)
mean(covermsummary)
colnames(coverm)
#let's reorganize to Florida, Panama, Belize
 
sorted<-c("jk18x8","jk18x11","jk18x7","jk18x9","p16x33","p16x34","p16x57","p16x60","p16x63","p16x65","p16x43","p16x44","jk18x27","jk18x34" ,"jk18x22","jk18x25","jk18x1bz","jk18xA","jk18x40","jk18x44")


covermsorted<-coverm[,sorted]    
#write.csv(covermsorted, "covermsorted.csv")   

covermsortedmapped<-covermsorted[-1,]

goodrownamesMAGscoverm<-row.names(covermsortedmapped)%in%tree$tip.label
covermsortedmappedgood<-covermsortedmapped[goodrownamesMAGscoverm,]
covermsortedmappedgoodlog<-log(covermsortedmappedgood+1)
 
colnames(covermsortedmappedgoodlog)

colnames(covermsortedmappedgoodlog)<-c("ICam (jk18x8)",   "ICam (jk18x11)",  "IcfR (jk18x7)" ,  "IcfR (jk18x9)",   "IRad (p16x33)" ,  "IRad (p16x34)"  , "ILae (p16x57)"  , "ILae (p16x60)"  , "IBoc (p16x63)"  ,
  "IBoc (p16x65)",   "ILow (p16x43)" ,  "ILow (p16x44)"  , "IVan (jk18x27)",  "IVan (jk18x34)" , "IRue (jk18x22)",  "IRue (jk18x25)" , "IFel (jk18x1bz)" ,"IFel (jk18xA)"  ,
  "IStr (jk18x40)" , "IStr (jk18x44)")
 
 
p <- ggtree(tree) + 
  geom_tiplab(size=1.5, align=TRUE, linesize=.5) +
  geom_label2(aes(subset=(node==alphaenc),x=branch, label="Alphaproteobacteria"), fill='grey80', size=4.5)+
  geom_label2(aes(subset=(node==Gammaenc),x=branch, label="Gammaproteobacteria"), fill='grey80', size=4.5)+
  geom_label2(aes(subset=(node==Proteoenc),x=branch, label="Proteobacteria"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Gemmatimonadotaenc),x=branch, label="Gemmatimonadota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Latescibacterotaenc),x=branch, label="Latescibacterota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Nitrospirotaenc),x=branch, label="Nitrospirota"), fill='white', size=4.5,nudge_x = 0.075)+
  geom_label2(aes(subset=(node==Poribacteriaenc),x=branch, label="Poribacteria"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==UBA8248magsenc),x=branch, label="UBA8248"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Acidoenc),x=branch, label="Acidobacteriota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Actinoenc),x=branch, label="Actinobacteriota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Binaenc),x=branch, label="Binatota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==chloroenc),x=branch, label="Chloroflexota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==CyanoBenc),x=branch, label="Cyanobacteria"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Bacenc),x=branch, label="Bacteroidota"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==Dadaenc),x=branch, label="Dadabacteria"), fill='white', size=4.5)+
  geom_label2(aes(subset=(node==130),x=branch, label="Verrucomicrobiota"), fill='white', size=4.5,nudge_x = -0.2)


data<-read.csv("dataforplottinggeneredundancy_git.csv", row.names = 1) #created this file by replacing '1s' in the pres/abs matrix with the name of the function


library(ggnewscale)
p1<-gheatmap(p, covermsortedmappedgoodlog, offset=0.085, width=0.175, font.size=1.5,colnames_angle=45, hjust=1, colnames=TRUE)+
  scale_fill_continuous(low="white",high="black")
p2 <- p1 + new_scale_fill()
q<-gheatmap(p2, data, offset=0.45, width=4, font.size=1.5,colnames_angle=45, hjust=1,
            colnames=TRUE) + scale_fill_manual(values=c("white","paleturquoise3","lavenderblush2","dark green","dark blue","red","khaki3",
                                                                                "indianred4","yellow2","mediumorchid3","green"), name="presence/absence")



svg("./Figure_S1.svg",height=30, width=80)
q
dev.off()



#now figure S2

covermsortedmapped<-read.csv("covermsortednoarchaea.csv", row.names=1)
covermsortedmappedgood<-covermsortedmapped
covermsortedmappedgoodlog<-log(covermsortedmappedgood+1)

reorderedcovermt<-t(covermsortedmappedgoodlog)

max(reorderedcovermt)
reorderedcovermt[,1]



res.pca <- prcomp(reorderedcovermt, scale = F)
pca.log1<-prcomp(reorderedcovermt,scale=F)
eigen1<-fviz_eig(res.pca)
scores1<-scores(pca.log1)[,1:2]

row.names(reorderedcovermt)



hostcovermcolors<-c("red",   "red",  "blue",   "blue"  ,
                    "maroon1",   "maroon1" ,  "green" ,  "green"  ,
                    "orange3" ,  "orange3" ,  "gray"   ,"gray"  ,
                    "turquoise1",  "turquoise1",  "purple1"  ,"purple1" ,
                    "orange", "orange" ,  "black" , "black" ) 

hostcovermpch<-c(15,   15,   15 ,15,   16,   16 ,  16 , 
                 16 ,  16  , 16 ,  16,  16  ,17 , 17 ,
                 17 , 17  ,17 ,  17,   17,  17)


myplot<-autoplot(pca.log1 , size=6, loadings = T,loadings.label = T, label.size = 1,loadings.label.size=3.5, cex=1, colour=hostcovermcolors, shape=hostcovermpch)


svg("FigureS2C.svg", height=10,width=10)
myplot + theme(panel.background = element_rect(fill = 'transparent', colour = 'black'))
dev.off()


#now the other panels for figure S2

reorderedcovermtcampana<-colSums(reorderedcovermt[1:2,]) 
reorderedcovermtramose<-colSums(reorderedcovermt[3:4,]) 
reorderedcovermtpink<-colSums(reorderedcovermt[5:6,] )
reorderedcovermtgreen<-colSums(reorderedcovermt[7:8,])
reorderedcovermtB<-colSums(reorderedcovermt[9:10,] )
reorderedcovermdataencrusting<-colSums(reorderedcovermt[11:12,]) 
reorderedcovermt1<-colSums(reorderedcovermt[13:14,] )
reorderedcovermt2<-colSums(reorderedcovermt[15:16,] )
reorderedcovermtfelix<-colSums(reorderedcovermt[17:18,] )
reorderedcovermtstrob<-colSums(reorderedcovermt[19:20,])

reorderedcovermtFLORIDA<-colSums(reorderedcovermt[1:4,]) 
reorderedcovermtPANAMA<-colSums(reorderedcovermt[5:12,]) 
reorderedcovermtBELIZE<-colSums(reorderedcovermt[13:20,]) 


presabscamp<-as.matrix((reorderedcovermtcampana > 0) + 0)
presabsramose<-as.matrix((reorderedcovermtramose > 0) + 0)
presabspink<-as.matrix((reorderedcovermtpink > 0) + 0)
presabsgreen<-as.matrix((reorderedcovermtgreen > 0) + 0)
presabsB<-as.matrix((reorderedcovermtB > 0) + 0)
presabsencrusting<-as.matrix((reorderedcovermdataencrusting > 0) + 0)
presabs1<-as.matrix((reorderedcovermt1 > 0) + 0)
presabs2<-as.matrix((reorderedcovermt2 > 0) + 0)
presabsfelix<-as.matrix((reorderedcovermtfelix > 0) + 0)
presabsstrob<-as.matrix((reorderedcovermtstrob > 0) + 0)

presabscovermdatatFLORIDA<-as.matrix((reorderedcovermtFLORIDA > 0) + 0)
presabscovermdatatPANAMA<-as.matrix((reorderedcovermtPANAMA > 0) + 0)
presabscovermdatatBELIZE<-as.matrix((reorderedcovermtBELIZE > 0) + 0)

newmatrixofcovermbypresabsbysite<-cbind(presabscovermdatatFLORIDA,presabscovermdatatPANAMA,presabscovermdatatBELIZE)

rowsumscovermbysite<-rowSums(newmatrixofcovermbypresabsbysite)

bysite<-c(length(which(rowsumscovermbysite==3))/424,
          length(which(rowsumscovermbysite==2))/424,
          length(which(rowsumscovermbysite==1))/424)


newmatrixofcovermbypresabs<-cbind(presabscamp,presabsramose,presabspink,presabsgreen,presabsB,
                                  presabsencrusting,presabs1,presabs2,presabsfelix, presabsstrob)

rowsumscoverm<-rowSums(newmatrixofcovermbypresabs)
vectorofrowsums<-c(length(which(rowsumscoverm==10))/424,
                   length(which(rowsumscoverm==9))/424,
                   length(which(rowsumscoverm==8))/424,
                   length(which(rowsumscoverm==7))/424,
                   length(which(rowsumscoverm==6))/424,
                   length(which(rowsumscoverm==5))/424,
                   length(which(rowsumscoverm==4))/424,
                   length(which(rowsumscoverm==3))/424,
                   length(which(rowsumscoverm==2))/424,
                   length(which(rowsumscoverm==1))/424)
namesofrowsums<-c("10","9","8","7","6","5","4","3","2","1")



allproportionsforcoverm<-c(vectorofrowsums,bysite)

proportionsofhostspeciesbysite<-matrix(ncol=3,bysite, byrow=T)
proportionsofhostspecies<-matrix(ncol=10,vectorofrowsums, byrow=T)

colnames(proportionsofhostspecies)=c("10","9","8","7","6","5","4","3","2","1")
colnames(proportionsofhostspeciesbysite)=c("3","2","1")

proportionsofhostspecies<-proportionsofhostspecies*100
proportionsofhostspeciesbysite<-proportionsofhostspeciesbysite*100

svg("./Fig_S2A.svg", height=10,width=7)
barplot(proportionsofhostspecies, ylim = c(0, 20), col="dark blue", cex.axis = 2, cex.names = 2)
dev.off()
svg("./Fig_S2B.svg", height=10,width=3)
barplot(proportionsofhostspeciesbysite, ylim = c(0, 70), col="orange", cex.axis = 2, cex.names = 2)
dev.off()


sum(proportionsofhostspecies[1:6]) #77.4 percent of MAGs in 5 species or more
proportionsofhostspecies[10] #3.3 % of MAGs found in only one host species
proportionsofhostspeciesbysite[3] #7.07% of MAGs found in only one site






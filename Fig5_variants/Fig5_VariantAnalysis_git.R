setwd("/path/to/wd")

#read in bed files. Important fields are:
#V1: contig of gene
#V2: start of gene
#V3: end of gene
#V4: KO annotation
#V7: the number of variants in the VCF file that overlap the gene position
#V8: same as V7, but without INDELs
#V9: length of the gene in bps
#V10: proportion of gene that is variant


jk18x11bins_2<-read.table("jk18x11bins.2.snpcount.bed")
jk18x11bins_34<-read.table("jk18x11bins.34.snpcount.bed")
jk18x11bins_48<-read.table("jk18x11bins.48.snpcount.bed")
jk18x11bins_53<-read.table("jk18x11bins.53.snpcount.bed")
jk18x1bzbins_11<-read.table("jk18x1bzbins.11.snpcount.bed")
jk18x1bzbins_41<-read.table("jk18x1bzbins.41.snpcount.bed")
jk18x1bzbins_61<-read.table("jk18x1bzbins.61.snpcount.bed")
jk18x1bzbins_75<-read.table("jk18x1bzbins.75.snpcount.bed")
jk18x22bins_18<-read.table("jk18x22bins.18.snpcount.bed")
jk18x22bins_30<-read.table("jk18x22bins.30.snpcount.bed")
jk18x22bins_31<-read.table("jk18x22bins.31.snpcount.bed")
jk18x22bins_41<-read.table("jk18x22bins.41.snpcount.bed")
jk18x22bins_52<-read.table("jk18x22bins.52.snpcount.bed")
jk18x22bins_60<-read.table("jk18x22bins.60.snpcount.bed")
jk18x22bins_66<-read.table("jk18x22bins.66.snpcount.bed")
jk18x22bins_71<-read.table("jk18x22bins.71.snpcount.bed")
jk18x22bins_80<-read.table("jk18x22bins.80.snpcount.bed")
jk18x25bins_5<-read.table("jk18x25bins.5.snpcount.bed")
jk18x25bins_61<-read.table("jk18x25bins.61.snpcount.bed")
jk18x25bins_87<-read.table("jk18x25bins.87.snpcount.bed")
jk18x25bins_9<-read.table("jk18x25bins.9.snpcount.bed")
jk18x27bins_25<-read.table("jk18x27bins.25.snpcount.bed")
jk18x27bins_37<-read.table("jk18x27bins.37.snpcount.bed")
jk18x34bins_13<-read.table("jk18x34bins.13.snpcount.bed")
jk18x34bins_15<-read.table("jk18x34bins.15.snpcount.bed")
jk18x34bins_59<-read.table("jk18x34bins.59.snpcount.bed")
jk18x40bins_30<-read.table("jk18x40bins.30.snpcount.bed")
jk18x40bins_47<-read.table("jk18x40bins.47.snpcount.bed")
jk18x7bins_13<-read.table("jk18x7bins.13.snpcount.bed")
jk18x7bins_32<-read.table("jk18x7bins.32.snpcount.bed")
jk18x7bins_52<-read.table("jk18x7bins.52.snpcount.bed")
jk18x7bins_59<-read.table("jk18x7bins.59.snpcount.bed")
jk18x7bins_8<-read.table("jk18x7bins.8.snpcount.bed")
jk18x8bins_7<-read.table("jk18x8bins.7.snpcount.bed")
jk18x9bins_84<-read.table("jk18x9bins.84.snpcount.bed")
jk18x9bins_88<-read.table("jk18x9bins.88.snpcount.bed")
jk18xAbins_10<-read.table("jk18xAbins.10.snpcount.bed")
jk18xAbins_40<-read.table("jk18xAbins.40.snpcount.bed")
jk18xAbins_42<-read.table("jk18xAbins.42.snpcount.bed")
jk18xAbins_44<-read.table("jk18xAbins.44.snpcount.bed")
jk18xAbins_48<-read.table("jk18xAbins.48.snpcount.bed")
p16x33bins_2<-read.table("p16x33bins.2.snpcount.bed")
p16x34bins_24<-read.table("p16x34bins.24.snpcount.bed")
p16x34bins_31<-read.table("p16x34bins.31.snpcount.bed")
p16x34bins_64<-read.table("p16x34bins.64.snpcount.bed")
p16x43bins_4<-read.table("p16x43bins.4.snpcount.bed")
p16x43bins_5<-read.table("p16x43bins.5.snpcount.bed")
p16x44bins_53<-read.table("p16x44bins.53.snpcount.bed")
p16x44bins_6<-read.table("p16x44bins.6.snpcount.bed")
p16x57bins_39<-read.table("p16x57bins.39.snpcount.bed")
p16x60bins_11<-read.table("p16x60bins.11.snpcount.bed")
p16x60bins_27<-read.table("p16x60bins.27.snpcount.bed")
p16x60bins_63<-read.table("p16x60bins.63.snpcount.bed")
p16x60bins_8<-read.table("p16x60bins.8.snpcount.bed")
p16x63bins_1<-read.table("p16x63bins.1.snpcount.bed")
p16x63bins_32<-read.table("p16x63bins.32.snpcount.bed")
p16x63bins_52<-read.table("p16x63bins.52.snpcount.bed")
p16x65bins_14<-read.table("p16x65bins.14.snpcount.bed")
p16x65bins_51<-read.table("p16x65bins.51.snpcount.bed")




#pull out the CSGs
jk18x11bins_2CSG<-jk18x11bins_2[grep("k141_564894jk", jk18x11bins_2$V1),][which(jk18x11bins_2[grep("k141_564894jk", jk18x11bins_2$V1),]$V2>=25738 & jk18x11bins_2[grep("k141_564894jk", jk18x11bins_2$V1),]$V3<=31785),]						
jk18x11bins_34CSG<-jk18x11bins_34[grep("k141_1457093jk", jk18x11bins_34$V1),][which(jk18x11bins_34[grep("k141_1457093jk", jk18x11bins_34$V1),]$V2>=10385 & jk18x11bins_34[grep("k141_1457093jk", jk18x11bins_34$V1),]$V3<=16254),]			
jk18x11bins_48CSG<-jk18x11bins_48[grep("k141_268402jk", jk18x11bins_48$V1),][which(jk18x11bins_48[grep("k141_268402jk", jk18x11bins_48$V1),]$V2>=21803 & jk18x11bins_48[grep("k141_268402jk", jk18x11bins_48$V1),]$V3<=27756),]				
jk18x11bins_53CSG<-jk18x11bins_53[grep("k141_1181153jk", jk18x11bins_53$V1),][which(jk18x11bins_53[grep("k141_1181153jk", jk18x11bins_53$V1),]$V2>=940 & jk18x11bins_53[grep("k141_1181153jk", jk18x11bins_53$V1),]$V3<=6925),]				
jk18x1bzbins_11CSG<-jk18x1bzbins_11[grep("k141_728604jk", jk18x1bzbins_11$V1),][which(jk18x1bzbins_11[grep("k141_728604jk", jk18x1bzbins_11$V1),]$V2>=63839 & jk18x1bzbins_11[grep("k141_728604jk", jk18x1bzbins_11$V1),]$V3<=69750),]		
jk18x1bzbins_41CSG<-jk18x1bzbins_41[grep("k141_709338jk", jk18x1bzbins_41$V1),][which(jk18x1bzbins_41[grep("k141_709338jk", jk18x1bzbins_41$V1),]$V2>=14876 & jk18x1bzbins_41[grep("k141_709338jk", jk18x1bzbins_41$V1),]$V3<=20751),]		
jk18x1bzbins_61CSG<-jk18x1bzbins_61[grep("k141_1020486jk", jk18x1bzbins_61$V1),][which(jk18x1bzbins_61[grep("k141_1020486jk", jk18x1bzbins_61$V1),]$V2>=13557 & jk18x1bzbins_61[grep("k141_1020486jk", jk18x1bzbins_61$V1),]$V3<=19455),]	
jk18x1bzbins_75CSG<-jk18x1bzbins_75[grep("k141_772177jk", jk18x1bzbins_75$V1),][which(jk18x1bzbins_75[grep("k141_772177jk", jk18x1bzbins_75$V1),]$V2>=7549 & jk18x1bzbins_75[grep("k141_772177jk", jk18x1bzbins_75$V1),]$V3<=13481),]		
jk18x22bins_18CSG<-jk18x22bins_18[grep("k141_782728jk", jk18x22bins_18$V1),][which(jk18x22bins_18[grep("k141_782728jk", jk18x22bins_18$V1),]$V2>=9092 & jk18x22bins_18[grep("k141_782728jk", jk18x22bins_18$V1),]$V3<=14911),]				
jk18x22bins_30CSG<-jk18x22bins_30[grep("k141_14565jk", jk18x22bins_30$V1),][which(jk18x22bins_30[grep("k141_14565jk", jk18x22bins_30$V1),]$V2>=27218 & jk18x22bins_30[grep("k141_14565jk", jk18x22bins_30$V1),]$V3<=33126),]				
jk18x22bins_31CSG<-jk18x22bins_31[grep("k141_125794jk", jk18x22bins_31$V1),][which(jk18x22bins_31[grep("k141_125794jk", jk18x22bins_31$V1),]$V2>=245 & jk18x22bins_31[grep("k141_125794jk", jk18x22bins_31$V1),]$V3<=6174),]				
jk18x22bins_41CSG<-jk18x22bins_41[grep("k141_813214jk", jk18x22bins_41$V1),][which(jk18x22bins_41[grep("k141_813214jk", jk18x22bins_41$V1),]$V2>=2732 & jk18x22bins_41[grep("k141_813214jk", jk18x22bins_41$V1),]$V3<=8685),]				
jk18x22bins_52CSG<-jk18x22bins_52[grep("k141_262221jk", jk18x22bins_52$V1),][which(jk18x22bins_52[grep("k141_262221jk", jk18x22bins_52$V1),]$V2>=14234 & jk18x22bins_52[grep("k141_262221jk", jk18x22bins_52$V1),]$V3<=20091),]				
jk18x22bins_60CSG<-jk18x22bins_60[grep("k141_221542jk", jk18x22bins_60$V1),][which(jk18x22bins_60[grep("k141_221542jk", jk18x22bins_60$V1),]$V2>=21103 & jk18x22bins_60[grep("k141_221542jk", jk18x22bins_60$V1),]$V3<=27029),]				
jk18x22bins_66CSG<-jk18x22bins_66[grep("k141_431527jk", jk18x22bins_66$V1),][which(jk18x22bins_66[grep("k141_431527jk", jk18x22bins_66$V1),]$V2>=6375 & jk18x22bins_66[grep("k141_431527jk", jk18x22bins_66$V1),]$V3<=12376),]				
jk18x22bins_71CSG<-jk18x22bins_71[grep("k141_622183jk", jk18x22bins_71$V1),][which(jk18x22bins_71[grep("k141_622183jk", jk18x22bins_71$V1),]$V2>=45332 & jk18x22bins_71[grep("k141_622183jk", jk18x22bins_71$V1),]$V3<=51240),]				
jk18x22bins_80CSG<-jk18x22bins_80[grep("k141_124588jk", jk18x22bins_80$V1),][which(jk18x22bins_80[grep("k141_124588jk", jk18x22bins_80$V1),]$V2>=51867 & jk18x22bins_80[grep("k141_124588jk", jk18x22bins_80$V1),]$V3<=57734),]				
jk18x25bins_5CSG<-jk18x25bins_5[grep("k141_720591jk", jk18x25bins_5$V1),][which(jk18x25bins_5[grep("k141_720591jk", jk18x25bins_5$V1),]$V2>=10343 & jk18x25bins_5[grep("k141_720591jk", jk18x25bins_5$V1),]$V3<=16263),]						
jk18x25bins_61CSG<-jk18x25bins_61[grep("k141_375062jk", jk18x25bins_61$V1),][which(jk18x25bins_61[grep("k141_375062jk", jk18x25bins_61$V1),]$V2>=37015 & jk18x25bins_61[grep("k141_375062jk", jk18x25bins_61$V1),]$V3<=42963),]				
jk18x25bins_87CSG<-jk18x25bins_87[grep("k141_574470jk", jk18x25bins_87$V1),][which(jk18x25bins_87[grep("k141_574470jk", jk18x25bins_87$V1),]$V2>=55392 & jk18x25bins_87[grep("k141_574470jk", jk18x25bins_87$V1),]$V3<=61289),]				
jk18x25bins_9CSG<-jk18x25bins_9[grep("k141_642020jk", jk18x25bins_9$V1),][which(jk18x25bins_9[grep("k141_642020jk", jk18x25bins_9$V1),]$V2>=187284 & jk18x25bins_9[grep("k141_642020jk", jk18x25bins_9$V1),]$V3<=193142),]					
jk18x27bins_25CSG<-jk18x27bins_25[grep("k141_1081210jk", jk18x27bins_25$V1),][which(jk18x27bins_25[grep("k141_1081210jk", jk18x27bins_25$V1),]$V2>=19680 & jk18x27bins_25[grep("k141_1081210jk", jk18x27bins_25$V1),]$V3<=25612),]			
jk18x27bins_37CSG<-jk18x27bins_37[grep("k141_433742jk", jk18x27bins_37$V1),][which(jk18x27bins_37[grep("k141_433742jk", jk18x27bins_37$V1),]$V2>=198592 & jk18x27bins_37[grep("k141_433742jk", jk18x27bins_37$V1),]$V3<=204521),]			
jk18x34bins_13CSG<-jk18x34bins_13[grep("k141_70560jk", jk18x34bins_13$V1),][which(jk18x34bins_13[grep("k141_70560jk", jk18x34bins_13$V1),]$V2>=136245 & jk18x34bins_13[grep("k141_70560jk", jk18x34bins_13$V1),]$V3<=142108),]				
jk18x34bins_15CSG<-jk18x34bins_15[grep("k141_530823jk", jk18x34bins_15$V1),][which(jk18x34bins_15[grep("k141_530823jk", jk18x34bins_15$V1),]$V2>=42882 & jk18x34bins_15[grep("k141_530823jk", jk18x34bins_15$V1),]$V3<=48824),]				
jk18x34bins_59CSG<-jk18x34bins_59[grep("k141_187261jk", jk18x34bins_59$V1),][which(jk18x34bins_59[grep("k141_187261jk", jk18x34bins_59$V1),]$V2>=389239 & jk18x34bins_59[grep("k141_187261jk", jk18x34bins_59$V1),]$V3<=395127),]			
jk18x40bins_30CSG<-jk18x40bins_30[grep("k141_1301043jk", jk18x40bins_30$V1),][which(jk18x40bins_30[grep("k141_1301043jk", jk18x40bins_30$V1),]$V2>=4518 & jk18x40bins_30[grep("k141_1301043jk", jk18x40bins_30$V1),]$V3<=10405),]			
jk18x40bins_47CSG<-jk18x40bins_47[grep("k141_2697227jk", jk18x40bins_47$V1),][which(jk18x40bins_47[grep("k141_2697227jk", jk18x40bins_47$V1),]$V2>=28599 & jk18x40bins_47[grep("k141_2697227jk", jk18x40bins_47$V1),]$V3<=34511),]			
jk18x7bins_13CSG<-jk18x7bins_13[grep("k141_104580jk", jk18x7bins_13$V1),][which(jk18x7bins_13[grep("k141_104580jk", jk18x7bins_13$V1),]$V2>=23106 & jk18x7bins_13[grep("k141_104580jk", jk18x7bins_13$V1),]$V3<=29002),]						
jk18x7bins_32CSG<-jk18x7bins_32[grep("k141_2384862jk", jk18x7bins_32$V1),][which(jk18x7bins_32[grep("k141_2384862jk", jk18x7bins_32$V1),]$V2>=20660 & jk18x7bins_32[grep("k141_2384862jk", jk18x7bins_32$V1),]$V3<=26562),]					
jk18x7bins_52CSG<-jk18x7bins_52[grep("k141_2166396jk", jk18x7bins_52$V1),][which(jk18x7bins_52[grep("k141_2166396jk", jk18x7bins_52$V1),]$V2>=10909 & jk18x7bins_52[grep("k141_2166396jk", jk18x7bins_52$V1),]$V3<=16808),]					
jk18x7bins_59CSG<-jk18x7bins_59[grep("k141_469756jk", jk18x7bins_59$V1),][which(jk18x7bins_59[grep("k141_469756jk", jk18x7bins_59$V1),]$V2>=-1 & jk18x7bins_59[grep("k141_469756jk", jk18x7bins_59$V1),]$V3<=5295),]							
jk18x7bins_8CSG<-jk18x7bins_8[grep("k141_1123632jk", jk18x7bins_8$V1),][which(jk18x7bins_8[grep("k141_1123632jk", jk18x7bins_8$V1),]$V2>=42743 & jk18x7bins_8[grep("k141_1123632jk", jk18x7bins_8$V1),]$V3<=48672),]							
jk18x8bins_7CSG<-jk18x8bins_7[grep("k141_98869jk", jk18x8bins_7$V1),][which(jk18x8bins_7[grep("k141_98869jk", jk18x8bins_7$V1),]$V2>=3231 & jk18x8bins_7[grep("k141_98869jk", jk18x8bins_7$V1),]$V3<=9172),]									
jk18x9bins_84CSG<-jk18x9bins_84[grep("k141_221715jk", jk18x9bins_84$V1),][which(jk18x9bins_84[grep("k141_221715jk", jk18x9bins_84$V1),]$V2>=16844 & jk18x9bins_84[grep("k141_221715jk", jk18x9bins_84$V1),]$V3<=22776),]						
jk18x9bins_88CSG<-jk18x9bins_88[grep("k141_431810jk", jk18x9bins_88$V1),][which(jk18x9bins_88[grep("k141_431810jk", jk18x9bins_88$V1),]$V2>=5759 & jk18x9bins_88[grep("k141_431810jk", jk18x9bins_88$V1),]$V3<=11670),]						
jk18xAbins_10CSG<-jk18xAbins_10[grep("k141_203481jk", jk18xAbins_10$V1),][which(jk18xAbins_10[grep("k141_203481jk", jk18xAbins_10$V1),]$V2>=142875 & jk18xAbins_10[grep("k141_203481jk", jk18xAbins_10$V1),]$V3<=148723),]					
jk18xAbins_40CSG<-jk18xAbins_40[grep("k141_1094237jk", jk18xAbins_40$V1),][which(jk18xAbins_40[grep("k141_1094237jk", jk18xAbins_40$V1),]$V2>=1192 & jk18xAbins_40[grep("k141_1094237jk", jk18xAbins_40$V1),]$V3<=7100),]					
jk18xAbins_42CSG<-jk18xAbins_42[grep("k141_212514jk", jk18xAbins_42$V1),][which(jk18xAbins_42[grep("k141_212514jk", jk18xAbins_42$V1),]$V2>=20531 & jk18xAbins_42[grep("k141_212514jk", jk18xAbins_42$V1),]$V3<=26532),]						
jk18xAbins_44CSG<-jk18xAbins_44[grep("k141_952455jk", jk18xAbins_44$V1),][which(jk18xAbins_44[grep("k141_952455jk", jk18xAbins_44$V1),]$V2>=1448 & jk18xAbins_44[grep("k141_952455jk", jk18xAbins_44$V1),]$V3<=7365),]						
jk18xAbins_48CSG<-jk18xAbins_48[grep("k141_1014631jk", jk18xAbins_48$V1),][which(jk18xAbins_48[grep("k141_1014631jk", jk18xAbins_48$V1),]$V2>=9041 & jk18xAbins_48[grep("k141_1014631jk", jk18xAbins_48$V1),]$V3<=14994),]					
p16x33bins_2CSG<-p16x33bins_2[grep("k141_288253p16", p16x33bins_2$V1),][which(p16x33bins_2[grep("k141_288253p16", p16x33bins_2$V1),]$V2>=26642 & p16x33bins_2[grep("k141_288253p16", p16x33bins_2$V1),]$V3<=32517),]							
p16x34bins_24CSG<-p16x34bins_24[grep("k141_677482p16", p16x34bins_24$V1),][which(p16x34bins_24[grep("k141_677482p16", p16x34bins_24$V1),]$V2>=349041 & p16x34bins_24[grep("k141_677482p16", p16x34bins_24$V1),]$V3<=354961),]				
p16x34bins_31CSG<-p16x34bins_31[grep("k141_18952p16", p16x34bins_31$V1),][which(p16x34bins_31[grep("k141_18952p16", p16x34bins_31$V1),]$V2>=1083 & p16x34bins_31[grep("k141_18952p16", p16x34bins_31$V1),]$V3<=6994),]						
p16x34bins_64CSG<-p16x34bins_64[grep("k141_119330p16", p16x34bins_64$V1),][which(p16x34bins_64[grep("k141_119330p16", p16x34bins_64$V1),]$V2>=2324 & p16x34bins_64[grep("k141_119330p16", p16x34bins_64$V1),]$V3<=8226),]					
p16x43bins_4CSG<-p16x43bins_4[grep("k141_886407p16", p16x43bins_4$V1),][which(p16x43bins_4[grep("k141_886407p16", p16x43bins_4$V1),]$V2>=7502 & p16x43bins_4[grep("k141_886407p16", p16x43bins_4$V1),]$V3<=13506),]							
p16x43bins_5CSG<-p16x43bins_5[grep("k141_63071p16", p16x43bins_5$V1),][which(p16x43bins_5[grep("k141_63071p16", p16x43bins_5$V1),]$V2>=17790 & p16x43bins_5[grep("k141_63071p16", p16x43bins_5$V1),]$V3<=23737),]								
p16x44bins_53CSG<-p16x44bins_53[grep("k141_524271p16", p16x44bins_53$V1),][which(p16x44bins_53[grep("k141_524271p16", p16x44bins_53$V1),]$V2>=81872 & p16x44bins_53[grep("k141_524271p16", p16x44bins_53$V1),]$V3<=87825),]					
p16x44bins_6CSG<-p16x44bins_6[grep("k141_205191p16", p16x44bins_6$V1),][which(p16x44bins_6[grep("k141_205191p16", p16x44bins_6$V1),]$V2>=2361 & p16x44bins_6[grep("k141_205191p16", p16x44bins_6$V1),]$V3<=8215),]							
p16x57bins_39CSG<-p16x57bins_39[grep("k141_792160p16", p16x57bins_39$V1),][which(p16x57bins_39[grep("k141_792160p16", p16x57bins_39$V1),]$V2>=23624 & p16x57bins_39[grep("k141_792160p16", p16x57bins_39$V1),]$V3<=29541),]					
p16x60bins_11CSG<-p16x60bins_11[grep("k141_339243p16", p16x60bins_11$V1),][which(p16x60bins_11[grep("k141_339243p16", p16x60bins_11$V1),]$V2>=47803 & p16x60bins_11[grep("k141_339243p16", p16x60bins_11$V1),]$V3<=53788),]					
p16x60bins_27CSG<-p16x60bins_27[grep("k141_367051p16", p16x60bins_27$V1),][which(p16x60bins_27[grep("k141_367051p16", p16x60bins_27$V1),]$V2>=2019 & p16x60bins_27[grep("k141_367051p16", p16x60bins_27$V1),]$V3<=7923),]					
p16x60bins_63CSG<-p16x60bins_63[grep("k141_212548p16", p16x60bins_63$V1),][which(p16x60bins_63[grep("k141_212548p16", p16x60bins_63$V1),]$V2>=16098 & p16x60bins_63[grep("k141_212548p16", p16x60bins_63$V1),]$V3<=21995),]					
p16x60bins_8CSG<-p16x60bins_8[grep("k141_854900p16", p16x60bins_8$V1),][which(p16x60bins_8[grep("k141_854900p16", p16x60bins_8$V1),]$V2>=11536 & p16x60bins_8[grep("k141_854900p16", p16x60bins_8$V1),]$V3<=17435),]							
p16x63bins_1CSG<-p16x63bins_1[grep("k141_102737p16", p16x63bins_1$V1),][which(p16x63bins_1[grep("k141_102737p16", p16x63bins_1$V1),]$V2>=17427 & p16x63bins_1[grep("k141_102737p16", p16x63bins_1$V1),]$V3<=23437),]							
p16x63bins_32CSG<-p16x63bins_32[grep("k141_835293p16", p16x63bins_32$V1),][which(p16x63bins_32[grep("k141_835293p16", p16x63bins_32$V1),]$V2>=43043 & p16x63bins_32[grep("k141_835293p16", p16x63bins_32$V1),]$V3<=48951),]					
p16x63bins_52CSG<-p16x63bins_52[grep("k141_621503p16", p16x63bins_52$V1),][which(p16x63bins_52[grep("k141_621503p16", p16x63bins_52$V1),]$V2>=449090 & p16x63bins_52[grep("k141_621503p16", p16x63bins_52$V1),]$V3<=454959),]				
p16x65bins_14CSG<-p16x65bins_14[grep("k141_913059p16", p16x65bins_14$V1),][which(p16x65bins_14[grep("k141_913059p16", p16x65bins_14$V1),]$V2>=48746 & p16x65bins_14[grep("k141_913059p16", p16x65bins_14$V1),]$V3<=54621),]					
p16x65bins_51CSG<-p16x65bins_51[grep("k141_66698p16", p16x65bins_51$V1),][which(p16x65bins_51[grep("k141_66698p16", p16x65bins_51$V1),]$V2>=18261 & p16x65bins_51[grep("k141_66698p16", p16x65bins_51$V1),]$V3<=24208),]						


CSGscombined<-rbind(jk18x11bins_2CSG,
                     jk18x11bins_34CSG,
                     jk18x11bins_48CSG,
                     jk18x11bins_53CSG,
                     jk18x1bzbins_11CSG,
                     jk18x1bzbins_41CSG,
                     jk18x1bzbins_61CSG,
                     jk18x1bzbins_75CSG,
                     jk18x22bins_18CSG,
                     jk18x22bins_30CSG,
                     jk18x22bins_31CSG,
                     jk18x22bins_41CSG,
                     jk18x22bins_52CSG,
                     jk18x22bins_60CSG,
                     jk18x22bins_66CSG,
                     jk18x22bins_71CSG,
                     jk18x22bins_80CSG,
                     jk18x25bins_5CSG,
                     jk18x25bins_61CSG,
                     jk18x25bins_87CSG,
                     jk18x25bins_9CSG,
                     jk18x27bins_25CSG,
                     jk18x27bins_37CSG,
                     jk18x34bins_13CSG,
                     jk18x34bins_15CSG,
                     jk18x34bins_59CSG,
                     jk18x40bins_30CSG,
                     jk18x40bins_47CSG,
                     jk18x7bins_13CSG,
                     jk18x7bins_32CSG,
                     jk18x7bins_52CSG,
                     jk18x7bins_59CSG,
                     jk18x7bins_8CSG,
                     jk18x8bins_7CSG,
                     jk18x9bins_84CSG,
                     jk18x9bins_88CSG,
                     jk18xAbins_10CSG,
                     jk18xAbins_40CSG,
                     jk18xAbins_42CSG,
                     jk18xAbins_44CSG,
                     jk18xAbins_48CSG,
                     p16x33bins_2CSG,
                     p16x34bins_24CSG,
                     p16x34bins_31CSG,
                     p16x34bins_64CSG,
                     p16x43bins_4CSG,
                     p16x43bins_5CSG,
                     p16x44bins_53CSG,
                     p16x44bins_6CSG,
                     p16x57bins_39CSG,
                     p16x60bins_11CSG,
                     p16x60bins_27CSG,
                     p16x60bins_63CSG,
                     p16x60bins_8CSG,
                     p16x63bins_1CSG,
                     p16x63bins_32CSG,
                     p16x63bins_52CSG,
                     p16x65bins_14CSG,
                     p16x65bins_51CSG)

write.csv(CSGscombined, "CSGscombined.csv")


meansgenomewide<-c(sum(jk18x11bins_2$V8)/sum(jk18x11bins_2$V9),
                       sum(jk18x11bins_34$V8)/sum(jk18x11bins_34$V9),
                       sum(jk18x11bins_48$V8)/sum(jk18x11bins_48$V9),
                       sum(jk18x11bins_53$V8)/sum(jk18x11bins_53$V9),
                       sum(jk18x1bzbins_11$V8)/sum(jk18x1bzbins_11$V9),
                       sum(jk18x1bzbins_41$V8)/sum(jk18x1bzbins_41$V9),
                       sum(jk18x1bzbins_61$V8)/sum(jk18x1bzbins_61$V9),
                       sum(jk18x1bzbins_75$V8)/sum(jk18x1bzbins_75$V9),
                       sum(jk18x22bins_18$V8)/sum(jk18x22bins_18$V9),
                       sum(jk18x22bins_30$V8)/sum(jk18x22bins_30$V9),
                       sum(jk18x22bins_31$V8)/sum(jk18x22bins_31$V9),
                       sum(jk18x22bins_41$V8)/sum(jk18x22bins_41$V9),
                       sum(jk18x22bins_52$V8)/sum(jk18x22bins_52$V9),
                       sum(jk18x22bins_60$V8)/sum(jk18x22bins_60$V9),
                       sum(jk18x22bins_66$V8)/sum(jk18x22bins_66$V9),
                       sum(jk18x22bins_71$V8)/sum(jk18x22bins_71$V9),
                       sum(jk18x22bins_80$V8)/sum(jk18x22bins_80$V9),
                       sum(jk18x25bins_5$V8)/sum(jk18x25bins_5$V9),
                       sum(jk18x25bins_61$V8)/sum(jk18x25bins_61$V9),
                       sum(jk18x25bins_87$V8)/sum(jk18x25bins_87$V9),
                       sum(jk18x25bins_9$V8)/sum(jk18x25bins_9$V9),
                       sum(jk18x27bins_25$V8)/sum(jk18x27bins_25$V9),
                       sum(jk18x27bins_37$V8)/sum(jk18x27bins_37$V9),
                       sum(jk18x34bins_13$V8)/sum(jk18x34bins_13$V9),
                       sum(jk18x34bins_15$V8)/sum(jk18x34bins_15$V9),
                       sum(jk18x34bins_59$V8)/sum(jk18x34bins_59$V9),
                       sum(jk18x40bins_30$V8)/sum(jk18x40bins_30$V9),
                       sum(jk18x40bins_47$V8)/sum(jk18x40bins_47$V9),
                       sum(jk18x7bins_13$V8)/sum(jk18x7bins_13$V9),
                       sum(jk18x7bins_32$V8)/sum(jk18x7bins_32$V9),
                       sum(jk18x7bins_52$V8)/sum(jk18x7bins_52$V9),
                       sum(jk18x7bins_59$V8)/sum(jk18x7bins_59$V9),
                       sum(jk18x7bins_8$V8)/sum(jk18x7bins_8$V9),
                       sum(jk18x8bins_7$V8)/sum(jk18x8bins_7$V9),
                       sum(jk18x9bins_84$V8)/sum(jk18x9bins_84$V9),
                       sum(jk18x9bins_88$V8)/sum(jk18x9bins_88$V9),
                       sum(jk18xAbins_10$V8)/sum(jk18xAbins_10$V9),
                       sum(jk18xAbins_40$V8)/sum(jk18xAbins_40$V9),
                       sum(jk18xAbins_42$V8)/sum(jk18xAbins_42$V9),
                       sum(jk18xAbins_44$V8)/sum(jk18xAbins_44$V9),
                       sum(jk18xAbins_48$V8)/sum(jk18xAbins_48$V9),
                       sum(p16x33bins_2$V8)/sum(p16x33bins_2$V9),
                       sum(p16x34bins_24$V8)/sum(p16x34bins_24$V9),
                       sum(p16x34bins_31$V8)/sum(p16x34bins_31$V9),
                       sum(p16x34bins_64$V8)/sum(p16x34bins_64$V9),
                       sum(p16x43bins_4$V8)/sum(p16x43bins_4$V9),
                       sum(p16x43bins_5$V8)/sum(p16x43bins_5$V9),
                       sum(p16x44bins_53$V8)/sum(p16x44bins_53$V9),
                       sum(p16x44bins_6$V8)/sum(p16x44bins_6$V9),
                       sum(p16x57bins_39$V8)/sum(p16x57bins_39$V9),
                       sum(p16x60bins_11$V8)/sum(p16x60bins_11$V9),
                       sum(p16x60bins_27$V8)/sum(p16x60bins_27$V9),
                       sum(p16x60bins_63$V8)/sum(p16x60bins_63$V9),
                       sum(p16x60bins_8$V8)/sum(p16x60bins_8$V9),
                       sum(p16x63bins_1$V8)/sum(p16x63bins_1$V9),
                       sum(p16x63bins_32$V8)/sum(p16x63bins_32$V9),
                       sum(p16x63bins_52$V8)/sum(p16x63bins_52$V9),
                       sum(p16x65bins_14$V8)/sum(p16x65bins_14$V9),
                       sum(p16x65bins_51$V8)/sum(p16x65bins_51$V9))

meansgenomewidedoublefordataframe<-cbind(meansgenomewide,meansgenomewide)

row.names(meansgenomewidedoublefordataframe)<-c("jk18x11bins_2",
                                 "jk18x11bins_34",
                                 "jk18x11bins_48",
                                 "jk18x11bins_53",
                                 "jk18x1bzbins_11",
                                 "jk18x1bzbins_41",
                                 "jk18x1bzbins_61",
                                 "jk18x1bzbins_75",
                                 "jk18x22bins_18",
                                 "jk18x22bins_30",
                                 "jk18x22bins_31",
                                 "jk18x22bins_41",
                                 "jk18x22bins_52",
                                 "jk18x22bins_60",
                                 "jk18x22bins_66",
                                 "jk18x22bins_71",
                                 "jk18x22bins_80",
                                 "jk18x25bins_5",
                                 "jk18x25bins_61",
                                 "jk18x25bins_87",
                                 "jk18x25bins_9",
                                 "jk18x27bins_25",
                                 "jk18x27bins_37",
                                 "jk18x34bins_13",
                                 "jk18x34bins_15",
                                 "jk18x34bins_59",
                                 "jk18x40bins_30",
                                 "jk18x40bins_47",
                                 "jk18x7bins_13",
                                 "jk18x7bins_32",
                                 "jk18x7bins_52",
                                 "jk18x7bins_59",
                                 "jk18x7bins_8",
                                 "jk18x8bins_7",
                                 "jk18x9bins_84",
                                 "jk18x9bins_88",
                                 "jk18xAbins_10",
                                 "jk18xAbins_40",
                                 "jk18xAbins_42",
                                 "jk18xAbins_44",
                                 "jk18xAbins_48",
                                 "p16x33bins_2",
                                 "p16x34bins_24",
                                 "p16x34bins_31",
                                 "p16x34bins_64",
                                 "p16x43bins_4",
                                 "p16x43bins_5",
                                 "p16x44bins_53",
                                 "p16x44bins_6",
                                 "p16x57bins_39",
                                 "p16x60bins_11",
                                 "p16x60bins_27",
                                 "p16x60bins_63",
                                 "p16x60bins_8",
                                 "p16x63bins_1",
                                 "p16x63bins_32",
                                 "p16x63bins_52",
                                 "p16x65bins_14",
                                 "p16x65bins_51")
write.csv(meansgenomewidedoublefordataframe,"meansgenomewidedoublefordataframe.csv")
#this data frame was then modified manually by adding the proportion of variant sites information (V10) from the CSGscombined.csv file. Taxonomic annotations were then added. Classes that had less than three MAGs were omitted.

CSGsbygeneplusgenome<-read.csv("Dataframe_of_proportions_variant.csv")



colorsforplottingCSGsandMAGs<-rep("gray50",228)
colorsforplottingCSGsandMAGs[1:15]<-"#3C5488FF"
colorsforplottingCSGsandMAGs[16:30]<-"#E64B35FF"
colorsforplottingCSGsandMAGs[31:45]<-"#91D1C2FF"

colorsforplottingCSGsandMAGs[61:82]<-"#3C5488FF"
colorsforplottingCSGsandMAGs[83:104]<-"#E64B35FF"
colorsforplottingCSGsandMAGs[105:126]<-"#91D1C2FF"

colorsforplottingCSGsandMAGs[149:151]<-"#3C5488FF"
colorsforplottingCSGsandMAGs[152:154]<-"#E64B35FF"
colorsforplottingCSGsandMAGs[155:157]<-"#91D1C2FF"

colorsforplottingCSGsandMAGs[161:165]<-"#3C5488FF"
colorsforplottingCSGsandMAGs[166:170]<-"#E64B35FF"
colorsforplottingCSGsandMAGs[171:175]<-"#91D1C2FF"

colorsforplottingCSGsandMAGs[181:184]<-"#3C5488FF"
colorsforplottingCSGsandMAGs[185:188]<-"#E64B35FF"
colorsforplottingCSGsandMAGs[189:192]<-"#91D1C2FF"

colorsforplottingCSGsandMAGs[197:204]<-"#3C5488FF"
colorsforplottingCSGsandMAGs[205:212]<-"#E64B35FF"
colorsforplottingCSGsandMAGs[213:220]<-"#91D1C2FF"


p3ocgcol<- ggplot(CSGsbygeneplusgenome, aes(class_fraction, means)) +
  geom_boxplot(outlier.shape = NA,colour = "grey50") +
  geom_line(aes(group = paired))+
  geom_dotplot(fill=colorsforplottingCSGsandMAGs, binaxis='y',binwidth = 0.001, stackdir="center",dotsize=5)


svg("./Figure_5.svg",height=6, width=15)
p3ocgcol  + theme(panel.background = element_rect(fill = "transparent", colour = NA), 
               panel.grid.minor = element_blank(), 
               panel.grid.major = element_blank(),
               axis.text.y=element_text(size=22),
               axis.title.y=element_text(size=26,face="bold"),
               axis.title.x=element_text(size=26,face="bold"),
               text=element_text(family="Arial"))+
  labs(y="Mean proportion of variable sites", x = "Class")+
  theme(axis.line = element_line(size = 0.25, color = 'black'))+
  theme(axis.ticks.y = element_line(colour = "black", size = 0.5))+
  theme(axis.ticks.x=element_blank())
dev.off()



#permutational t-tests with BH correction

library(Deducer)

Deducer::perm.t.test(CSGsbygeneplusgenome$means[1:15],CSGsbygeneplusgenome$means[16:30],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #alpha K00222 p-value < 2.2e-16
Deducer::perm.t.test(CSGsbygeneplusgenome$means[1:15],CSGsbygeneplusgenome$means[46:60],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #alpha K05917 p-value < 2.2e-16
Deducer::perm.t.test(CSGsbygeneplusgenome$means[1:15],CSGsbygeneplusgenome$means[31:45],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #alpha K01852 p-value < 2.2e-16

p.adjust(c(2.2e-16,2.2e-16,2.2e-16),method="BH",n=3)

Deducer::perm.t.test(CSGsbygeneplusgenome$means[61:82],CSGsbygeneplusgenome$means[83:104],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #gamma K00222 p-value = 0.0036
Deducer::perm.t.test(CSGsbygeneplusgenome$means[61:82],CSGsbygeneplusgenome$means[127:148],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #gamma K05917 p-value = 0.003
Deducer::perm.t.test(CSGsbygeneplusgenome$means[61:82],CSGsbygeneplusgenome$means[105:126],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #gamma K01852 p-value = 0.2516

p.adjust(c(0.0036,0.003,0.2516),method="BH",n=3)

Deducer::perm.t.test(CSGsbygeneplusgenome$means[149:151],CSGsbygeneplusgenome$means[152:154],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #binatia K00222 p-value = 0.0491
Deducer::perm.t.test(CSGsbygeneplusgenome$means[149:151],CSGsbygeneplusgenome$means[158:160],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #binatia K05917 p-value = 0.0477
Deducer::perm.t.test(CSGsbygeneplusgenome$means[149:151],CSGsbygeneplusgenome$means[155:157],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #binatia K01852 p-value = 0.0468

p.adjust(c(0.0491,0.0477,0.0468),method="BH",n=3)

Deducer::perm.t.test(CSGsbygeneplusgenome$means[161:165],CSGsbygeneplusgenome$means[166:170],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #Dehalococcoidia K00222 p-value = 0.3717
Deducer::perm.t.test(CSGsbygeneplusgenome$means[161:165],CSGsbygeneplusgenome$means[176:180],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #Dehalococcoidia K05917 p-value = 0.2256
Deducer::perm.t.test(CSGsbygeneplusgenome$means[161:165],CSGsbygeneplusgenome$means[171:175],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #Dehalococcoidia K01852 p-value = 0.012

p.adjust(c(0.3717,0.2256,0.012),method="BH",n=3)

Deducer::perm.t.test(CSGsbygeneplusgenome$means[181:184],CSGsbygeneplusgenome$means[185:188],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #Nitrospiria K00222 p-value = 0.012
Deducer::perm.t.test(CSGsbygeneplusgenome$means[181:184],CSGsbygeneplusgenome$means[193:196],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #Nitrospiria K05917 p-value = 0.0128
Deducer::perm.t.test(CSGsbygeneplusgenome$means[181:184],CSGsbygeneplusgenome$means[189:192],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #Nitrospiria K01852 p-value = 0.0129

p.adjust(c(0.012,0.0128,0.0129),method="BH",n=3)

Deducer::perm.t.test(CSGsbygeneplusgenome$means[197:204],CSGsbygeneplusgenome$means[205:212],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #Rhodothermia K00222 p-value = 9.999e-05
Deducer::perm.t.test(CSGsbygeneplusgenome$means[197:204],CSGsbygeneplusgenome$means[221:228],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #Rhodothermia K05917 p-value < 2.2e-16
Deducer::perm.t.test(CSGsbygeneplusgenome$means[197:204],CSGsbygeneplusgenome$means[213:220],statistic="t",alternative=c("two.sided"), midp=TRUE, B=10000) #Rhodothermia K01852 p-value < 2.2e-16

p.adjust(c(9.999e-05,2.2e-16,2.2e-16),method="BH",n=3)





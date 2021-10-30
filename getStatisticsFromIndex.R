directory<-"/Users/huebler/Desktop/DatabasePipeline/bacteria/list/index.txt"
directory<-"/Users/huebler/Desktop/index-backup.txt"
directory<-"/Users/huebler/Desktop/index.txt"
data<-as.data.frame(read.delim2(directory,stringsAsFactors = FALSE),rownames=1)

data$wantToDownload
data[as.numeric(data$OffPathRelaxed)>0,]
length(withTax$taxID)
length(withTax[withTax$OnPathRelaxed>=0.5,]$taxID)
length(withTax[withTax$OnPathRelaxed<0.5,]$taxID)
length(withTax[withTax$OnPathRelaxed<0.1,]$OnPathRelaxed)
length(withTax[withTax$OnPathStrict<0.1,]$OnPathRelaxed)

#if species tax id == tax iD this does not work? nope that is not the problem
# also assembly level does not seem to influence this
#reference not reference does not unfluence this as well
# if we say less than 0.5 gets removed strict 7184 kept 4103, Relaxed 8074 kept 3213 removed
# check if there are even species in the genus



ID<-withTax[withTax$taxID==withTax$speciesTaxID, ]
table(withTax[withTax$OnPathRelaxed<0.1, ]$assembly_level)
table(withTax[withTax$OnPathRelaxed>0.5, ]$assembly_level)
table(withTax[withTax$OnPathRelaxed<0.1, ]$reference)
table(withTax[withTax$OnPathRelaxed>0.5, ]$reference)
table(withTax$assembly_level)

withTax[ID,]
length(withTax[ID,]$taxID)
ID<-withTax[withTax$OnPathStrict>0.6,]$taxID==withTax[withTax$OnPathStrict>0.6,]$speciesTaxID
length(ID)

hist(as.numeric(withTax$OnPathStrict))
hist(as.numeric(withTax$OnPathRelaxed))
as.factor(data$wantToKeep)
length(data$wantToDownload)
assemblies<-data[c(data$assembly_level!="COMPLETE"),]
sort(assemblies$asm_name)
plot(x=as.numeric(assemblies$moleculeCount), y=assemblies$assembly_level)
downloaded<-data[c(data$wantToKeep=="true"),]

length(downloaded$assembly_level)
table(downloaded$assembly_level)
table(data$assembly_level)
#molecule count
median(as.numeric(assemblies$moleculeCount))
mean(as.numeric(assemblies$moleculeCount))
hist(as.numeric(assemblies$moleculeCount))
max(as.numeric(assemblies$moleculeCount))

#component count
median(as.numeric(assemblies$componentCount))
mean(as.numeric(assemblies$componentCount))
hist(as.numeric(assemblies$componentCount))


#contig L50
median(as.numeric(assemblies$contigCountL50))
mean(as.numeric(assemblies$contigCountL50))
hist(as.numeric(assemblies$contigCountL50))

#contig N50
median(as.numeric(assemblies$contigCountN50))
mean(as.numeric(assemblies$contigCountN50))
hist(as.numeric(assemblies$contigCountN50))


#uneven distribution mean and median are different
#

median(as.numeric(assemblies$topLevelCount))
mean(as.numeric(assemblies$topLevelCount))
hist(as.numeric(assemblies$topLevelCount))


##TotalGaplength uneven distribution. SO I think contigs we use should not contain gaps
median(as.numeric(assemblies$totalGapLength))
mean(as.numeric(assemblies$totalGapLength))
hist(as.numeric(assemblies$totalGapLength))

#adapter DB funden
PHYx
Vector sequenzen
length(data$PercentageKept)
sort(data$taxID)
data[data$PercentageKept==0,]<-1
percentage<-c(as.double(data$PercentageKept)*100)
hist(percentage,xlim= range(85,100))
c(as.numeric(data$PercentageKept)*100)
non_null<-data[as.numeric(data$PercentageKept)>0,]
non_null<-non_null[as.numeric(non_null$PercentageKept)<1,]
length(non_null$taxID)
table(data$assembly_level)
hist(as.numeric(non_null$PercentageKept))
mean(as.numeric(non_null$PercentageKept))
median(as.numeric(non_null$PercentageKept))
table(nulls$assembly_level)
hist(data$NumberTotalContigs)
length(data$NumberTotalContigs)
length(unique(data$taxID))
length(data$speciesTaxID)
assemblies<-data[c(data$assembly_level!='COMPLETE'),]
as.character.numeric_version()
data$OnPathRelaxed
library(ggplot2)
data$OnPathStrict<-as.numeric(data$OnPathStrict) 
p<-ggplot(data, aes(x=OnPathStrict)) +   geom_histogram(color="black", fill="grey")+ 
  xlab("Fraction of assigned reads\non correct path")+ ylab("")+theme_bw()+
  ggtitle("Fraction of reads assigned\n to correct path (strict)")+ theme(text=element_text(size=16, face="bold"))
p
ggsave(filename = "OPS.pdf",path = "/Users/huebler/Desktop/")
data$OnPathRelaxed<-as.numeric(data$OnPathRelaxed)

p<-ggplot(data, aes(x=OnPathRelaxed)) +   geom_histogram(color="black", fill="red")+ 
  xlab("Fraction of assigned reads\non correct path")+ ylab("")+theme_bw()+
  ggtitle("Fraction of reads assigned\n to correct path (relaxed)")+ theme(text=element_text(size=16, face="bold"))
ggsave(filename = "OPR.pdf",path = "/Users/huebler/Desktop/")

data$PercentageKept<-as.numeric(data$PercentageKept)
p<-ggplot(data, aes(x=PercentageKept)) +   geom_histogram(color="black", fill="grey")+ 
  xlab("Fraction kept of assemply")+ ylab("")+theme_bw()+
  ggtitle("Fraction kept of assemply")+ theme(text=element_text(size=16, face="bold"))
p
ggsave(filename = "Perk.pdf",path = "/Users/huebler/Desktop/")
data$OnPathRelaxed<-as.numeric(data$OnPathRelaxed)

onPath<-as.double(data$OnPath)
onPath[is.na(onPath)]<--1
hist(onPath)
hist(as.double(data$OnPathStrict))
hist(as.double(data$OnPathRelaxed))

offPathStrict.<-as.double(data$OffPathStrict.)
offPath[is.na(offPath)]<--1
unavailable<-data[data$OffPathStrict=="-1",]$FileName
length(unavailable)
unavailable
data$speciesTaxID
treponema
#prefiltering sort by values on NCBI 
# NCBI diversitaets tree phylogeny
# plot histograms of number of kept and removed contigs here 100000 bp
hist(assemblies$NumberTotalContigs)
hist(assemblies$NumberKeptContigs)
hist(assemblies$NumberRemovedContigs)
#how much is kept 
#wie ist das bei S enterica 



hist(data$assembly_level)
str(as.factor(data$assembly_level))

data[data$Adapter==TRUE,]$Name
as.factor(data$Adapter)
table(data[data$Adapter=='true',]$assembly_level)
hist(data[data$Adapter==TRUE,]$AdapterOccurance)
missing_assemblies<-read.csv("/Users/huebler/Desktop/missing_assemblies.txt",stringsAsFactors = FALSE)

mean(as.numeric(data$AdapterOccurance))
median(as.numeric(data$AdapterOccurance))
max(data$AdapterOccurance)
data[data$AdapterOccurance>5,]$Name
df<-data[data$asm_name==missing_assemblies[1,],]

length(missing_assemblies)
while(i<=length(missing_assemblies$missing.assemblies)){
 df<- rbind(df,data[data$asm_name==missing_assemblies[i,],])
  i<-i+1
}
table(data$wantToDownload)
data.

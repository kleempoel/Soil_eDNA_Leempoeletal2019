SamplesFile='Miseq4_M12sSoil_Samples.txt' #List of samples in tab delimited columns.
ShortFile=substr(SamplesFile,1,nchar(SamplesFile)-4)
EcotagFile='Miseq4_M12sSoil_Samples/Combined_Uniq_Sorted_Ecotag_CleanH.c10.txt'
source("../../R_Scripts/R_Obitools.r") #Location of the Obitools R functions
param=data.frame(mincount=1) #Minimum count of same sequences, pre-blast filtering (obitools/obigrep)
param$minlength=150 #Minimum length of sequence (obitools)
param$maxlength=192 #Maximum length of sequence (obitools)



#####1. Load Data
Samples_Processed=read.table(paste(substr(SamplesFile,1,nchar(SamplesFile)-4),"_filecheck.txt",sep=""),header = T,sep = "\t")
Samples=read.table(SamplesFile,header = T,sep="\t")
Data_Ecotag_raw=read.table(EcotagFile,header = T,sep="\t") #Must be sorted by count number
Data_Ecotag=ecotag_fieldchange(Data_Ecotag_raw) #change column names



#####2. Compute Statistics, need to update index_files.txt first
Obistats_commands(Samples_Processed) #Then input commands in obitools
Stats_Samples=Stats_Processing(Samples)

#Plot percentage of reads lost per sample at each processing step 
pdf(paste(ShortFile,"/",ShortFile,"_ReadsPC_Processing.pdf",sep=""),width = 10,height = 5)
matplot(t(Stats_Samples[,(ncol(Stats_Samples)-3):ncol(Stats_Samples)]), type = "l",xaxt="n",ylab = "% of reads")
title(main=paste(ShortFile," Reads% Processing",sep=""))
axis(side=1, at=seq(1, 6, by=1), labels = F); text(x=seq(1, 6, by=1),  par("usr")[3], labels = c("PairAli","PairEnd","Cut","Filtered"), pos =1, xpd = T)
dev.off()
#Histogram of the number of reads per sample after processing
pdf(paste(ShortFile,"/",ShortFile,"_Histogram_Reads_AfterFiltering.pdf",sep=""),width = 10,height = 5)
hist(Stats_Samples$Nr_Filtered,breaks = 20,main=paste(ShortFile," - Reads after filtering",sep=""))
dev.off()

# Output list of samples that have lost at least X percent of their reads during processing
print(cbind(as.character(Stats_Samples$SampleName[which(Stats_Samples$pc_Nr_Filtered<0.6)]),
  round(Stats_Samples$pc_Nr_Filtered[which(Stats_Samples$pc_Nr_Filtered<0.6)],digits = 3)),quote = F)
# Output list of samples that have less than X reads after processing
print(cbind(as.character(Stats_Samples$SampleName[which(Stats_Samples$Nr_Filtered<50000)]),
            Stats_Samples$Nr_Filtered[which(Stats_Samples$Nr_Filtered<50000)]),quote = F)

#Combine different tables
FA=read.table("FA_Miseq4.txt",header = T,sep="\t")
Samples_AllData=Add2Table(Stats_Samples,FA,"SampleName","SampleName")
Samples_AllData=Add2Table(Samples_AllData,Samples_Processed,"SampleName","SampleName")
write.table(Samples_AllData,file = paste(ShortFile,"/",ShortFile,"_AllData.txt",sep=""),quote = F,row.names = F,col.names = T,sep="\t")
rm(FA)

#Unique sequences count and reads
#obistat -c count Miseq4_M12sSoil_Samples/Combined_Uniq.fasta > Data_Stats/Miseq4_M12sSoil_Samples_Combined_Uniq.count
seq_counts <- read.table(paste('Data_Stats/',ShortFile,'_Combined_Uniq.count',sep=""), h=T)
names(seq_counts) <- c('count','seqs','reads'); nbSeqs = sum(seq_counts$seqs); nbReads = sum(seq_counts$reads);
pdf(paste(ShortFile,"/",ShortFile,"_ReadCount_per_sequence.pdf",sep=""),width = 10,height = 5)
plot(seq_counts[,c('count','seqs')], log='xy', main = paste(ShortFile,":", nbSeqs, "uniq sequences for", nbReads, "reads"))
dev.off()

#Count singletons
print(paste("#Singletons: ",seq_counts$seqs[seq_counts$count==1]," (",round(seq_counts$seqs[seq_counts$count==1]/nbSeqs*100,0),
            "% of sequences)",sep = ""),quote = F)
rm(seq_counts,nbSeqs,nbReads)



##### 3. Filtering by identity
#Histogram of best identity
pdf(paste(ShortFile,"/",ShortFile,"_Histogram of best identity.pdf",sep=""),width = 10,height = 5)
hist(Data_Ecotag$best_identity,main = paste(ShortFile, "  Histogram of best identity (n=",nrow(Data_Ecotag),")",sep=""),breaks = 50)
dev.off()
#How many sequences only match at max 95%
param$min_identity=0.95
length(which(Data_Ecotag$best_identity<param$min_identity))
#Filtering data based on that threshold
Data_Ecotag_Removed=Data_Ecotag[which(Data_Ecotag$best_identity<param$min_identity),]
Data_Ecotag_Removed$Removal="Identity_thrs"
Data_Ecotag=Data_Ecotag[which(Data_Ecotag$best_identity>=param$min_identity),]



##### 4. ObiClean
#update the dataset
Data_Ecotag=CleanH_count(Data_Ecotag,Samples)
Data_Ecotag_Removed=CleanH_count(Data_Ecotag_Removed,Samples)
#filter count below 10
length(which(Data_Ecotag$count<10))
tData_Ecotag_Removed=Data_Ecotag[which(Data_Ecotag$count<10),]; tData_Ecotag_Removed$Removal="Count"; Data_Ecotag_Removed=rbind(Data_Ecotag_Removed,tData_Ecotag_Removed); rm(tData_Ecotag_Removed)
Data_Ecotag=Data_Ecotag[which(Data_Ecotag$count>=10),]
#Removing sequences that are not head or singleton
length(which(Data_Ecotag$obiclean_head=="False"))
Data_Ecotag=Data_Ecotag[which(Data_Ecotag$obiclean_head=="True"),]

#Any sequence must be head and singleton more than internal
Data_Ecotag_Clean=cbind(Data_Ecotag$obiclean_headcount,Data_Ecotag$obiclean_internalcount,Data_Ecotag$obiclean_singletoncount)
Data_Ecotag_Clean=cbind(Data_Ecotag_Clean,Data_Ecotag_Clean[,1]+Data_Ecotag_Clean[,3]-Data_Ecotag_Clean[,2])
length(which(Data_Ecotag_Clean[,4]<1))
tData_Ecotag_Removed=Data_Ecotag[which(Data_Ecotag_Clean[,4]<1),]; tData_Ecotag_Removed$Removal="Obiclean"; 
Data_Ecotag_Removed=rbind(Data_Ecotag_Removed,tData_Ecotag_Removed); rm(tData_Ecotag_Removed)
Data_Ecotag=Data_Ecotag[which(Data_Ecotag_Clean[,4]>=1),]





##### 5. Group data by reference sequence and calculate relative read abundance
RefSeq_Data=Group_Ecotag_by_RefSeq(Data_Ecotag,Samples)
RefSeq_Data_Removed=Group_Ecotag_by_RefSeq(Data_Ecotag_Removed,Samples)
write.table(RefSeq_Data_Removed,file = paste(substr(EcotagFile,1,nchar(EcotagFile)-4),"_RefSeq_Data_Removed.txt",sep=""),
            col.names = T,row.names = F,quote = F,sep = "\t")
write.table(Data_Ecotag_Removed,file = paste(substr(EcotagFile,1,nchar(EcotagFile)-4),"_Data_Ecotag_Removed.txt",sep=""),
            col.names = T,row.names = F,quote = F,sep = "\t")
round(nrow(RefSeq_Data)/nrow(Data_Ecotag),digits = 3)
#Species to genus correction
param$sp_identity=0.99
RefSeq_Data=sp2genus_correction(RefSeq_Data,param)



##### 6. Remove reference sequences that did not obtain a given taxonomic rank
UniqueRank=UniqueCount(RefSeq_Data$rank)
write.table(UniqueRank, "clipboard", sep="\t", row.names=FALSE)
RefSeq_Data=RefSeq_Data[which(RefSeq_Data$rank!="order"),]
RefSeq_Data=RefSeq_Data[which(RefSeq_Data$rank!="suborder"),]

#Save raw list of MOTUs
write.table(RefSeq_Data,file = paste(substr(EcotagFile,1,nchar(EcotagFile)-4),"_Data_Ecotag_RefSeq.txt",sep=""),
            col.names = T,row.names = F,quote = F,sep = "\t")



##### 6b. Taxonomic correction
colSamples=unique(match(Samples$SampleName,colnames(RefSeq_Data))) #getting column idx of samples
# RefSeq_Data=taxa_correction(RefSeq_Data,Seq_ID,action,NewTaxa=NA,colSamples=NA). Action is either Remove, ChangeTaxa or Group
RefSeq_Data_G=RefSeq_Data
RefSeq_Data_G=taxa_correction(RefSeq_Data_G,Seq_ID=c("GQ264618","AY164517","AF407087","EF153719","KR732817","JX516068","KM078763","AF447214"),action="Remove") #Removing birds
RefSeq_Data_G=taxa_correction(RefSeq_Data_G,Seq_ID=c("AB048589"),action="Remove") #Removing sequences from unclear origin
RefSeq_Data_G=taxa_correction(RefSeq_Data_G,Seq_ID=c("JN632672"),action="Remove") #Removing sequences that should have reach species level
RefSeq_Data_G=taxa_correction(RefSeq_Data_G,Seq_ID=c("KT368731","AP013078"),action="Group",colSamples=colSamples) #Grouping horse sequences
RefSeq_Data_G=taxa_correction(RefSeq_Data_G,Seq_ID=c("AY012150","KM224285"),action="Remove") #Removing sequences that should have reach species level, unexpected felinae
RefSeq_Data_G=taxa_correction(RefSeq_Data_G,Seq_ID=c("U59174","U67289"),action="Group",colSamples=colSamples) #Grouping scurius niger sequences
RefSeq_Data_G=taxa_correction(RefSeq_Data_G,Seq_ID=c("AC008434","AC197159","D38116"),action="Remove") #Removing hominidae sequences not as species level
RefSeq_Data_G=taxa_correction(RefSeq_Data_G,Seq_ID=c("AB196722","AC002087","AC021914","AC062033","AC074322","AC090204"),action="Group",colSamples=colSamples) #Grouping human sequences
RefSeq_Data_G=taxa_correction(RefSeq_Data_G,Seq_ID=c("JF261174","KY707300","KY707303","KY707310"),action="Remove") #Removing likely PCR errors in Cricetidae
RefSeq_Data_G=taxa_correction(RefSeq_Data_G,Seq_ID=c("AJ005780"),action="Remove") #Removing second rattus sequence. Unlikely to have both in JRBP



##### 7. Filtering based on negatives. OTU is removed if it has a maximal average relative abundance in negative controls. 

##Calculate relative read abundance
RefSeq_DataRA=RefSeq_Data_G;
Samples$count=NA; Samples$countPC=NA
for (is in 1:nrow(Samples)){
  Samples$count[is]=sum(RefSeq_Data_G[,match(Samples$SampleName[is],colnames(RefSeq_Data_G))])
  RefSeq_DataRA[,match(Samples$SampleName[is],colnames(RefSeq_Data_G))]=round(RefSeq_Data_G[,match(Samples$SampleName[is],colnames(RefSeq_Data_G))]/Samples$count[is],digits = 4)
}
RefSeq_DataRA$count=round(RefSeq_DataRA$count/sum(RefSeq_DataRA$count),digits = 4)
RefSeq_DataRA$nb_match=round(RefSeq_DataRA$nb_match/sum(RefSeq_DataRA$nb_match),digits = 4)
Samples$countPC=round(Samples$count/sum(Samples$count),digits = 3)

##Output sequences found in negative extractions and PCR
NEG_Samples=Samples[which(Samples$Site=="NEG"),]
colNEG=match(NEG_Samples$SampleName,colnames(RefSeq_Data_G)) #getting column idx of negative samples
RefSeq_Data_NEG=RefSeq_Data_G[,c(1:which(colnames(RefSeq_Data_G)=="sequence"),colNEG)]
colNEG2=match(NEG_Samples$SampleName,colnames(RefSeq_Data_NEG)) #getting column idx of negative extractions
RefSeq_Data_NEG$countN=rowSums(RefSeq_Data_NEG[,colNEG2]) #Adding count for negatives
#RRA of each species in each negative extraction
for(is in 1:nrow(NEG_Samples)){
  RefSeq_Data_NEG$newvar=round(RefSeq_Data_NEG[,which(colnames(RefSeq_Data_NEG)==NEG_Samples$SampleName[is])]/RefSeq_Data_G$count,digits = 7)
  colnames(RefSeq_Data_NEG)[ncol(RefSeq_Data_NEG)]=paste(NEG_Samples$SampleName[is],"_RRA",sep="")
}
RefSeq_Data_NEG=RefSeq_Data_NEG[RefSeq_Data_NEG$countN>0,] #Keeping only detected MOTUs in negative extractions 
write.table(RefSeq_Data_NEG, "clipboard", sep="\t", row.names=FALSE)

#Calculate average relative abundance in negatives
RefSeq_DataRA$AVG_RA_NEG=rowSums(RefSeq_Data_G[,colNEG])/sum(RefSeq_Data_G[,colNEG])
TRUE_Samples=Samples[which(Samples$Site!="NEG"),]
colSamples=match(TRUE_Samples$SampleName,colnames(RefSeq_Data_G)) #getting column idx of samples
RefSeq_DataRA$AVG_RA_TRUE=rowSums(RefSeq_Data_G[,colSamples])/sum(RefSeq_Data_G[,colSamples])
RefSeq_DataRA$Contaminant=RefSeq_DataRA$AVG_RA_NEG>RefSeq_DataRA$AVG_RA_TRUE
#Output reference sequences considered as contaminant
RefSeq_Contaminants=RefSeq_Data_G[RefSeq_DataRA$Contaminant==T,]
write.table(RefSeq_Contaminants[,1:8], "clipboard", sep="\t", row.names=FALSE)

## IF Filter out reference sequences considered as contaminant and deleting negative samples from the dataset
# RefSeq_MOTUs=RefSeq_Data_G[RefSeq_DataRA$Contaminant==F,-match(NEG_Samples$SampleName,colnames(RefSeq_Data_G))]
# RefSeq_MOTUs$count=rowSums(RefSeq_MOTUs[,match(TRUE_Samples$SampleName,colnames(RefSeq_MOTUs))])

## IF not filter out contaminant
RefSeq_MOTUs_G=RefSeq_Data_G[,-colNEG]
RefSeq_MOTUs_G$count=rowSums(RefSeq_MOTUs_G[,match(TRUE_Samples$SampleName,colnames(RefSeq_MOTUs_G))])
RefSeq_MOTUs=RefSeq_Data[,-colNEG]
RefSeq_MOTUs$count=rowSums(RefSeq_MOTUs[,match(TRUE_Samples$SampleName,colnames(RefSeq_MOTUs))])



##### X. Output table TagJumpF
write.table(RefSeq_MOTUs_G,file = paste(substr(EcotagFile,1,nchar(EcotagFile)-4),"_RefSeq_Data_G.txt",sep=""),
            col.names = T,row.names = F,quote = F,sep = "\t")
RefSeq2Fasta(file_location=paste(substr(EcotagFile,1,nchar(EcotagFile)-4),"_RefSeq_Data_G.txt",sep=""),TRUE_Samples)
write.table(RefSeq_MOTUs,file = paste(substr(EcotagFile,1,nchar(EcotagFile)-4),"_RefSeq_Data.txt",sep=""),
            col.names = T,row.names = F,quote = F,sep = "\t")

##### 9. Filtering for tag jump. Any OTU abundance representing less than X% of the total OTU abundance across samples
# #if by species (like zinger et barbat et eDNA book)
param$tagjumpRA=0.0005 #0.05%
RefSeq_MOTUs_G_TJ=RefSeq_MOTUs_G
RefSeq_MOTUs_G_TJ$TJ_countTJ=RefSeq_MOTUs_G$count*param$tagjumpR
colsamples=match(TRUE_Samples$SampleName,colnames(RefSeq_MOTUs_G))
for (iseq in 1:nrow(RefSeq_MOTUs_G_TJ)){
  for (is in 1:length(colsamples)){
    if(RefSeq_MOTUs_G_TJ[iseq,colsamples[is]]<RefSeq_MOTUs_G_TJ$TJ_countTJ[iseq]){
      RefSeq_MOTUs_G_TJ[iseq,colsamples[is]]=0
    }
  }
}

# #if by sample (like port)
# param$tagjumpRA=0.0006
# RefSeq_MOTUs_G_TJ=RefSeq_MOTUs_G
# TRUE_Samples$countTJ=TRUE_Samples$count*param$tagjumpR
# colsamples=match(TRUE_Samples$SampleName,colnames(RefSeq_MOTUs_G_TJ))
# for (iseq in 1:nrow(RefSeq_MOTUs_G_TJ)){
#   for (is in 1:length(colsamples)){
#     if(RefSeq_MOTUs_G_TJ[iseq,colsamples[is]]<TRUE_Samples$countTJ[is]){
#       RefSeq_MOTUs_G_TJ[iseq,colsamples[is]]=0
#     }
#   }
# }

##### X. Output table TagJumpJ
write.table(RefSeq_MOTUs_G_TJ,file = paste(substr(EcotagFile,1,nchar(EcotagFile)-4),"_RefSeq_Data_G_TJ.txt",sep=""),
            col.names = T,row.names = F,quote = F,sep = "\t")

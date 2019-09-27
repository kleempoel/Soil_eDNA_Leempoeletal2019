FileCheck=function(Samples,param){
  fcdata=as.data.frame(Samples$SampleName)
  colnames(fcdata)=c("SampleName")
  fcdata$MiseqName=Samples$MiseqName
  fcdata$gz_r1=NA
  fcdata$gz_r2=NA
  if(!file.exists("index_files.txt")){
    print("can't find index_files.txt") 
  }
  indexfile=read.table("index_files.txt",stringsAsFactors = F)
  indexfile <- data.frame(lapply(indexfile, function(x) {substr(x,3,length(x))}))
  
  #Check if gzfiles exist for each sample and if all gzfiles are used
  gzfiles=as.data.frame(indexfile$V1[grep(".gz",indexfile$V1)])
  colnames(gzfiles)=c("path")
  gzfiles$sample=NA
  for (is in 1:nrow(Samples)){
    tgz=grep(paste("/",Samples$MiseqName[is],sep=""),gzfiles$path)
    if(length(tgz)>0){
      gzfiles$sample[tgz]=as.character(Samples$MiseqName[is])
      tgz2=gzfiles$path[tgz]
      tgz3=grep("R1",tgz2)
      fcdata$gz_r1[is]=as.character(tgz2[tgz3[1]])
      tgz4=grep("R2",tgz2)
      fcdata$gz_r2[is]=as.character(tgz2[tgz4[1]])
      rm(tgz2,tgz3,tgz4)
    }
    rm(tgz)
  }
  UnusedGZ=gzfiles[which(is.na(gzfiles$sample)==T),]
  for (i in nrow(UnusedGZ)){
    print(paste("Warning: File ",UnusedGZ$path," does not correspond to any samples",sep=""))
  }
  
  fastfiles=as.data.frame(indexfile$V1[grep(".fast",indexfile$V1)])
  colnames(fastfiles)=c("path")
  #Check if gzfiles have been unzipped in the temp folder, as well as paired aligned and joined files
  fcdata$r1unzip=NA
  fcdata$r2unzip=NA
  fcdata$pairedend=NA
  fcdata$pairedali=NA
  for (is in 1:nrow(Samples)){
    samplename=Samples$MiseqName[is]
    if (length(grep(paste('tmp/',samplename, '_R1.fastq',sep=''),fastfiles$path))>0){
      fcdata$r1unzip[is]=as.character(fastfiles$path[grep(paste('tmp/',samplename, '_R1.fastq',sep=''),fastfiles$path)]) }
    if (length(grep(paste('tmp/',samplename, '_R2.fastq',sep=''),fastfiles$path))>0){
      fcdata$r2unzip[is]=as.character(fastfiles$path[grep(paste('tmp/',samplename, '_R2.fastq',sep=''),fastfiles$path)]) }
    if (length(grep(paste('tmp/',samplename,'_PairedEnd.fastq',sep=''),fastfiles$path))>0){
      fcdata$pairedend[is]=as.character(fastfiles$path[grep(paste('tmp/',samplename,'_PairedEnd.fastq',sep=''),fastfiles$path)]) }
    if (length(grep(paste('Data_PairedAlign/',samplename,'_PairedEnd.ali.fastq',sep=''),fastfiles$path))>0){
      fcdata$pairedali[is]=as.character(fastfiles$path[grep(paste('Data_PairedAlign/',samplename,'_PairedEnd.ali.fastq',sep=''),fastfiles$path)]) }
  }
  
  #Check if paired files have been cut and filtered
  fcdata$cut=NA
  fcdata$filtered=NA
  for (is in 1:nrow(Samples)){
    samplename=Samples$SampleName[is]
    if (length(grep(paste("Data_Cut/",samplename,"_PairedEnd.ali.cut.fastq",sep=''),fastfiles$path))>0){
      fcdata$cut[is]=as.character(fastfiles$path[grep(paste("Data_Cut/",samplename,"_PairedEnd.ali.cut.fastq",sep=''),fastfiles$path)]) }
    if (length(grep(paste("Data_Filtered/",samplename, "_Filtered.fastq",sep=''),fastfiles$path))>0){
      fcdata$filtered[is]=as.character(fastfiles$path[grep(paste("Data_Filtered/",samplename, "_Filtered.fastq",sep=''),fastfiles$path)])  }
  }
  return(fcdata)
}

Pair_Align=function(Samples,SamplesFile,datafolder,score_min,overwrite=F){
  #Samples  Dataframe containing only original sample dataframe to process 
  #SamplesFile  Name of file containing list of samples.
  #datafolder Directory Containing all gzip files (fine if files are in subfolders)
  #score_min  Minimum paired alignment score (obitools/illuminapairedend)
  
  if (file.exists("File GZ check.txt")){
    print("Loading existing .gz checking file")
    gzall=read.table("File GZ check.txt",header = T,sep = '\t',stringsAsFactors = F)
  }else{
    print("No existing .gz checking file, creating one")
    gzall <- as.data.frame(list.files(path=datafolder,pattern = "gz",recursive = T))
    colnames(gzall)="FileLocation"
    gzall$SampleName=NA
    gzall$check=0
  }
  dir.create('tmp')
  dir.create('Data_PairedAlign')
  Samples$R1_Location=NA
  Samples$R2_Location=NA
  
  Ns=nrow(Samples)
  print(paste(Ns,' samples to process',sep=""))
  appT=0
  SamplesFile2=gsub('.txt','',SamplesFile)
  for (i in 1:Ns){
    samplename=Samples$MiseqName[i]
    print(paste('Processing sample ', samplename, ' (',i,'/',Ns,')...',sep=''))
    samploc=grep(paste('/',samplename,sep=''),gzall$FileLocation)
    gzsample=gzall$FileLocation[samploc]
    gzall$check[samploc]=gzall$check[samploc]+1
    gzall$SampleName[samploc]=as.character(Samples$SampleName[i])
    Samples$R1_Location[i]=as.character(gzsample[1])
    Samples$R2_Location[i]=as.character(gzsample[2])
    fileConn=paste('PairAlign_',SamplesFile2,'.txt',sep="")
    if (length(samploc)>0){
      if(overwrite==F){
        if (file.exists(paste("Data_PairedAlign/",samplename,"_PairedEnd.ali.fastq",sep=''))){
          print(paste('file Data_PairedAlign/',samplename,'_PairedEnd.ali.fastq already exists, no overwrite',sep=''))
        }else{
          if (appT==0){
            write(paste('gunzip -k ', datafolder, '/', gzsample[1], ' tmp/',samplename, '/', samplename, '_R1.fastq',sep=''),fileConn,append=F)
          }else{
            write(paste('gunzip -k ', datafolder, '/', gzsample[1], ' tmp/',samplename, '/', samplename, '_R1.fastq',sep=''),fileConn,append=T)
          }
          write(paste('gunzip -k ', datafolder, '/', gzsample[2], ' tmp/',samplename, '/', samplename, '_R2.fastq',sep=''),fileConn,append=T)
          write(paste('illuminapairedend --score-min=',score_min,' -r ', 
                      'tmp/',samplename,'/',samplename,'_R1.fastq',
                      ' tmp/',samplename,'/',samplename,'_R2.fastq',
                      ' > ','tmp/',samplename,'/',samplename,'_PairedEnd.fastq',sep=''),fileConn,append=T)
          sink(paste('PairAlign_',SamplesFile2,'.txt',sep=""),append = T)
          cat("obigrep -p 'mode!="); cat('"joined"'); cat(paste("' tmp/",samplename,"/",samplename,"_PairedEnd.fastq > ",
                                                                "Data_PairedAlign/",samplename,"_PairedEnd.ali.fastq",sep=""))
          cat('\n')
          cat('\n')
          sink()
          appT=1
        }
      }else{
        if (appT==0){
          write(paste('gunzip -k ', datafolder, '/', gzsample[1], ' tmp/',samplename, '/', samplename, '_R1.fastq',sep=''),fileConn,append=F)
        }else{
          write(paste('gunzip -k ', datafolder, '/', gzsample[1], ' tmp/',samplename, '/', samplename, '_R1.fastq',sep=''),fileConn,append=T)
        }
        write(paste('gunzip -k ', datafolder, '/', gzsample[2], ' tmp/',samplename, '/', samplename, '_R2.fastq',sep=''),fileConn,append=T)
        write(paste('illuminapairedend --score-min=',score_min,' -r ',
                    'tmp/',samplename,'/',samplename,'_R1.fastq',
                    ' tmp/',samplename,'/',samplename,'_R2.fastq',
                    ' > ','tmp/',samplename,'/',samplename,'_PairedEnd.fastq',sep=''),fileConn,append=T)
        sink(paste('PairAlign_',SamplesFile2,'.txt',sep=""),append = T)
        cat("obigrep -p 'mode!="); cat('"joined"'); cat(paste("' tmp/",samplename,"/",samplename,"_PairedEnd.fastq > ",
                                                              "Data_PairedAlign/",samplename,"_PairedEnd.ali.fastq",sep=""))
        cat('\n')
        cat('\n')
        sink()
        appT=1
      }
    }
    
  }
  print(paste(length(which(gzall$check>0)),"/",nrow(gzall)," of Miseq files in folder ",datafolder," processed"))
  print(paste("Copy paste the content of ",'PairAlign_',SamplesFile2,'.txt',
              " in obitools to execute the script",sep=""))
  return(list(gzall,Samples))
}

Obitools_Filters=function(SamplesFile,Samples,PrimerFOR,PrimerREV,mincount,minlength,pcrerror,qualthrs,StatsOnly=F){
  dir.create('Data_Cut')
  dir.create('Data_Stats')
  Ns=nrow(Samples)
  SamplesFile2=gsub('.txt','',SamplesFile)
  for (i in 1:Ns){
    samplename=Samples$SampleName[i]
    MiseqName=Samples$MiseqName[i]
    if (!is.na(Samples$R1_Location[i])){
      if(!StatsOnly){
        fileConn=paste('Cutadapt_',SamplesFile2,'.txt',sep="")
        tcom=paste('cutadapt -a ',PrimerREV,' -q ',qualthrs, ' -m ',minlength,' -M ',minlength+100,
                   ' -o ',"Data_Cut/",samplename,"_PairedEnd.ali.cutREV.fastq",
                   " Data_PairedAlign/",MiseqName,"_PairedEnd.ali.fastq",
                   " --untrimmed-output Data_Cut/",samplename,"_PairedEnd.trimmedREV.fastq",
                   " > ","Data_Cut/",samplename,'_PairedEnd.ali.cutREV.report.txt',sep="")
        # system(tcom)
        if (i==1){
          write(tcom,fileConn,append=F)
        }else{
          write(tcom,fileConn,append=T)
        }
        write(paste('cutadapt -g ',PrimerFOR,' -o ',"Data_Cut/",samplename,"_PairedEnd.ali.cut.fastq",
                    " Data_Cut/",samplename,"_PairedEnd.ali.cutREV.fastq",
                    " --untrimmed-output Data_Cut/",samplename,"_PairedEnd.trimmedFOR.fastq",
                    " > ","Data_Cut/",samplename,'_PairedEnd.ali.cutFOR.report.txt',sep=""),fileConn,append=T)
        write(paste('obiuniq -m sample'," Data_Cut/",samplename,"_PairedEnd.ali.cut.fastq > Data_Uniq/",
                    samplename, '_uniq.unsorted.fasta', sep="" ),fileConn,append = T)
        write(paste("obigrep -l ",minlength," -p 'count>=", mincount, "' Data_Uniq/", samplename, 
                    '_uniq.unsorted.fasta >  Data_Uniq/',samplename, '_uniq.c',mincount,".l",minlength,'.fasta', sep="" ),
              fileConn,append = T)
        write(paste("obiclean -r ",pcrerror," -H ", 
                    'Data_Uniq/',samplename, '_uniq.c',mincount,".l",minlength,'.fasta',
                    ' > Data_Uniq/',samplename, '_uniq.c',mincount,".l",minlength,'.clean.fasta', sep="" ),
              fileConn,append = T)
        write('\n',fileConn,append = T)
        if(file.exists(paste('tmp/',MiseqName,'/',MiseqName,'_R1.fastq',sep=""))){
          write(paste('obicount ','tmp/',MiseqName,'/',MiseqName,'_R1.fastq > Data_Stats/',
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ','tmp/',MiseqName,'/',MiseqName,'_R2.fastq >> Data_Stats/',
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ','tmp/',MiseqName,'/',MiseqName,"_PairedEnd.fastq >> Data_Stats/",
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ',"Data_PairedAlign/",MiseqName,"_PairedEnd.ali.fastq >> Data_Stats/",
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ',"Data_Cut/",samplename,"_PairedEnd.ali.cut.fastq >> Data_Stats/",
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ',"Data_Uniq/",samplename, "_uniq.unsorted.fasta >> Data_Stats/",
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ','Data_Uniq/',samplename, '_uniq.c',mincount,".l",minlength,".fasta >> Data_Stats/",
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ','Data_Uniq/',samplename, '_uniq.c',mincount,".l",minlength,".clean.fasta >> Data_Stats/",
                      samplename,".txt",sep=""),fileConn,append = T)
          write('\n',fileConn,append = T)
        }
      }else{
        fileConn=paste('Cutadapt_',SamplesFile2,'.txt',sep="")
        if(file.exists(paste('tmp/',MiseqName,'/',MiseqName,'_R1.fastq',sep=""))){
          if (i==1){
            write(paste('obicount ','tmp/',MiseqName,'/',MiseqName,'_R1.fastq > Data_Stats/',
                        samplename,".txt",sep=""),fileConn,append = F)
          }else{
            write(paste('obicount ','tmp/',MiseqName,'/',MiseqName,'_R1.fastq > Data_Stats/',
                        samplename,".txt",sep=""),fileConn,append = T)
          }
          write(paste('obicount ','tmp/',MiseqName,'/',MiseqName,'_R2.fastq >> Data_Stats/',
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ','tmp/',MiseqName,'/',MiseqName,"_PairedEnd.fastq >> Data_Stats/",
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ',"Data_PairedAlign/",MiseqName,"_PairedEnd.ali.fastq >> Data_Stats/",
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ',"Data_Cut/",samplename,"_PairedEnd.ali.cut.fastq >> Data_Stats/",
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ',"Data_Uniq/",samplename, "_uniq.unsorted.fasta >> Data_Stats/",
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ','Data_Uniq/',samplename, '_uniq.c',mincount,".l",minlength,".fasta >> Data_Stats/",
                      samplename,".txt",sep=""),fileConn,append = T)
          write(paste('obicount ','Data_Uniq/',samplename, '_uniq.c',mincount,".l",minlength,".clean.fasta >> Data_Stats/",
                      samplename,".txt",sep=""),fileConn,append = T)
          write('\n',fileConn,append = T)
        }
      }
    } 
    write('\n',fileConn,append = T)
    rm(samplename,MiseqName)
  }
}

Stats_Processing=function(Samples){
  Ns=nrow(Samples)
  
  Stats_Samples=Samples
  Stats_Samples$Nr_R1=NA
  Stats_Samples$Nr_R2=NA
  Stats_Samples$Nr_Paired=NA
  Stats_Samples$Nr_Alignedjoined=NA
  Stats_Samples$Nr_Cut=NA
  Stats_Samples$Nr_Filtered=NA
  Stats_Samples$pc_Nr_Paired=NA
  Stats_Samples$pc_Nr_Alignedjoined=NA
  Stats_Samples$pc_Nr_Cut=NA
  Stats_Samples$pc_Nr_Filtered=NA

  for (i in 1:Ns){
    samplename=Samples$SampleName[i]
    if(file.exists(paste("Data_Stats/",samplename,".txt",sep=""))){
      f <- file(paste("Data_Stats/",samplename,".txt",sep=""), open="rb")
      nlines <- 0L
      while (length(chunk <- readBin(f, "raw", 65536)) > 0) {
        nlines <- nlines + sum(chunk == as.raw(10L))
      }
      close(f)
      if (nlines>0){
        SampleStats=read.table(paste("Data_Stats/",samplename,".txt",sep=""))
        if(nrow(SampleStats)==6){
          Stats_Samples$Nr_R1[i]=SampleStats[1,2]
          Stats_Samples$Nr_R2[i]=SampleStats[2,2]
          Stats_Samples$Nr_Paired[i]=SampleStats[3,2]
          Stats_Samples$Nr_Alignedjoined[i]=SampleStats[4,2]
          Stats_Samples$Nr_Cut[i]=SampleStats[5,2]
          Stats_Samples$Nr_Filtered[i]=SampleStats[6,2]
          Stats_Samples$pc_Nr_Paired[i]=Stats_Samples$Nr_Paired[i]/Stats_Samples$Nr_R1[i]
          Stats_Samples$pc_Nr_Alignedjoined[i]=Stats_Samples$Nr_Alignedjoined[i]/Stats_Samples$Nr_R1[i]
          Stats_Samples$pc_Nr_Cut[i]=Stats_Samples$Nr_Cut[i]/Stats_Samples$Nr_R1[i]
          Stats_Samples$pc_Nr_Filtered[i]=Stats_Samples$Nr_Filtered[i]/Stats_Samples$Nr_R1[i]
        }
        rm(SampleStats)
      }
      rm(f)
    }
  }
  return(Stats_Samples)
}

Add2Table=function(ReferenceTable,Table,fieldID_REF,fieldID_Table){
  #find common rows
  fieldIDx_REF=as.numeric(which(colnames(ReferenceTable)==fieldID_REF))
  fieldIDx_Table=as.numeric(which(colnames(Table)==fieldID_Table))
  CommonRef=match(ReferenceTable[,fieldIDx_REF[1]],Table[,fieldIDx_Table])
  Table_sort=Table[CommonRef,]
  NewTable=cbind(ReferenceTable[,fieldIDx_REF[1]],ReferenceTable,Table_sort)
  colnames(NewTable)[1] <- fieldID_REF
  colremove=vector()
  for (ic1 in 1:ncol(NewTable)){
    for (ic2 in ic1:ncol(NewTable)){
      if(ic1!=ic2){
        if(colnames(NewTable[ic1])==colnames(NewTable[ic2])){
          if (length(levels(NewTable[,ic1]))==length(levels(NewTable[,ic2]))){
            if(length(which((NewTable[,ic1]==NewTable[,ic2])==T))==nrow(NewTable)){
              print(paste("column ",colnames(NewTable[ic1]),"[",ic1,"] and ",colnames(NewTable[ic2]),"[",ic2,
                          "] have same name and content, removing redundant column",sep=""),quote = F)
              colremove=c(colremove,ic2)
            }
          }
        }
      }
    }
  }
  colremove=(unique(colremove))
  NewTable2 <- NewTable[ -colremove ]
  NewTable2 <- subset( NewTable, select = -c(colremove) )
  return(NewTable2)
}


#Merge samples
Sequence_Occurence=function(Samples,mincount,minlength){
  library("bio3d")
  tFasta=read.fasta(paste('Data_Uniq/',Samples$SampleName[1], '_uniq.c',mincount,".l",minlength,'.clean.fasta', sep="" ))
  tSeq=matrix()
  for (iseq in 1:nrow(tFasta$ali)){
    tSeq[iseq]=gsub(", ","",toString(tFasta$ali[iseq,]))
  }
  Seq=as.matrix(tSeq)
  for (i in 2:nrow(Samples)){
    samplename=Samples$SampleName[i]
    if(file.exists(paste('Data_Uniq/',samplename, '_uniq.c',mincount,".l",minlength,'.clean.fasta', sep="" ))){
      rm(tFasta)
      tFasta=read.fasta(paste('Data_Uniq/',samplename, '_uniq.c',mincount,".l",minlength,'.clean.fasta', sep="" ))
      tSeq=matrix()
      for (iseq in 1:nrow(tFasta$ali)){
        tSeq[iseq]=gsub(", ","",toString(tFasta$ali[iseq,]))
      }
      Seq=rbind(Seq,as.matrix(tSeq))
    }
  }
  Seq=gsub("-","",Seq)
  SeqU=unique(Seq)
  SeqU_samples=as.data.frame(SeqU)
  colnames(SeqU_samples)=c("Sequence")
  SeqU_samplesM=matrix(0,nrow = nrow(SeqU_samples),ncol = nrow(Samples))
  for (i in 1:nrow(Samples)){
    samplename=Samples$SampleName[i]
    if(file.exists(paste('Data_Uniq/',samplename, '_uniq.c',mincount,".l",minlength,'.clean.fasta', sep="" ))){
      tFasta=read.fasta(paste('Data_Uniq/',samplename, '_uniq.c',mincount,".l",minlength,'.clean.fasta', sep="" ))
      tSeq=matrix()
      for (iseq in 1:nrow(tFasta$ali)){
        tSeq[iseq]=gsub(", ","",toString(tFasta$ali[iseq,]))
      }
      tSeq=gsub("-","",tSeq)
      # tcount=grep("count", readLines(paste('Data_Uniq/',samplename, '_uniq.c',mincount,".l",minlength,'.clean.fasta', sep="" )), value = TRUE)
      for (iseq in 1:nrow(tFasta$ali)){
        # t2=gregexpr(pattern ="obiclean_count",tcount[iseq])
        SeqU_samplesM[which(duplicated(rbind(SeqU, tSeq[iseq]), fromLast = TRUE)),i]=1
        # as.numeric(gsub(".*obiclean_count=\\{'XXX': s*|\\}; sminR.*", "", tcount[iseq]))
      }
    }
  }
  SeqU_samples$Seq_Length=NA
  for (iseq in 1:nrow(SeqU_samples)){
    SeqU_samples$Seq_Length[iseq]=nchar(as.character(SeqU_samples$Sequence[iseq]))
  }
  SeqU_samples$Total=rowSums(SeqU_samplesM)
  SeqU_samples=cbind(SeqU_samples, SeqU_samplesM) 
  colnames(SeqU_samples)=c("Sequence","Seq_Length","Total",as.character(Samples$SampleName))
  return(SeqU_samples)
}

GetCommands_Obitools=function(Samples_Processed,SamplesFile,param){
  fileConn=paste('OBITools_Commands_',substr(SamplesFile,1,nchar(SamplesFile)-4),"_",as.character(format(Sys.time(), "%Y-%b-%d")),'.txt',sep="")
  write("#Commands to process missing files",fileConn,append=F)
  write("mkdir tmp",fileConn,append=F)
  write("mkdir Data_PairedAlign",fileConn,append=T)
  write("mkdir Data_Cut",fileConn,append=T)
  write("mkdir Data_Filtered",fileConn,append=T)
  for (i in 1:nrow(Samples_Processed)){
    #gunzip r1
    if(is.na(Samples_Processed$r1unzip[i])){
      if(!is.na(Samples_Processed$gz_r1[i])){
        write("\n",fileConn,append = T)
        comout=GetObitoolsCommand(SampleName=Samples_Processed$SampleName[i],MiseqName=Samples_Processed$MiseqName[i],filepath=Samples_Processed$gz_r1[i],obicom = "gunzip_r1")
        write(comout,fileConn,append = T)
        rm(comout)
      }else{
        write(paste("#r1 gz file for ",Samples_Processed$MiseqName[i]," does not exist",sep=""),fileConn,append = T)
      }
    }
    #gunzip r2
    if(is.na(Samples_Processed$r2unzip[i])){
      if(!is.na(Samples_Processed$gz_r2[i])){
        comout=GetObitoolsCommand(SampleName=Samples_Processed$SampleName[i],MiseqName=Samples_Processed$MiseqName[i],filepath=Samples_Processed$gz_r2[i],obicom = "gunzip_r2")
        write(comout,fileConn,append = T)
        rm(comout)
      }else{
        write(paste("#r2 gz file for ",Samples_Processed$MiseqName[i]," does not exist",sep=""),fileConn,append = T)
      }
    }
    #paired end
    if(is.na(Samples_Processed$pairedend[i])){
      if(!is.na(Samples_Processed$gz_r1[i])){
        comout=GetObitoolsCommand(SampleName=Samples_Processed$SampleName[i],MiseqName=Samples_Processed$MiseqName[i],param,obicom = "pairedend")
        write(comout,fileConn,append = T)
        rm(comout)
      }else{
        write(paste("#r1 gz file for ",Samples_Processed$MiseqName[i]," does not exist, will not write paired-end command",sep=""),fileConn,append = T)
      }
    }
    #paired aligned
    if(is.na(Samples_Processed$pairedali[i])){
      if(!is.na(Samples_Processed$gz_r1[i])){
        comout=GetObitoolsCommand(SampleName=Samples_Processed$SampleName[i],MiseqName=Samples_Processed$MiseqName[i],obicom = "pairedali")
        # write(comout,fileConn,append = T)
        # rm(comout)
        sink(fileConn,append = T)
        cat(comout)
        cat("\n")
        sink()
        rm(comout)
      }else{
        write(paste("#r2 gz file for ",Samples_Processed$MiseqName[i]," does not exist",sep=""),fileConn,append = T)
      }
    }
    #cutadapt
    if(is.na(Samples_Processed$cut[i])){
      if(!is.na(Samples_Processed$gz_r1[i])){
        comout=GetObitoolsCommand(SampleName=Samples_Processed$SampleName[i],MiseqName=Samples_Processed$MiseqName[i],param,obicom = "cutREV")
        write(comout,fileConn,append = T)
        rm(comout)
        comout=GetObitoolsCommand(SampleName=Samples_Processed$SampleName[i],MiseqName=Samples_Processed$MiseqName[i],param,obicom = "cutFOR")
        write(comout,fileConn,append = T)
        rm(comout)
      }else{
        write(paste("#r1 gz file for ",Samples_Processed$MiseqName[i]," does not exist, will not write cutadapt",sep=""),fileConn,append = T)
      }
    }
    
    #Filter out sequences which length is out of boundn...
    if(is.na(Samples_Processed$filtered[i])){
      if(!is.na(Samples_Processed$gz_r1[i])){
        comout=GetObitoolsCommand(SampleName=Samples_Processed$SampleName[i],MiseqName=Samples_Processed$MiseqName[i],param,obicom = "obigreplength")
        write(comout,fileConn,append = T)
        rm(comout)
      }else{
        write(paste("#r1 gz file for ",Samples_Processed$MiseqName[i]," does not exist, will not write obigreplength command",sep=""),fileConn,append = T)
      }
    }
    
  }
  write('\n',fileConn,append = T)
  write('\n',fileConn,append = T)
  write(paste("mkdir ", substr(SamplesFile,1,nchar(SamplesFile)-4),sep=""),fileConn,append=T)
  for (i in 1:nrow(Samples)){
    if (i==1){
      write(paste('obiannotate --set-tag=sample:',Samples$SampleName[i],' Data_Filtered/',Samples$SampleName[i],'_Filtered.fastq > ',
                  substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_annot.fastq',sep=''),fileConn,append = T)
    }else{
      write(paste('obiannotate --set-tag=sample:',Samples$SampleName[i],' Data_Filtered/',Samples$SampleName[i],'_Filtered.fastq >> ',
                  substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_annot.fastq',sep=''),fileConn,append = T)
    }
  }
  write('\n',fileConn,append = T)
  write(paste('obiuniq -m sample ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_annot.fastq',
              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.fasta',sep=""),fileConn,append = T)
  write(paste('obigrep -p ',"'count>", param$mincount, "' ", 
              substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.fasta',
              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.c', param$mincount,'.fasta',sep=""),fileConn,append = T)
  write(paste('obisort -r -k count ',
              substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.c', param$mincount,'.fasta',
              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Sorted.c', param$mincount,'.fasta',sep=""),fileConn,append = T)
  write(paste("ecotag -d ",param$ecopcr_database," -R ", param$ecopcr_fasta," ",
              substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Sorted.c', param$mincount,'.fasta',
              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4),'/Combined_Uniq_Sorted_Ecotag.c', param$mincount,'.fasta', sep="" ),fileConn,append = T)
  write(paste('obitab -d -o ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Sorted_Ecotag.c', param$mincount,'.fasta',
              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4),'/Combined_Uniq_Sorted_Ecotag.c', param$mincount,'.txt',sep="" ),fileConn,append = T)
  write(paste('obiclean -r 1 -d 1 -H -s merged_sample ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Sorted_Ecotag.c', param$mincount,'.fasta',
              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4),'/Combined_Uniq_Sorted_Ecotag_CleanH.c', param$mincount,'.fasta',sep="" ),fileConn,append = T)
  write(paste('obitab -d -o ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Sorted_Ecotag_CleanH.c', param$mincount,'.fasta',
              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4),'/Combined_Uniq_Sorted_Ecotag_CleanH.c', param$mincount,'.txt',sep="" ),fileConn,append = T)
  
}

GetObitoolsCommand=function(SampleName,MiseqName,param,filepath=NA,obicom){
  if (obicom=="gunzip_r1"){
    comout=paste('gunzip -k -c ', filepath, ' > tmp/', MiseqName, '_R1.fastq',sep='')
  }
  if (obicom=="gunzip_r2"){
    comout=paste('gunzip -k -c ', filepath, ' > tmp/', MiseqName, '_R2.fastq',sep='')
  }
  if (obicom=="pairedend"){
    comout=paste('illuminapairedend --score-min=',param$pairing_score_min,' -r ', 'tmp/',MiseqName,'_R1.fastq',' tmp/',MiseqName,'_R2.fastq',' > ','tmp/',MiseqName,'_PairedEnd.fastq',sep='')
  }
  if (obicom=="pairedali"){
    t1="obigrep -p 'mode!="
    t2='"joined"'
    t3=paste("' tmp/",MiseqName,"_PairedEnd.fastq > ","Data_PairedAlign/",MiseqName,"_PairedEnd.ali.fastq",sep="")
    comout=paste(t1,t2,t3,sep='')
    # comout=gsub("\\","",comout, fixed=TRUE)
  }
  if (obicom=="cutREV"){
    comout=paste('cutadapt -a ',param$PrimerREV,' -q ',param$qualthrs,
          ' -o ',"Data_Cut/",SampleName,"_PairedEnd.ali.cutREV.fastq",
          " Data_PairedAlign/",MiseqName,"_PairedEnd.ali.fastq",
          " --untrimmed-output Data_Cut/",SampleName,"_PairedEnd.trimmedREV.fastq",
          " > ","Data_Cut/",SampleName,'_PairedEnd.ali.cutREV.report.txt',sep="")
  }
  if (obicom=="cutFOR"){
    comout=paste('cutadapt -g ',param$PrimerFOR,' -o ',"Data_Cut/",SampleName,"_PairedEnd.ali.cut.fastq",
                 " Data_Cut/",SampleName,"_PairedEnd.ali.cutREV.fastq",
                 " --untrimmed-output Data_Cut/",SampleName,"_PairedEnd.trimmedFOR.fastq",
                 " > ","Data_Cut/",SampleName,'_PairedEnd.ali.cutFOR.report.txt',sep="")
  }
  if (obicom=="obiannotate"){
    comout=paste('obiannotate --set-tag=sample:',SampleName,' Data_Uniq/', SampleName, '_uniq.unsorted.fasta > Data_Uniq/',
                 SampleName, '_PairedEnd.ali.cut.annot.fastq',sep='')
  }
  if (obicom=="obiuniq"){
    comout=paste('obiuniq -m sample'," Data_Cut/",SampleName,"_PairedEnd.ali.cut.fastq > Data_Uniq/",
                 SampleName, '_uniq.unsorted.fasta', sep="" )
  }
  if (obicom=="obigrepcount"){
    comout=paste("obigrep -l ",param$minlength," -p 'count>=", param$mincount, "' Data_Uniq/", SampleName, 
                 '_uniq.unsorted.fasta >  Data_Uniq/',SampleName, '_uniq.c',param$mincount,".l",param$minlength,'.fasta', sep="" )
  }
  if (obicom=="obigreplength"){
    comout=paste("obigrep --lmin ",param$minlength," --lmax ",param$maxlength,
                 " Data_Cut/",SampleName,"_PairedEnd.ali.cut.fastq ",
                 "> Data_Filtered/",SampleName,"_Filtered.fastq", sep="" )
  }
  if (obicom=="obiclean"){
    comout=paste("obiclean -r ",param$pcrerror," -H ", 
                 'Data_Uniq/',SampleName, '_uniq.c',param$mincount,".l",param$minlength,'.fasta',
                 ' > Data_Uniq/',SampleName, '_uniq.c',param$mincount,".l",param$minlength,'.clean.fasta', sep="" )
  }
  if (obicom=="ecotag"){
    comout=paste("ecotag -d ",param$ecopcr_database," -R ", param$ecopcr_fasta," Data_Uniq/",SampleName, '_uniq.c',param$mincount,".l",param$minlength,
                 '.clean.fasta', " > Data_Ecotag/",SampleName, '_ecotag_uniq.c',param$mincount,".l", param$minlength,'.clean.fasta', sep="" )
  }
  if (obicom=="ecotab"){
    comout=paste("obitab -d -o Data_Ecotag/",SampleName, '_ecotag_uniq.c',param$mincount,".l", param$minlength,'.clean.fasta',
                 " > Data_Ecotag/",SampleName, '_ecotag_uniq.c',param$mincount,".l", param$minlength,'.clean.txt',sep="" )
  }
  return(comout)
}

Obistats_commands=function(Samples_Processed){
  fileConn=paste('Obistats_commands_',as.character(format(Sys.time(), "%Y-%b-%d")),'.txt',sep="")
  write("mkdir Data_Stats",fileConn,append=F)
  for (i in 1:nrow(Samples_Processed)){
    write(paste('obicount ',Samples_Processed$r1unzip[i],' > Data_Stats/',
              Samples_Processed$SampleName[i],".txt",sep=""),fileConn,append = T)
    write(paste('obicount ',Samples_Processed$r2unzip[i],' >> Data_Stats/',
                Samples_Processed$SampleName[i],".txt",sep=""),fileConn,append = T)
    write(paste('obicount ',Samples_Processed$pairedend[i],' >> Data_Stats/',
                Samples_Processed$SampleName[i],".txt",sep=""),fileConn,append = T)
    write(paste('obicount ',Samples_Processed$pairedali[i],' >> Data_Stats/',
                Samples_Processed$SampleName[i],".txt",sep=""),fileConn,append = T)
    write(paste('obicount ',Samples_Processed$cut[i],' >> Data_Stats/',
                Samples_Processed$SampleName[i],".txt",sep=""),fileConn,append = T)
    write(paste('obicount ',Samples_Processed$filtered[i],' >> Data_Stats/',
                Samples_Processed$SampleName[i],".txt",sep=""),fileConn,append = T)
    write('\n',fileConn,append = T)
  }
  write('\n',fileConn,append = T)
  # for (i in 1:nrow(Samples_Processed)){
  #   write(paste('obistat -c score ',Samples_Processed$pairedend[i],' > Data_Stats/',
  #               Samples_Processed$SampleName[i],"_pairedend.scores",sep=""),fileConn,append = T)
  #   
  # }
  write(paste('obicount ',substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_annot.fastq',
            ' > Data_Stats/Stats_', substr(SamplesFile,1,nchar(SamplesFile)-4), '.txt',sep=""),fileConn,append = T)
  write(paste('obicount ',substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.fasta',
              ' >> Data_Stats/Stats_', substr(SamplesFile,1,nchar(SamplesFile)-4), '.txt',sep=""),fileConn,append = T)
  write(paste('obicount ',substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.c', param$mincount,'.fasta',
              ' >> Data_Stats/Stats_', substr(SamplesFile,1,nchar(SamplesFile)-4), '.txt',sep=""),fileConn,append = T)
  # write(paste('obicount ',substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Clean.l',param$minlength,'.L',param$maxlength, '.c', param$mincount,'.fasta',sep=""),fileConn,append = T)
  
}

Ecotag=function(Samples,param,SamplesFile){
  fileConn=paste('Ecotag_commands_',as.character(format(Sys.time(), "%Y-%b-%d")),'.txt',sep="")
  write(paste("mkdir ", substr(SamplesFile,1,nchar(SamplesFile)-4),sep=""),fileConn,append=F)
  for (i in 1:nrow(Samples)){
    # write(paste('obiannotate --set-tag=sample:',Samples$SampleName[i],' Data_Uniq/', Samples$SampleName[i], '_uniq.c', param$mincount, '.l', param$minlength, '.clean.fasta > Data_Uniq/',
    #             Samples$SampleName[i], '_uniq.c', param$mincount, '.l', param$minlength, '.clean.annot.fasta',sep=''),fileConn,append = T)
    if (i==1){
      write(paste('obiannotate --set-tag=sample:',Samples$SampleName[i],' Data_Filtered/',Samples$SampleName[i],'_Filtered.fastq > ',
                  substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_annot.fastq',sep=''),fileConn,append = T)
    }else{
      write(paste('obiannotate --set-tag=sample:',Samples$SampleName[i],' Data_Filtered/',Samples$SampleName[i],'_Filtered.fastq >> ',
                  substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_annot.fastq',sep=''),fileConn,append = T)
    }
  }
  write('\n',fileConn,append = T)
  # for (i in 1:nrow(Samples)){
  #   if (i==1){
  #     # write(paste('cat Data_Uniq/',Samples$SampleName[i], '_uniq.c', param$mincount, '.l', param$minlength, '.clean.annot.fasta > Data_Uniq/Combined_',
  #               # substr(SamplesFile,1,nchar(SamplesFile)-4),'.c', param$mincount, '.l', 
  #               # param$minlength,'.fasta',sep=''),fileConn,append = T)
  #     write(paste('cat Data_Uniq/',Samples$SampleName[i], '_uniq.unsorted.annot.fasta > Data_Uniq/Combined_',
  #                 substr(SamplesFile,1,nchar(SamplesFile)-4),'.unsorted.fasta',sep=''),fileConn,append = T)
  #     write(paste('cat Data_Uniq/',Samples$SampleName[i], '_uniq.unsorted.annot.fasta > Data_Uniq/Combined_',
  #                 substr(SamplesFile,1,nchar(SamplesFile)-4),'.unsorted.fasta',sep=''),fileConn,append = T)
  #   }else{
  #     # write(paste('cat Data_Uniq/',Samples$SampleName[i], '_uniq.c', param$mincount, '.l', param$minlength, '.clean.annot.fasta >> Data_Uniq/Combined_',
  #                 # substr(SamplesFile,1,nchar(SamplesFile)-4),'.c', param$mincount, '.l', 
  #                 # param$minlength,'.fasta',sep=''),fileConn,append = T)
  #     write(paste('cat Data_Uniq/',Samples$SampleName[i], '_uniq.unsorted.annot.fasta >> Data_Uniq/Combined_',
  #                 substr(SamplesFile,1,nchar(SamplesFile)-4),'.unsorted.fasta',sep=''),fileConn,append = T)
  #   }
  # }
  # write('\n',fileConn,append = T)
  # write(paste('obiuniq -m sample Data_Uniq/Combined_', substr(SamplesFile,1,nchar(SamplesFile)-4),'.c', param$mincount, '.l',param$minlength,'.fasta',
      # ' > Data_Uniq/Combined_',substr(SamplesFile,1,nchar(SamplesFile)-4),'_Uniq.c', param$mincount, '.l',param$minlength,'.fasta',sep=""),fileConn,append = T)
  # write(paste('obigrep --lmin ',param$minlength,' --lmax ',param$maxlength, ' ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_annot.fastq ',
              # ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_lengthcut.l',param$minlength,'.L',param$maxlength,'.fasta',sep=""),fileConn,append = T)
  write(paste('obiuniq -m sample ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_annot.fastq',
              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.fasta',sep=""),fileConn,append = T)
  write(paste('obigrep -p ',"'count>", param$mincount, "' ", 
              substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.fasta',
              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.c', param$mincount,'.fasta',sep=""),fileConn,append = T)
  write(paste('obisort -r -k count ',
              substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.c', param$mincount,'.fasta',
              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Sorted.c', param$mincount,'.fasta',sep=""),fileConn,append = T)
  write(paste("ecotag -d ",param$ecopcr_database," -R ", param$ecopcr_fasta," ",
              substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Sorted.c', param$mincount,'.fasta',
              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4),'/Combined_Uniq_Sorted_Ecotag.c', param$mincount,'.fasta', sep="" ),fileConn,append = T)
  write(paste('obitab -d -o ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Sorted_Ecotag.c', param$mincount,'.fasta',
              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4),'/Combined_Uniq_Sorted_Ecotag.c', param$mincount,'.txt',sep="" ),fileConn,append = T)
  
  # write(paste("obiclean -r ",param$pcrerror," -H ",  
  #             substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.l',param$minlength,'.L',param$maxlength, '.c', param$mincount,'.fasta',
  #             ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Clean.c', param$mincount,'.fasta',sep=""),fileConn,append = T)
  # write(paste('obisort -r -k count ',
  #             substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Clean.l',param$minlength,'.L',param$maxlength, '.c', param$mincount,'.fasta',
  #       ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Clean_Sorted.l',param$minlength,'.L',param$maxlength, '.c', param$mincount,'.fasta',sep=""),fileConn,append = T)
  # write(paste("ecotag -d ",param$ecopcr_database," -R ", param$ecopcr_fasta," ",
  #             substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Clean_Sorted.l',param$minlength,'.L',param$maxlength, '.c', param$mincount,'.fasta',
  #             ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Clean_Sorted_Ecotag.l',param$minlength,'.L',param$maxlength, '.c', param$mincount,'.fasta', sep="" ),fileConn,append = T)
  # write(paste('obitab -d -o ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Clean_Sorted_Ecotag.l',param$minlength,'.L',param$maxlength, '.c', param$mincount,'.fasta',
  #              ' > ', substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Clean_Sorted_Ecotag.l',param$minlength,'.L',param$maxlength, '.c', param$mincount,'.txt',sep="" ),fileConn,append = T)
  # write(paste('obicount ',substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_annot.fastq ',sep=""),fileConn,append = T)
  # write(paste('obicount ',substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_lengthcut.l',param$minlength,'.L',param$maxlength,'.fasta',sep=""),fileConn,append = T)
  # write(paste('obicount ',substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.l',param$minlength,'.L',param$maxlength,'.fasta',sep=""),fileConn,append = T)
  # write(paste('obicount ',substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq.l',param$minlength,'.L',param$maxlength, '.c', param$mincount,'.fasta',sep=""),fileConn,append = T)
  # write(paste('obicount ',substr(SamplesFile,1,nchar(SamplesFile)-4), '/Combined_Uniq_Clean.l',param$minlength,'.L',param$maxlength, '.c', param$mincount,'.fasta',sep=""),fileConn,append = T)
  # 
}

ecotag_fieldchange=function(Data_Ecotag){
  ecotag_colnames1=colnames(Data_Ecotag)
  pos_bi=grep("best_identity",ecotag_colnames1)
  db_name=sub("best_identity","",ecotag_colnames1[pos_bi])
  ecotag_colnames2=sub(db_name,"",ecotag_colnames1)
  pos_samples=grep("sample.",ecotag_colnames1)
  ecotag_colnames2[pos_samples]=gsub(x = ecotag_colnames2[pos_samples], pattern = "\\.", replacement = "-") 
  ecotag_colnames2[pos_samples]=gsub(x = ecotag_colnames2[pos_samples], pattern = "sample-", replacement = "")
  colnames(Data_Ecotag)=ecotag_colnames2
  return(Data_Ecotag)
}

Group_Ecotag_by_RefSeq=function(Data_Ecotag,Samples){
  RefSeq_Data=data.frame(RefSeq=unique(Data_Ecotag$best_match))
  #For each reference sequence...
  RefSeq_Data$rank=NA; RefSeq_Data$scientific_name=NA; RefSeq_Data$genus_name=NA; RefSeq_Data$family_name=NA; 
  RefSeq_Data$nb_match=NA; RefSeq_Data$count=NA; RefSeq_Data$max_identity=NA; RefSeq_Data$sequence=NA;
  print("Grouping reference sequences...",quote = F)
  pb <- txtProgressBar(min = 0, max = nrow(RefSeq_Data), style = 3)
  for (iseq in 1:nrow(RefSeq_Data)){
    t_data=Data_Ecotag[which(Data_Ecotag$best_match==RefSeq_Data$RefSeq[iseq]),]
    #Get scientific name, rank, genus, family
    RefSeq_Data$scientific_name[iseq]=as.character(t_data$scientific_name[1])
    RefSeq_Data$rank[iseq]=as.character(t_data$rank[1])
    RefSeq_Data$genus_name[iseq]=as.character(t_data$genus_name[1])
    RefSeq_Data$family_name[iseq]=as.character(t_data$family_name[1])
    RefSeq_Data$order_name[iseq]=as.character(t_data$order_name[1])
    #Count number of unique sequences
    RefSeq_Data$nb_match[iseq]=nrow(t_data)
    #Count total number of reads
    RefSeq_Data$count[iseq]=sum(t_data$count)
    #max best identity
    RefSeq_Data$max_identity[iseq]=max(t_data$best_identity)
    RefSeq_Data$sequence[iseq]=as.character(t_data$sequence[1])
    rm(t_data)
    setTxtProgressBar(pb, iseq)
  }
  close(pb)
  
  #Add count for each sample
  print("Updating count per sample for each reference sequence...",quote = F)
  pb <- txtProgressBar(min = 0, max = nrow(Samples), style = 3)
  for (is in 1:nrow(Samples)){
    RefSeq_Data$new_sample=NA
    tdata=Data_Ecotag[,match(c("best_match",as.character(Samples$SampleName[is])),colnames(Data_Ecotag))]
    for (iseq in 1:nrow(RefSeq_Data)){
      RefSeq_Data$new_sample[iseq]=sum(tdata[which(tdata$best_match==RefSeq_Data$RefSeq[iseq]),2])
    }
    names(RefSeq_Data)[names(RefSeq_Data) == 'new_sample'] <- as.character(Samples$SampleName[is])
    rm(tdata)
    setTxtProgressBar(pb, is)
  }
  close(pb)
  return(RefSeq_Data)
}

UniqueCount=function(data_field){
  data_unique=data.frame(values=unique(data_field))
  data_unique$count=NA
  for (i in 1:nrow(data_unique)){
    data_unique$count[i]=length(which(data_field==data_unique$values[i]))
  }
  return(data_unique)
}

CleanH_count=function(Data_Ecotag,Samples){
  #find and delete obiclean cluster columns
  oc_clusters=grep("obiclean_cluster",colnames(Data_Ecotag))
  Data_Ecotag=Data_Ecotag[,-oc_clusters]
  
  #find cleancount columns and update count
  oc_count=grep("obiclean_count",colnames(Data_Ecotag))
  oc_count_names=gsub(x = colnames(Data_Ecotag)[oc_count], pattern = "\\.", replacement = "-")
  oc_count_names=gsub(x = oc_count_names, pattern = "obiclean_count-", replacement = "")
  
  for (is in 1:length(oc_count_names)){
    col_count=which(colnames(Data_Ecotag)==oc_count_names[is])
    clean_count=Data_Ecotag[,oc_count[is]]
    clean_count[ is.na(clean_count)] <- 0
    Data_Ecotag[,col_count]=clean_count
  }
  Data_Ecotag=Data_Ecotag[,-oc_count]
  colsamples=match(Samples$SampleName,colnames(Data_Ecotag))
  Data_Ecotag$count=rowSums(Data_Ecotag[,colsamples])
  return(Data_Ecotag)
}

RefSeq2Fasta=function(file_location,TRUE_Samples){
  #load data
  RefSeq=read.table(file_location,sep = "\t",header = T)
  colnames(RefSeq)=gsub(x = colnames(RefSeq), pattern = "\\.", replacement = "-") 
  #prepare writing
  RefSeq=RefSeq[order(RefSeq$family),]
  fileConn=paste(file_location,'.fasta',sep="")
  #get col sequence
  colseq=which(colnames(RefSeq)=="sequence")
  colsamples=match(TRUE_Samples$SampleName,colnames(RefSeq))
  colmetadata=seq(1,ncol(RefSeq),1); colmetadata=colmetadata[-(c(colseq,colsamples))]
  for (is in 1:nrow(RefSeq)){
    tmetadata=paste(">",RefSeq$RefSeq[is],"_",RefSeq$scientific_name[is],sep="")
    tmetadata=gsub(x = tmetadata, pattern = " ", replacement = "_")
      for (im in 1:length(colmetadata)){
        tm2=paste(colnames(RefSeq)[colmetadata[im]],"=",RefSeq[is,colmetadata[im]],sep="")
        tmetadata=paste(tmetadata,tm2,sep = "; ")
      }
      tms=paste("'",colnames(RefSeq)[colsamples[1]],"': ",RefSeq[is,colsamples[1]],sep="")
      for (im in 2:length(colsamples)){
        tms2=paste("'",colnames(RefSeq)[colsamples[im]],"': ",RefSeq[is,colsamples[im]],sep="")
        tms=paste(tms,tms2,sep = ", ")
        rm(tms2)
      }
      tms=paste("merged_sample={",tms,"}",sep="")
      tmetadata2=paste(tmetadata,tms,"; ")
      if (is==1){
        write(tmetadata2,fileConn,append = F)
        write(as.character(RefSeq$sequence[is]),fileConn,append = T)
      }else{
        write(tmetadata2,fileConn,append = T)
        write(as.character(RefSeq$sequence[is]),fileConn,append = T)
      }
  }
  
}

sp2genus_correction=function(RefSeq_Data,param){
  #Get all sequences ranked at species and subspecies level
  idx_sp=which(RefSeq_Data$rank %in% c("species", "subspecies"))
  #Get all those for which identity is between given thresholds
  idx_for_genus=which(RefSeq_Data$max_identity[idx_sp]>param$min_identity & RefSeq_Data$max_identity[idx_sp]<param$sp_identity)
  RefSeq_Data$rank[idx_sp[idx_for_genus]]="genus"
  RefSeq_Data$scientific_name[idx_sp[idx_for_genus]]=RefSeq_Data$genus_name[idx_sp[idx_for_genus]]
  return(RefSeq_Data)
}

taxa_correction=function(RefSeq_Data,Seq_ID,action,NewTaxa=NA,colSamples=NA){
  # NewTaxa should be expressed as c(rank, scientific_name, genus (optional), family (optional))
  seq_idx=match(Seq_ID,RefSeq_Data$RefSeq)
  if(is.na(seq_idx)!=T){
    if (action=="Remove"){
      print(paste("Sequence",Seq_ID,RefSeq_Data$rank[seq_idx], RefSeq_Data$scientific_name[seq_idx],RefSeq_Data$genus_name[seq_idx],
                  RefSeq_Data$family[seq_idx],"count=",RefSeq_Data$count[seq_idx],"Removed"),quote = F)
      RefSeq_Data=RefSeq_Data[-seq_idx,]
    }
    if(action=="ChangeTaxa"){
      trefseq=RefSeq_Data[seq_idx,]
      RefSeq_Data$rank[seq_idx]=NewTaxa[1]
      RefSeq_Data$scientific_name[seq_idx]=NewTaxa[2]
      if(length(NewTaxa)>2){
        RefSeq_Data$genus_name[seq_idx]=NewTaxa[3]
      }
      if(length(NewTaxa)>3){
        RefSeq_Data$family[seq_idx]=NewTaxa[4]
      }
      print(paste("Sequence",Seq_ID,trefseq$rank, trefseq$scientific_name,trefseq$genus_name,trefseq$family,"count=",RefSeq_Data$count[seq_idx],
                  " > changed to > ",RefSeq_Data$rank[seq_idx],RefSeq_Data$scientific_name[seq_idx],RefSeq_Data$genus_name[seq_idx],RefSeq_Data$family[seq_idx]),quote = F)
    }
    if(action=="Group"){
      #Keep metadata of the refseq with most reads
      #Update total and per sample count
      RefSeq_togroup=RefSeq_Data[seq_idx,]
      RefSeq_Data$count[seq_idx[1]]=sum(RefSeq_togroup$count)
      RefSeq_Data[seq_idx[1],colSamples]=colSums(RefSeq_togroup[,colSamples])
      #Remove rows that were grouped
      RefSeq_Data=RefSeq_Data[-seq_idx[2:length(seq_idx)],]
    }
  }
  return(RefSeq_Data)
}

group_by_name=function(RefSeq_MOTUs,filegroup){
  #load filegroup
  group_data=read.table(filegroup,sep="\t",header = T)
  
  #Add column to keep track of the number of sequences regrouped
  col_seq=which(colnames(RefSeq_MOTUs)=="sequence")
  tdata=RefSeq_MOTUs[,1:(col_seq-1)]
  tdata$nb_refseq=1
  RefSeq_MOTUs=cbind(tdata,RefSeq_MOTUs[,(col_seq-1):ncol(RefSeq_MOTUs)])
  rm(tdata,col_seq)
  #group sequences by scientific name
  col_seq=which(colnames(RefSeq_MOTUs)=="sequence")
  is=1
  while(is < nrow(RefSeq_MOTUs)){
    togroup=which(RefSeq_MOTUs$scientific_name==RefSeq_MOTUs$scientific_name[is])
    if (length(togroup)>1){
      togroup=togroup[which(togroup!=is)] #group only those that are not the main sequence
      RefSeq_MOTUs$nb_refseq[is]=RefSeq_MOTUs$nb_refseq[is]+length(togroup)
      #Keep metadata of the refseq with most reads
      #Update total and per sample count
      RefSeq_togroup=RefSeq_MOTUs[togroup,]
      RefSeq_MOTUs$count[is]=RefSeq_MOTUs$count[is]+sum(RefSeq_togroup$count)
      RefSeq_MOTUs[is,(col_seq+1):ncol(RefSeq_MOTUs)]=RefSeq_MOTUs[is,(col_seq+1):ncol(RefSeq_MOTUs)]+colSums(RefSeq_togroup[,(col_seq+1):ncol(RefSeq_MOTUs)])
      #Remove rows that were grouped
      RefSeq_MOTUs=RefSeq_MOTUs[-togroup,]
      }
    is=is+1
  }
  
  #delete records according to filegroup Remove
  todel=as.vector(group_data$scientific_name[which(group_data$grouping_action=="Remove")])
  if(length(todel)>0){
    RefSeq_MOTUs=RefSeq_MOTUs[-match(todel,RefSeq_MOTUs$scientific_name),]
  }
  
  #group with other sequence if action is needed
  if(length(which(group_data$grouping_action=="Remove"))>0){
    tochange=group_data[-which(group_data$grouping_action=="Remove"),];
  }else{
    tochange=group_data
  }
  tochange=tochange[-which(tochange$grouping_action=="Keep"),]
  for (is in 1:nrow(tochange)){
    togroup=which(RefSeq_MOTUs$scientific_name==tochange$scientific_name[is])
    mainref=which(RefSeq_MOTUs$scientific_name==tochange$grouping_action[is]) #find the main sequence
    RefSeq_MOTUs$nb_refseq[mainref]=RefSeq_MOTUs$nb_refseq[mainref]+length(togroup)
    #Keep metadata of the main refseq
    #Update total and per sample count
    RefSeq_togroup=RefSeq_MOTUs[togroup,]
    RefSeq_MOTUs$count[mainref]=RefSeq_MOTUs$count[mainref]+sum(RefSeq_togroup$count)
    RefSeq_MOTUs[mainref,(col_seq+1):ncol(RefSeq_MOTUs)]=RefSeq_MOTUs[mainref,(col_seq+1):ncol(RefSeq_MOTUs)]+colSums(RefSeq_togroup[,(col_seq+1):ncol(RefSeq_MOTUs)])
    #Remove rows that were grouped
    RefSeq_MOTUs=RefSeq_MOTUs[-togroup,]
  }
  return(RefSeq_MOTUs)
}

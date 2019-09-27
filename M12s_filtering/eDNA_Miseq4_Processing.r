##### Parameters
#BEFORE ANYTHING, the script needs a list of all files in the folder and subfolder. type 	find > index_files.txt 	in the terminal in the current directory and import the file in the current R directory
# See end of script for mandatory install before lauching script

# List of samples to process in tab delimited columns. Must contain MiseqName and SampleName
# SamplesFile='Miseq4_M16s_Samples.txt'
# SamplesFile='Miseq4_M12sSoil_Samples.txt'
# SamplesFile='Miseq4_M12sWater_Samples.txt'
SamplesFile='Miseq4_All_M12s_Samples.txt'
source("../../R_Scripts/R_Obitools.r") #Location of the Obitools R functions

### Obitools & Cutadapt Parameters
param=data.frame(pairing_score_min=40) #Minimum paired alignment score (obitools/illuminapairedend)
param$qualthrs=30 #Sequence quality threshold, used in cutadapt
param$mincount=10 #Singleton filtering Minimum count of same sequences, pre-blast filtering (obitools/obigrep)

### M16s
# param$PrimerFOR='tcactattttgcnacataga' #Forward primer, M16S mammals  16S A&M Rv2 short 16S 52 TCACTATTTTGCNACATAGA 68-71 Rasmussen et al. (2009)
# param$PrimerREV='tagctcgtctggtttcgggg' #Reverse primer, M16S mammals  16S A&M Fv2 69 16S 52 CCCCGAAACCAGACGAGCTA 68-71 Rasmussen et al. (2009)
# param$minlength=20 #Minimum length of sequence, pre-blast filtering (obitools and cutadapt). 20 for M16s, 50 for IN16STK and worms
# param$maxlength=52 #52 for M16s, 128 for IN16STK
# param$ecopcr_database='../Database/Vertebrates/Vertebrates_EMBL_181108' #EcoPCR Database
# param$ecopcr_fasta='../Database/Vertebrates/M16s_v4_clean_uniqID_uniq_length_stats.fasta'

### M12s
param$PrimerFOR='catagtggggtatctaatcccagtttg' #Forward primer Mimammal
param$PrimerREV='gctggcacgaaatttaccaaccc' #Reverse primer Mimammal
param$minlength=150 #Minimum length of sequence (obitools and cutadapt)
param$maxlength=192 #Maximum length of sequence (obitools and cutadapt)
param$ecopcr_database='../Database/Vertebrates/Vertebrates_EMBL_181108' #EcoPCR Database
param$ecopcr_fasta='../Database/Vertebrates/M12s_v2_clean_uniqID_uniq_length_stats.fasta'


#####1. Checking which files exist and which commands have been run so far
Samples=read.table(SamplesFile,header = T,sep="\t")
Samples_Processed=FileCheck(Samples,param)
write.table(Samples_Processed,file = paste(substr(SamplesFile,1,nchar(SamplesFile)-4),"_filecheck.txt",sep=""),quote = F,row.names = F,col.names = T,sep="\t")
#Opt. Load check table
Samples_Processed=read.table(paste(substr(SamplesFile,1,nchar(SamplesFile)-4),"_filecheck.txt",sep=""),header = T,sep = "\t")

#####2. Get commands
GetCommands_Obitools(Samples_Processed,SamplesFile,param)



##### Install
# For linux, install Python-dev packages
# sudo apt-get install python-dev
# 
# If problem with runit
# https://askubuntu.com/questions/765565/how-to-fix-processing-with-runit-and-git-daemon-run
# 
# Download obitools from
# https://pypi.python.org/pypi/OBITools
# Install Obitools
# http://metabarcoding.org//obitools/doc/welcome.html#installing-the-obitools
# python get-obitools.py
# 
# Activate/Desactivate Obitools with
# ./obitools
# or 
# obitools
# exit with   exit
# 
# Install cutadapt on Linux
# pip install --user --upgrade cutadapt
# 
#Packages to install in R
# install.packages("R.utils")
# install.packages("bio3d")

### To get the Blast database
# ./bin/makeblastdb -in Mammal_12s.fasta -dbtype 'nucl' -out Mammal_12s_Blastplus

### To get the ecopcr database
# obiconvert --fasta Genbank_M12S.fasta -t ./taxdump --ecopcrdb-output=M12S_last
# ecoPCR -d M12S_last -e 3 -l 150 -L 250 CATAGTGGGGTATCTAATCCCAGTTTG GGGTTGGTAAATTTCGTGCCAGC > v05.ecopcr
# obigrep -d M12S_last --require-rank=species --require-rank=genus --require-rank=family v05.ecopcr > M12S_last.fasta
# obiuniq -d M12S_last M12S_last.fasta > M12S_last_clean_uniq.fasta
# obigrep -d M12S_last --require-rank=family M12S_last_clean_uniq.fasta > M12S_last_clean_uniq_clean.fasta
# obiannotate --uniq-id M12S_last_clean_uniq_clean.fasta > db_v05.fasta
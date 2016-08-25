//macro code for Fig.2

dir="\\\\XXXXXX\\";
baseName="XXX";
OG="_w1Cube 2_s";
TMR="_w2Cube 5_s";
Dapi="_w1Cube 1_s";

//mask with Dapi
for(im=1; im<=5; im++) {
    open(dir+baseName+Dapi+im+".tif");
    setAutoThreshold("Default dark");
    //run("Threshold...");
    run("Convert to Mask");
    run("Set Measurements...", "area mean min display redirect=["+baseName+Dapi+im+".tif] decimal=0");
    
    //measure Dapi signal
    run("Analyze Particles...", "size=300-Infinity circularity=0.20-1.00 show=Outlines display add");
    roiManager("Show All with labels");
    saveAs("tiff", dir+"Drawing of_"+baseName+Dapi+im+".tif");
    roiManager("Show All");
    
    //measure OG signal
    open(dir+baseName+OG+im+".tif");
    run("Set Measurements...", "area mean min display redirect=["+baseName+OG+im+".tif] decimal=0");
    roiManager("Measure");
    
    roiManager("Show All with labels");
    saveAs("tiff", dir+"Drawing of_"+baseName+OG+im+".tif");
    close();
    
    //measure TMR signal
    open(dir+baseName+TMR+im+".tif");
    run("Set Measurements...", "area mean min display redirect=["+baseName+TMR+im+".tif] decimal=0");
    roiManager("Measure");
    
    //Close
    run("Close All");
    roiManager("Delete");
    
}



//macro code for Fig.4
dir="D:\\XXXXX\\";
baseName="XX";

for(im=1; im<=97; im++) {
    open(dir+baseName+".tif");
    print(im);
    Stack.setPosition(2,1,im);
    run("Split Channels");
    selectWindow("C2-"+baseName+".tif");
    run("Set Measurements...", "area mean min display redirect=None decimal=0");
    setAutoThreshold("Default dark");
    run("Convert to Mask", "method=Default background=Default calculate");
    run("Fill Holes", "slice");
    //  Option
    //  run("Watershed", "slice");
    
    run("Analyze Particles...", "size=1000-Infinity show=Outlines add");
    roiManager("Show All with labels");
    
    selectWindow("C3-"+baseName+".tif");
    roiManager("Measure");
    close();
    
    roiManager("Deselect");
    roiManager("Delete");
    
    run("Close All");
}

for(im=1; im<=97; im++) {
    open(dir+baseName+".tif");
    print(im);
    Stack.setPosition(2,1,im);
    run("Split Channels");
    selectWindow("C2-"+baseName+".tif");
    run("Set Measurements...", "area mean min display redirect=None decimal=0");
    setAutoThreshold("Default dark");
    run("Convert to Mask", "method=Default background=Default calculate");
    run("Fill Holes", "slice");
    //  Option
    //  run("Watershed", "slice");
    
    run("Analyze Particles...", "size=1000-Infinity show=Outlines add");
    roiManager("Show All with labels");
    
    selectWindow("C1-"+baseName+".tif");
    roiManager("Measure");
    close();
    
    roiManager("Deselect");
    roiManager("Delete");
    
    run("Close All");
}

for(im=1; im<=97; im++) {
    open(dir+baseName+".tif");
    print(im);
    Stack.setPosition(2,1,im);
    run("Split Channels");
    selectWindow("C2-"+baseName+".tif");
    run("Duplicate...", "duplicate");
    run("Set Measurements...", "area mean min display redirect=None decimal=0");
    setAutoThreshold("Default dark");
    run("Convert to Mask", "method=Default background=Default calculate");
    run("Fill Holes", "slice");
    //  Option
    //  run("Watershed", "slice");
    
    run("Analyze Particles...", "size=1000-Infinity show=Outlines add");
    roiManager("Show All with labels");
    selectWindow("C2-"+baseName+".tif");
    roiManager("Measure");
    close();
    
    roiManager("Deselect");
    roiManager("Delete");
    
    run("Close All");
}



for(im=1; im<=97; im++) {
    open(dir+baseName+".tif");
    print(im);
    Stack.setPosition(2,1,im);
    run("Split Channels");
    
    selectWindow("C2-"+baseName+".tif");
    run("Set Measurements...", "area mean min display redirect=None decimal=0");
    setAutoThreshold("Default dark");
    run("Convert to Mask", "method=Default background=Default calculate");
    run("Fill Holes", "slice");
    //  Option
    //  run("Watershed", "slice");
    
    run("Analyze Particles...", "size=1000-Infinity show=Outlines add");
    roiManager("Show All with labels");
    
    selectWindow("C4-"+baseName+".tif");
    roiManager("Measure");
    close();
    
    roiManager("Deselect");
    roiManager("Delete");
    
    run("Close All");
}

if (isOpen("Log")) {
    selectWindow("Log");
    run("Close");
}



##code for ChIP-Seq

##bwa pair-end alignment
module load bwa/0.7.10/gcc.4.4.7
bwa mem female.hg19.fasta J299_BC7YD7ACXX_Lane5_ACAGTG_Input.1_val_1.fq J299_BC7YD7ACXX_Lane5_ACAGTG_Input.2_val_2.fq > input_fpe.sam
bwa mem female.hg19.fasta J299_TH.merged.1_val_1.fq J299_TH.merged.2_val_2.fq > TH_fpe.sam
bwa mem female.hg19.fasta J299_NOC.merged.1_val_1.fq J299_NOC.merged.2_val_2.fq > NOC_fpe.sam

#sam to bam
module load samtools/1.2/gcc.4.4.7
samtools view -b ~/chipseq/input_fpe.sam > input_fpe.bam
samtools view -b ~/chipseq/TH_fpe.sam > TH_fpe.bam
samtools view -b ~/chipseq/NOC_fpe.sam > NOC_fpe.bam

#sort bam files
module load samtools/1.2/gcc.4.4.7
samtools sort ~/chipseq/input_fpe.bam input_sorted_fpe
samtools sort ~/chipseq/TH_fpe.bam TH_sorted_fpe
samtools sort ~/chipseq/NOC_fpe.bam NOC_sorted_fpe

##Remove duplicate
module load  picard-tools/1.119/java.1.7.0_67
java -Xmx2g -jar /public/apps/picard-tools/1.119/bin/MarkDuplicates.jar INPUT=Ninput_sorted_fpe OUTPUT=input_ND_fpe.bam REMOVE_DUPLICATES=true METRICS_FILE=input_metrics.txt
java -Xmx2g -jar /public/apps/picard-tools/1.119/bin/MarkDuplicates.jar INPUT=TH_sorted_fpe.bam OUTPUT=TH_ND_fpe.bam REMOVE_DUPLICATES=true METRICS_FILE=input_metrics.txt
java -Xmx2g -jar /public/apps/picard-tools/1.119/bin/MarkDuplicates.jar INPUT=NOC_sorted_fpe.bam OUTPUT=NOC_ND_fpe.bam REMOVE_DUPLICATES=true METRICS_FILE=input_metrics.txt

##uniq align, sort
module load samtools/1.2/gcc.4.4.7
samtools view -bq 1 input_ND_fpe.bam > input_uniq_fpe.bam
samtools view -bq 1 TH_ND_fpe.bam > TH_uniq_fpe.bam
samtools view -bq 1 NOC_ND_fpe.bam > NOC_uniq_fpe.bam

##down sampling
module load  picard-tools/1.119/java.1.7.0_67
java -Xmx2g -jar /public/apps/picard-tools/1.119/bin/DownsampleSam.jar INPUT=TH_uniq_fpe.bam OUTPUT=TH_127m.bam PROBABILITY=0.723
java -Xmx2g -jar /public/apps/picard-tools/1.119/bin/DownsampleSam.jar INPUT=NOC_uniq_fpe.bam OUTPUT=NOC_127m.bam PROBABILITY=0.892

##bam to bed
module load bedtools2/2.24.0/gcc.4.4.7
bedtools  bamtobed -i input_uniq_fpe.bam > input_127m.bed -sorted
bedtools  bamtobed -i NOC_127m.bam > NOC_127m.bed -sorted
bedtools  bamtobed -i TH_127m.bam > TH_127m.bed -sorted

##coverage in 1-1000kb window
module load bedtools2/2.24.0/gcc.4.4.7
bedtools makewindows -g ~/chipseq/female.hg19.fasta.fai -w 1000 > ~/chipseq/female.hg19_1kbwin.bed
bedtools makewindows -g ~/chipseq/female.hg19.fasta.fai -w 10000 > ~/chipseq/female.hg19_10kbwin.bed
bedtools makewindows -g ~/chipseq/female.hg19.fasta.fai -w 100000 > ~/chipseq/female.hg19_100kbwin.bed
bedtools makewindows -g ~/chipseq/female.hg19.fasta.fai -w 500000 > ~/chipseq/female.hg19_500kbwin.bed
bedtools makewindows -g ~/chipseq/female.hg19.fasta.fai -w 1000000 > ~/chipseq/female.hg19_1000kbwin.bed

bedtools coverage -a ~/chipseq/female.hg19_1kbwin.bed -b ~/chipseq/input_127m.bed > ~/chipseq2/input127m_1kb.bed -sorted
bedtools coverage -a ~/chipseq/female.hg19_1kbwin.bed -b ~/chipseq/NOC_127m.bed > ~/chipseq2/NOC127m_1kb.bed -sorted
bedtools coverage -a ~/chipseq/female.hg19_1kbwin.bed -b ~/chipseq/TH_127m.bed > ~/chipseq2/TH127m_1kb.bed -sorted

bedtools coverage -a ~/chipseq/female.hg19_10kbwin.bed -b ~/chipseq/input_127m.bed > ~/chipseq2/input127m_10kb.bed -sorted
bedtools coverage -a ~/chipseq/female.hg19_10kbwin.bed -b ~/chipseq/NOC_127m.bed > ~/chipseq2/NOC127m_10kb.bed -sorted
bedtools coverage -a ~/chipseq/female.hg19_10kbwin.bed -b ~/chipseq/TH_127m.bed > ~/chipseq2/TH127m_10kb.bed -sorted

bedtools coverage -a ~/chipseq/female.hg19_100kbwin.bed -b ~/chipseq/input_127m.bed > ~/chipseq2/input127m_100kb.bed -sorted
bedtools coverage -a ~/chipseq/female.hg19_100kbwin.bed -b ~/chipseq/NOC_127m.bed > ~/chipseq2/NOC127m_100kb.bed -sorted
bedtools coverage -a ~/chipseq/female.hg19_100kbwin.bed -b ~/chipseq/TH_127m.bed > ~/chipseq2/TH127m_100kb.bed -sorted

bedtools coverage -a ~/chipseq/female.hg19_500kbwin.bed -b ~/chipseq/input_127m.bed > ~/chipseq2/input127m_500kb.bed -sorted
bedtools coverage -a ~/chipseq/female.hg19_500kbwin.bed -b ~/chipseq/NOC_127m.bed > ~/chipseq2/NOC127m_500kb.bed -sorted
bedtools coverage -a ~/chipseq/female.hg19_500kbwin.bed -b ~/chipseq/TH_127m.bed > ~/chipseq2/TH127m_500kb.bed -sorted

bedtools coverage -a ~/chipseq/female.hg19_1000kbwin.bed -b ~/chipseq/input_127m.bed > ~/chipseq2/input127m_1000kb.bed -sorted
bedtools coverage -a ~/chipseq/female.hg19_1000kbwin.bed -b ~/chipseq/NOC_127m.bed > ~/chipseq2/NOC127m_1000kb.bed -sorted
bedtools coverage -a ~/chipseq/female.hg19_1000kbwin.bed -b ~/chipseq/TH_127m.bed > ~/chipseq2/TH127m_1000kb.bed -sorted

paste input127m_1kb.bed NOC127m_1kb.bed TH127m_1kb.bed | awk  '{print $1"\t", $2"\t", $3"\t", $4"\t", $11"\t", $18"\t"}' > 127m_1kbwin_reads.txt
paste input127m_10kb.bed NOC127m_10kb.bed TH127m_10kb.bed | awk  '{print $1"\t", $2"\t", $3"\t", $4"\t", $11"\t", $18"\t"}' > 127m_10kbwin_reads.txt
paste input127m_100kb.bed NOC127m_100kb.bed TH127m_100kb.bed | awk  '{print $1"\t", $2"\t", $3"\t", $4"\t", $11"\t", $18"\t"}' > 127m_100kbwin_reads.txt
paste input127m_500kb.bed NOC127m_500kb.bed TH127m_500kb.bed | awk  '{print $1"\t", $2"\t", $3"\t", $4"\t", $11"\t", $18"\t"}' > 127m_500kbwin_reads.txt
paste input127m_1000kb.bed NOC127m_1000kb.bed TH127m_1000kb.bed | awk  '{print $1"\t", $2"\t", $3"\t", $4"\t", $11"\t", $18"\t"}' > 127m_1000kbwin_reads.txt


##------>GO TO R

one <- read.delim("/Volumes/files/Hanae Sato/chipseq2/127m_1kbwin_reads.txt", header=F)
ten <- read.delim("/Volumes/files/Hanae Sato/chipseq2/127m_10kbwin_reads.txt", header=F)
hundred <- read.delim("/Volumes/files/Hanae Sato/chipseq2/127m_100kbwin_reads.txt", header=F)
fivehundred <- read.delim("/Volumes/files/Hanae Sato/chipseq2/127m_500kbwin_reads.txt", header=F)
thousand <- read.delim("/Volumes/files/Hanae Sato/chipseq2/127m_1000kbwin_reads.txt", header=F)

a1=one[which(one[,4]<quantile(one[,4], 0.99)),]
b1=a1[which(a1[,4]>quantile(a1[,4], 0.10)),]

a10=ten[which(ten[,4]<quantile(ten[,4], 0.99)),]
b10=a10[which(a10[,4]>quantile(a10[,4], 0.10)),]

a100=hundred[which(hundred[,4]<quantile(hundred[,4], 0.99)),]
b100=a100[which(a100[,4]>quantile(a100[,4], 0.05)),]

a500=fivehundred[which(fivehundred[,4]<quantile(fivehundred[,4], 0.99)),]
b500=a500[which(a500[,4]>quantile(a500[,4], 0.10)),]

a1000=thousand[which(thousand[,4]<quantile(thousand[,4], 0.99)),]
b1000=a1000[which(a1000[,4]>quantile(a1000[,4], 0.10)),]

dev1=data.frame(b1[,1:3], b1[,5:6]/b1[,4])
dev10=data.frame(b10[,1:3], b10[,5:6]/b10[,4])
dev100=data.frame(b100[,1:3], b100[,5:6]/b100[,4])
dev500=data.frame(b500[,1:3], b500[,5:6]/b500[,4])
dev1000=data.frame(b1000[,1:3], b1000[,5:6]/b1000[,4])

plot(density(dev1[,5]), xlim=c(0,32), ylim=c(0,0.8), lwd=3, col=2, main="1kb normalized")
lines(density(dev1[,4]),lwd=3, col=4)
legend("topright", inset=.05, c("G1", "S-G2"), lty=1,lwd=3, col=c(4,2))
plot(density(dev10[,5]), xlim=c(0,5), ylim=c(0,1.2), lwd=3, col=2, main="10kb normalized")
lines(density(dev10[,4]),lwd=3, col=4)
plot(density(dev100[,5]), xlim=c(0,12), ylim=c(0,1.2), lwd=3, col=2, main="100kb normalized")
lines(density(dev100[,4]),lwd=3, col=4)
plot(density(dev500[,5]), xlim=c(0,5), ylim=c(0,1.4), lwd=3, col=2, main="500kb normalized")
lines(density(dev500[,4]),lwd=3, col=4)
plot(density(dev1000[,5]), xlim=c(0,5), ylim=c(0,1.2),lwd=3, col=2, main="1000kb normalized")
lines(density(dev1000[,4]),lwd=3, col=4)


###make exclude window (0 in input window)
awk '{if($4==0) print $1"\t"$2"\t"$3}' input127m_500kb.bed > exclude500kb.bed

###Eliminate unnessesary coloum
awk '{print $1"\t"$2"\t"$3}' input_127m.bed > RD0.bed


bedtools shuffle -i RD0.bed -g genomehg19.txt -excl exclude500kb.bed > ~/chipseq2/500kb/RD0.bed
bedtools sort -i RD0.bed > RD0s.bed


##Shuffle x100
#make shuffle100kb.sh

#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -N sort
#$ -j y
#$ -cwd
#$ -m e
#$ -l h_vmem=100g
module load bedtools2/2.24.0/gcc.4.4.7
bedtools sort -i ~/chipseq2/female.hg19_500kbwin.bed > ~/chipseq2/female.hg19_500kbwins.bed

#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -N shuffle
#$ -j y
#$ -cwd
#$ -m e
#$ -l h_vmem=100g
module load bedtools2/2.24.0/gcc.4.4.7
echo "shuffle"
bedtools shuffle -i RD0.bed -g genomehg19.txt -excl exclude500kb.bed -seed ${seed} > ~/chipseq2/500kb/RD${seed}.bed
echo "sort"
bedtools sort -i ~/chipseq2/500kb/RD${seed}.bed > ~/chipseq2/500kb/RD${seed}s.bed
echo "cov"
bedtools coverage -a ~/chipseq2/female.hg19_500kbwins.bed -b ~/chipseq2/500kb/RD${seed}s.bed > ~/chipseq2/500kb/RD${seed}_500kb.bed -sorted


for i in {1..100}; do qsub -v seed=$i shuffle500kb.sh; done


module load bedtools2/2.24.0/gcc.4.4.7
bedtools sort -i input127m_500kb.bed > input127m_500kbs.bed
bedtools sort -i NOC127m_500kb.bed > NOC127m_500kbs.bed
bedtools sort -i TH127m_500kb.bed > TH127m_500kbs.bed

###---->go to R(1)
tmp<- read.delim("/Volumes/files/Hanae Sato/chipseq2/input127m_500kbs.bed",sep="", header=FALSE)
tmp1<- read.delim("/Volumes/files/Hanae Sato/chipseq2/NOC127m_500kbs.bed",sep="", header=FALSE)
tmp2<- read.delim("/Volumes/files/Hanae Sato/chipseq2/TH127m_500kbs.bed",sep="", header=FALSE)
Reads = cbind(tmp[,1:4], tmp1[,4], tmp2[,4])
for (i in 1:100)
{
    print(i)
    tmp_RD <- read.delim(paste("/Volumes/files/Hanae Sato/chipseq2/500kb/RD",i,"_500kb.bed",sep=""), header=FALSE)
    tmp_reads = tmp_RD[,4]
    Reads = cbind(Reads,tmp_reads)
}



File "Reads"
V1: chr name
V2: chr start
V3: chr end
V4: Input
V5: G1 IP
V6: S-G2 IP
V7-106:shuffle


#Remove rows contain top 1% and lower than 99%
a500s=Reads[which(Reads[,4]<quantile(Reads[,4], 0.99)),]
b500s=Reads[which(Reads[,4]>quantile(Reads[,4], 0.1)),]

#Normalize each row by Input
dev500s=data.frame(b500s[,1:3], b500s[,5:106]/b500s[,4])


##get the value of pit
library(pastecs)
devG1=dev500s[which(dev500s[,4]>1.2),]
devG1=devG1[which(devG1[,4]<1.5),]
tp=turnpoints(density(devG1[,4])$y)
summary(tp)
Turning points for: density(devG1[,4])$y

nbr observations  : 512
nbr ex-aequos     : 0
nbr turning points: 5 (first point is a peak)
E(p) = 340 Var(p) = 90.7 (theoretical)

point type         proba     info
1   119 peak  0.000000e+00      Inf
2   230  pit 4.739025e-194 642.2095
3   245 peak 9.058816e-220 727.6449
4   368  pit 1.599368e-239 793.2633
5   398 peak 2.044372e-221 733.1145

density(devG1[,4])$x[368]
[1] 1.440705
abline(v=1.440705, lwd=3, lty=3, col=1)

devS=dev500s[which(dev500s[,5]>1.2),]
devS=devS[which(devS[,5]<1.5),]
tp=turnpoints(density(devS[,5])$y)
summary(tp)

Turning points for: density(devS[, 5])$y

nbr observations  : 512
nbr ex-aequos     : 0
nbr turning points: 5 (first point is a peak)
E(p) = 340 Var(p) = 90.7 (theoretical)

point type         proba     info
1   117 peak 7.610484e-275 910.6022
2   177  pit 4.158042e-109 360.0343
3   202 peak 3.575974e-239 792.1025
4   328  pit 2.456206e-290 962.0627
5   385 peak 8.854510e-293 970.1785

density(devS[, 5])$x[328]
[1] 1.407514
abline(v=1.407514, lwd=3, lty=3, col=1)

##Density plot
plot(density(dev500s[,4]), xlim=c(0,5), ylim=c(0,1.5), axes=FALSE, ann=FALSE, col=4, lwd=5)
box()
axis(3, xaxp=c(0,5,5))
axis(2, xaxp=c(0,1.5,4))
abline(v=1.440705, lwd=3, lty=3, col=1)

plot(density(dev500s[,5]), xlim=c(0,5), ylim=c(0,1.5), axes=FALSE, ann=FALSE, col=2, lwd=5)
box()
axis(2, xaxp=c(0,1.5,4))
abline(v=1.407514, lwd=3, lty=3, col=1)

plot(density(dev500s[,6]), xlim=c(0,5), ylim=c(0,1.5), axes=FALSE, ann=FALSE, col="gray20", lwd=5)
for (i in 7:105){
    print(i)
    lines(density(dev500s[,i]), col="gray20", lwd=3)}
box()
axis(2, xaxp=c(0,1.5,4))



#Take top windows
G1tmp=dev500s[which(dev500s[,4]>1.440705),]
G1=data.frame(G1tmp[,1:3],G1tmp[,4])
Stmp=dev500s[which(dev500s[,5]>1.407514),]
S=data.frame(Stmp[,1:3],Stmp[,5])
write.table(G1,"~/Desktop/Shuffle/1.44_G1.bed", sep="\t",row.names=F,col.names=F,quote=F)
write.table(S,"~/Desktop/Shuffle/1.4_S.bed", sep="\t",row.names=F,col.names=F,quote=F)

for (i in 6:105){
    print(i)
    tmp=dev500s[which(dev500s[,i]>1.440705),]
    tmp2=data.frame(tmp[,1:3],tmp[,i])
    write.table(tmp2,paste("~/Desktop/Shuffle/1.44_",i-5,".bed", sep=""), sep="\t",row.names=F,col.names=F,quote=F)
}

for (i in 6:105){
    print(i)
    tmp=dev500s[which(dev500s[,i]>1.407514),]
    tmp2=data.frame(tmp[,1:3],tmp[,i])
    write.table(tmp2,paste("~/Desktop/Shuffle/1.4_",i-5,".bed", sep=""), sep="\t",row.names=F,col.names=F,quote=F)
}


##Get G band data
https://genome.ucsc.edu
Download
Human
Annotation datbase
Cytoband.txt

##------->Go TO R
write.table(cytoBand,"~/Desktop/Shuffle/cytoBand.bed", sep="\t",row.names=F,col.names=F,quote=F)
input500b=data.frame(b500s[,1:4])
write.table(input500b,"~/Desktop/Shuffle/input500b.bed", sep="\t",row.names=F,col.names=F,quote=F)


##------>GO TO TERMINAL(1)
scp -r ~/Desktop/Shuffle hsato@10.254.207.131:~/chipseq2/Shuffle

module load bedtools2/2.24.0/gcc.4.4.7
bedtools intersect -wo -a ~/chipseq2/Shuffle/1.44_G1.bed -b ~/chipseq2/Shuffle/cytoBand.bed > G1band.bed
bedtools intersect -wo -a ~/chipseq2/Shuffle/1.4_S.bed -b ~/chipseq2/Shuffle/cytoBand.bed > Sband.bed
bedtools intersect -wo -a ~/chipseq2/Shuffle/input500b.bed -b ~/chipseq2/Shuffle/cytoBand.bed > b500band.bed

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$9}' G1band.bed > G1band500.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$9}' Sband.bed > Sband500.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$9}' b500band.bed > band500.bed

#make intersect.sh
#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -N intersect
#$ -j y
#$ -cwd
#$ -m e
#$ -l h_vmem=100g
module load bedtools2/2.24.0/gcc.4.4.7
bedtools intersect -wo -a ~/chipseq2/Shuffle/1.44_${seed}.bed -b ~/chipseq2/Shuffle/cytoBand.bed > 1.44_${seed}band.bed
bedtools intersect -wo -a ~/chipseq2/Shuffle/1.4_${seed}.bed -b ~/chipseq2/Shuffle/cytoBand.bed > 1.4_${seed}band.bed


for i in {1..100}; do qsub -v seed=$i intersect.sh; done

------->Go TO R(2)


//data.einstein.yu.edu/users/Hanae Sato/chipseq2/500kb/


G1_band <- read.delim(paste("/Volumes/files/Hanae Sato/chipseq2/Shuffle/G1band500.bed",sep=""), header=FALSE)
band500 <- read.delim(paste("/Volumes/files/Hanae Sato/chipseq2/Shuffle/band500.bed",sep=""), header=FALSE)

band= c(length(which(G1_band[,5]=="gneg"))/length(which(band500[,5]=="gneg"))/978, length(which(G1_band[,5]=="gpos25"))/length(which(band500[,5]=="gpos25"))/978, length(which(G1_band[,5]=="gpos50"))/length(which(band500[,5]=="gpos50"))/978, length(which(G1_band[,5]=="gpos75"))/length(which(band500[,5]=="gpos75"))/978, length(which(G1_band[,5]=="gpos100"))/length(which(band500[,5]=="gpos100"))/978)


S_band <- read.delim(paste("/Volumes/files/Hanae Sato/chipseq2/Shuffle/Sband500.bed",sep=""), header=FALSE)
tmp_band =  c(length(which(S_band[,5]=="gneg"))/length(which(band500[,5]=="gneg"))/831, length(which(S_band[,5]=="gpos25"))/length(which(band500[,5]=="gpos25"))/831, length(which(S_band[,5]=="gpos50"))/length(which(band500[,5]=="gpos50"))/831, length(which(S_band[,5]=="gpos75"))/length(which(band500[,5]=="gpos75"))/831, length(which(S_band[,5]=="gpos100"))/length(which(band500[,5]=="gpos100"))/831)

band = cbind(band,tmp_band)

##get mean number of coloumn in Shuffle 100x files
RD <- read.delim(paste("//data.einstein.yu.edu/users/Hanae Sato/chipseq2/Shuffle/1.44_1.bed",sep=""), header=FALSE)
row=length(RD[,1])
for (i in 2:100){
    print(i)
    tmp_row <- read.delim(paste("//data.einstein.yu.edu/users/Hanae Sato/chipseq2/Shuffle/1.44_",i,".bed",sep=""), header=FALSE)
    tmp1 = length(tmp_row[,1])
    row = cbind(row,tmp1)
}
mean(row)
[1] 595.41


for (i in 1:100){
    print(i)
    tmp_band <- read.delim(paste("/Volumes/files/Hanae Sato/chipseq2/1.44_",i,"band.bed",sep=""), header=FALSE)
    tmp1 = c(length(which(tmp_band[,9]=="gneg"))/length(which(band500[,5]=="gneg"))/595.41, length(which(tmp_band[,9]=="gpos25"))/length(which(band500[,5]=="gpos25"))/595.41, length(which(tmp_band[,9]=="gpos50"))/length(which(band500[,5]=="gpos50"))/595.41, length(which(tmp_band[,9]=="gpos75"))/length(which(band500[,5]=="gpos75"))/595.41, length(which(tmp_band[,9]=="gpos100"))/length(which(band500[,5]=="gpos100"))/595.41)
    band = cbind(band,tmp1)
}


##plot
x=c(0,25,50,75,100)
plot(x, band[,2], pch=16, col=2, ylim=c(0,6e-04), axes=FALSE, ann=FALSE)
box()
axis(3, xaxp=c(0,100,4))
axis(2, xaxp=c(0,0.3,6))
for (i in 3:102){
    print(i)
    points(x,band[,i], pch=16, col="gray20")
}
points(x,band[,1], pch=16, col=4)
legend("topleft", inset=.05, c("G1", "S-G2", "Random"), pch=c(16,16,16), col=c(4,2,"gray20"))

#plot
plot(density(band[5,3:102]), xlim=c(1e-04,6e-04), xlab="Gband100")
polygon(density(band[5,3:102]), col="gray20", border="gray20")
abline(v=band[5,1], lwd=5, col=4)
abline(v=band[5,2], lwd=5, col=2)

plot(density(band[4,3:102]), xlim=c(1e-04,5e-04), xlab="Gband75")
polygon(density(band[4,3:102]), col="gray20", border="gray20")
abline(v=band[4,1], lwd=5, col=4)
abline(v=band[4,2], lwd=5, col=2)

plot(density(band[3,3:102]), xlim=c(1e-05, 2e-04), xlab="Gband50")
polygon(density(band[3,3:102]), col="gray20", border="gray20")
abline(v=band[3,1], lwd=5, col=4)
abline(v=band[3,2], lwd=5, col=2)

plot(density(band[2,3:102]), xlim=c(1e-05, 2e-04), xlab="Gband25")
polygon(density(band[2,3:102]), col="gray20", border="gray20")
abline(v=band[2,1], lwd=5, col=4)
abline(v=band[2,2], lwd=5, col=2)

plot(density(band[1,3:102]), xlim=c(6e-05, 2e-04), xlab="Gband0")
polygon(density(band[1,3:102]), col="gray20", border="gray20")
abline(v=band[1,1], lwd=5, col=4)
abline(v=band[1,2], lwd=5, col=2)


#Get number of overlap
module load bedtools2/2.24.0/gcc.4.4.7
bedtools coverage -a ~/chipseq2/Shuffle/1.44_G1.bed -b ~/chipseq2/Shuffle/1.4_S.bed > ~/chipseq2/G1inS.bed -sorted

bedtools coverage -a ~/chipseq2/Shuffle/1.44_G1.bed -b ~/chipseq2/Shuffle/1.44_1.bed > ~/chipseq2/G1in1.bed -sorted

awk '{if ($5==1) print}' G1inS.bed |wc -l
810
awk '{if ($5==1) print}' G1in1.bed |wc -l
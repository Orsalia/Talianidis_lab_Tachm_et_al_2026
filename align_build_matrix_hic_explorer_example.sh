Local Alignment with BWA
/home/orsalia/bwa/./bwa mem -A1 -B4  -E50 -L0 -t 8 /home/orsalia/Desktop/Bioinformatics/bwa_index/genomemm10.fa AS-894118-LR-66581_R1.fastq.gz 2>>mate_R1.log | samtools view -Shb - > mate_R1.bam
/home/orsalia/bwa/./bwa mem -A1 -B4  -E50 -L0 -t 8 /home/orsalia/Desktop/Bioinformatics/bwa_index/genomemm10.fa AS-894118-LR-66581_R2.fastq.gz 2>>mate_R2.log | samtools view -Shb - > mate_R2.bam

/home/orsalia/bwa/./bwa mem -A1 -B4  -E50 -L0 -t 8 /home/orsalia/Desktop/Bioinformatics/bwa_index/genomemm10.fa AS-937525-LR-68282_R1.fastq.gz 2>>mate_R1replarge.log | samtools view -Shb - > mate_R1replarge.bam
/home/orsalia/bwa/./bwa mem -A1 -B4  -E50 -L0 -t 8 /home/orsalia/Desktop/Bioinformatics/bwa_index/genomemm10.fa AS-937525-LR-68282_R2.fastq.gz 2>>mate_R2replarge.log | samtools view -Shb - > mate_R2replarge.bam

/home/orsalia/bwa/./bwa mem -A1 -B4  -E50 -L0 -t 8 /home/orsalia/Desktop/Bioinformatics/bwa_index/genomemm10.fa AS-937505-LR-68281_R1.fastq.gz 2>>mate_R1rep2.log | samtools view -Shb - > mate_R1rep2.bam
/home/orsalia/bwa/./bwa mem -A1 -B4  -E50 -L0 -t 8 /home/orsalia/Desktop/Bioinformatics/bwa_index/genomemm10.fa AS-937505-LR-68281_R2.fastq.gz 2>>mate_R2rep2.log | samtools view -Shb - > mate_R2rep2.bam

Find restriction enzyme cut points in genome mm10.fa
python3 digest_genome.py -r ^GATC G^ANTC -o restriction_sites.bed /home/orsalia/Desktop/Bioinformatics/mm10/genomemm10.fa

Heavy process Build contact matrixes for a conventional computer
 PATH=~/anaconda3/bin:$PATH
# Create matrix of interactions "dangling sequence GATC ANT" (10kb bin size)
hicBuildMatrix --samFiles mate_R1.bam mate_R2.bam --binSize 10000 --restrictionSequence GATC GANTC --danglingSequence GATC ANT --restrictionCutFile restriction_sites.bed --outFileName AS-894118-LR-66581.h5 --QCfolder AS-894118-LR-66581_QC --threads 8 --inputBufferSize 6000
 
hicBuildMatrix --samFiles mate_R1replarge.bam mate_R2replarge.bam --binSize 10000 --restrictionSequence GATC GANTC --danglingSequence GATC ANT --restrictionCutFile restriction_sites.bed --outFileName AS-937525-LR-68282.h5 --QCfolder AS-937525-LR-68282_QC --threads 8 --inputBufferSize 4000 --skipDuplicationCheck
 
hicBuildMatrix --samFiles mate_R1rep2.bam mate_R2rep2.bam --binSize 10000 --restrictionSequence GATC GANTC --danglingSequence GATC ANT --restrictionCutFile restriction_sites.bed --outFileName AS-937505-LR-68281.h5 --QCfolder AS-937505-LR-68281_QC --threads 8 --inputBufferSize 4000 --skipDuplicationCheck

#!/bin/bash

### Define values ---------------------------------------------------------

#define project folder
IN_DIR=/path/to/data/
#set database folder
db=path/to/db

#set virus for filtering
blast="$VIRUS"_refs

#set compute resources where m=GB t=CPU
MEMORY=2
THREADS=10

#change into project folder
cd $IN_DIR
mkdir 1_reads

#then, copy source AGRF reads into reads folder while renaming them
while read i; do
cp $IN_DIR/source/"$i"_*_R1_001.fastq.gz $IN_DIR/1_reads/"$i"_R1.fastq.gz
cp $IN_DIR/source/"$i"_*_R2_001.fastq.gz $IN_DIR/1_reads/"$i"_R2.fastq.gz
done < $IN_DIR/IDs.list


#quality trim reads with bbduk/37.99
while read i; do
bbduk.sh -Xmx"$MEMORY"g threads="$THREADS" in=$IN_DIR/1_reads/"$i"_R1.fastq.gz in2=$IN_DIR/1_reads/"$i"_R2.fastq.gz out=stdout.fq ref=phix k=31 hdist=1 | bbduk.sh -Xmx"$MEMORY"g threads="$THREADS" interleaved=true in=stdin.fq out=stdout.fq ref=adapters ktrim=r k=17 mink=3 | bbduk.sh -Xmx"$MEMORY"g threads="$THREADS" interleaved=true in=stdin.fq out=$IN_DIR/1_reads/"$i"_trim_1.fastq.gz out2=$IN_DIR/1_reads/"$i"_trim_2.fastq.gz minlen=50 qtrim=rl trimq=20 entropy=0.7
done < $IN_DIR/IDs.list

#normalise read coverage for denovo assembly


#QC check
while read i; do
fastqc -t "$THREADS" $IN_DIR/1_reads/"$i"_*.fastq.gz
done < $IN_DIR/IDs.list
multiqc $IN_DIR/1_reads/*_fastqc.zip --ignore *_norm* -o $IN_DIR/1_reads/multiqc_fastQC_report -n multiqc_fastQC_report
rm $IN_DIR/1_reads/*_fastqc.zip $IN_DIR/1_reads/*_fastqc.html; 

#CLIPPING

VIRUS=SARSCOV2.Eden
while read i; do
cutadapt -j "$THREADS" -a file:$db/primers/"$VIRUS".primers.clips.A.fa -A file:$db/primers/"$VIRUS".primers.clips.A.fa -o $IN_DIR/1_reads/"$i"_trim_temp1_1.fastq.gz -p $IN_DIR/1_reads/"$i"_trim_temp1_2.fastq.gz $IN_DIR/1_reads/"$i"_trim_1.fastq.gz $IN_DIR/1_reads/"$i"_trim_2.fastq.gz
cutadapt -j "$THREADS" -g file:$db/primers/"$VIRUS".primers.clips.B.fa -G file:$db/primers/"$VIRUS".primers.clips.B.fa -o $IN_DIR/1_reads/"$i"_trim_temp2_1.fastq.gz -p $IN_DIR/1_reads/"$i"_trim_temp2_2.fastq.gz $IN_DIR/1_reads/"$i"_trim_temp1_1.fastq.gz $IN_DIR/1_reads/"$i"_trim_temp1_2.fastq.gz
bbduk.sh -Xmx"$MEMORY"g threads="$THREADS" in=$IN_DIR/1_reads/"$i"_trim_temp2_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_temp2_2.fastq.gz out=$IN_DIR/1_reads/"$i"_trim_clip_1.fastq.gz out2=$IN_DIR/1_reads/"$i"_trim_clip_2.fastq.gz minlen=50
rm $IN_DIR/1_reads/"$i"_trim_temp*
done < $IN_DIR/IDs.list; 

#create folder for reference mapping
mkdir $IN_DIR/2_ref_map

#map trimmed read to reference genomes
while read i; do
   
bbmap.sh -Xmx"$MEMORY"g threads="$THREADS" in=$IN_DIR/1_reads/"$i"_trim_clip_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_clip_2.fastq.gz outm=$IN_DIR/2_ref_map/"$i".ref_mapped.bam ref=$db/MN908947.3.fasta;
samtools sort -@ "$THREADS" -o $IN_DIR/2_ref_map/"$i".ref_mapped.sorted.bam $IN_DIR/2_ref_map/"$i".ref_mapped.bam
rm $IN_DIR/2_ref_map/"$i".ref_mapped.bam
done < $IN_DIR/IDs.list;


mkdir $IN_DIR/3_vcf

while read i; do

bcftools mpileup --no-BAQ --count-orphans --min-BQ 0 -f $db/MN908947.3.fasta $IN_DIR/2_ref_map/"$i".ref_mapped.sorted.bam | bcftools call -mA --ploidy 1 - > $IN_DIR/3_vcf/"$i".vcf
done < $IN_DIR/IDs.list;




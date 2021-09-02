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
done < $IN_DIR/IDs.EDEN.list; 

#create folder for reference mapping
mkdir $IN_DIR/2_ref_map

#map trimmed read to reference genomes
while read i; do
   
bbmap.sh -Xmx"$MEMORY"g threads="$THREADS" in=$IN_DIR/1_reads/"$i"_trim_clip_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_clip_2.fastq.gz outm=$IN_DIR/2_ref_map/"$i".ref_mapped.bam ref=$db/MN908947.3.fasta;
samtools sort -@ "$THREADS" -o $IN_DIR/2_ref_map/"$i".ref_mapped.sorted.bam $IN_DIR/2_ref_map/"$i".ref_mapped.bam
rm $IN_DIR/2_ref_map/"$i".ref_mapped.bam
done < $IN_DIR/IDs.EDEN.list;


mkdir $IN_DIR/3_vcf

while read i; do

bcftools mpileup --no-BAQ --count-orphans --min-BQ 0 -f $db/MN908947.3.fasta $IN_DIR/2_ref_map/"$i".ref_mapped.sorted.bam | bcftools call -mA --ploidy 1 - > $IN_DIR/3_vcf/"$i".vcf
done < $IN_DIR/IDs.EDEN.list;


#CLIPPING

VIRUS=SARSCOV2.Midnight
while read i; do
cutadapt -j "$THREADS" -a file:$db/primers/"$VIRUS".primers.clips.A.fa -A file:$db/primers/"$VIRUS".primers.clips.A.fa -o $IN_DIR/1_reads/"$i"_trim_temp1_1.fastq.gz -p $IN_DIR/1_reads/"$i"_trim_temp1_2.fastq.gz $IN_DIR/1_reads/"$i"_trim_1.fastq.gz $IN_DIR/1_reads/"$i"_trim_2.fastq.gz
cutadapt -j "$THREADS" -g file:$db/primers/"$VIRUS".primers.clips.B.fa -G file:$db/primers/"$VIRUS".primers.clips.B.fa -o $IN_DIR/1_reads/"$i"_trim_temp2_1.fastq.gz -p $IN_DIR/1_reads/"$i"_trim_temp2_2.fastq.gz $IN_DIR/1_reads/"$i"_trim_temp1_1.fastq.gz $IN_DIR/1_reads/"$i"_trim_temp1_2.fastq.gz
bbduk.sh -Xmx"$MEMORY"g threads="$THREADS" in=$IN_DIR/1_reads/"$i"_trim_temp2_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_temp2_2.fastq.gz out=$IN_DIR/1_reads/"$i"_trim_clip_1.fastq.gz out2=$IN_DIR/1_reads/"$i"_trim_clip_2.fastq.gz minlen=50
rm $IN_DIR/1_reads/"$i"_trim_temp*
done < $IN_DIR/IDs.MIDNIGHT.list; 

#create folder for reference mapping
mkdir $IN_DIR/2_ref_map

#map trimmed read to reference genomes
while read i; do
   
bbmap.sh -Xmx"$MEMORY"g threads="$THREADS" in=$IN_DIR/1_reads/"$i"_trim_clip_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_clip_2.fastq.gz outm=$IN_DIR/2_ref_map/"$i".ref_mapped.bam ref=$db/MN908947.3.fasta;
samtools sort -@ "$THREADS" -o $IN_DIR/2_ref_map/"$i".ref_mapped.sorted.bam $IN_DIR/2_ref_map/"$i".ref_mapped.bam
rm $IN_DIR/2_ref_map/"$i".ref_mapped.bam
done < $IN_DIR/IDs.MIDNIGHT.list; 


mkdir $IN_DIR/3_vcf

while read i; do

bcftools mpileup --no-BAQ --count-orphans --min-BQ 0 -f $db/MN908947.3.fasta $IN_DIR/2_ref_map/"$i".ref_mapped.sorted.bam | bcftools call -mA --ploidy 1  - > $IN_DIR/3_vcf/"$i".vcf
done < $IN_DIR/IDs.MIDNIGHT.list; 




mkdir $IN_DIR/3_vcf

while read i; do

samtools mpileup -f $db/MN908947.3.fasta $IN_DIR/2_ref_map/"$i".ref_mapped.sorted.bam --no-BAQ --count-orphans --min-BQ 0 | sed -n '500,515p'
done < $IN_DIR/IDs.EDEN.list;

#generate coverage maps
cd $IN_DIR/2_ref_map
while read i; do
qualimap bamqc -bam "$i".ref_mapped.sorted.bam  -outformat PDF -outfile "$i".ref_mapped.coverage.pdf --java-mem-size="$MEMORY"G
done < $IN_DIR/IDs.list
cd $IN_DIR

#cleanup files
while read i; do
mv $IN_DIR/2_ref_map/"$i".ref_mapped.sorted_stats/"$i".ref_mapped.coverage.pdf $IN_DIR/2_ref_map/
rm -r $IN_DIR/2_ref_map/"$i".ref_mapped.sorted_stats
done < $IN_DIR/IDs.list























#create folder for de novo assembly
mkdir $IN_DIR/3_contigs

#denovo assemble trimmed reads with megahit and rename output contigs with library name
while read i; do
megahit -t "$THREADS" -1 $IN_DIR/1_reads/"$i"_trim_norm_1.fastq.gz -2 $IN_DIR/1_reads/"$i"_trim_norm_2.fastq.gz -o $IN_DIR/3_contigs/megahit_"$i" --out-prefix "$i" --min-contig-len 500;
sed "s/>k/>Draft_"$i"_k/g" $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.fa > $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.fa;
done < $IN_DIR/IDs.list

#export denovo assembly contig stats
while read i; do
grep ">" $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.fa | cut -d ">" -f2 | tr ' ' '\t' | sed 's/flag\=//g' | sed 's/multi\=//g' | sed 's/len\=//g' > $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.stats
done < $IN_DIR/IDs.list

#extract high coverage denovo assembly contigs
while read i; do
awk -F$'\t' '{OFS=FS}{if ($3>50) print $1}' $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.stats > $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.hicov.list
seqtk subseq $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.fa $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.hicov.list > $IN_DIR/3_contigs/megahit_"$i"/"$i".contigs.megahit.hicov.fa
done < $IN_DIR/IDs.list

#create folder for blasting
mkdir $IN_DIR/4_filter

#combine denovo contigs into single file
cat $IN_DIR/3_contigs/megahit_*/*.contigs.megahit.hicov.fa > $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.fasta

#blast to local database
export BLASTDB=$db/blast
blastn -query $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.fasta -db "$db"/blast/"$VIRUS"_refs -evalue 1E-10 -num_threads 1 -out $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.blastn_"$VIRUS".txt -outfmt "6 qseqid qlen stitle sstart send pident length evalue sstrand"

#get top virus blast results for each contig
awk -F$'\t' '!seen[$1]++' $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.blastn_"$VIRUS".txt > $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.blastn_"$VIRUS".top.txt

#take column containing contig names
cut -f1 $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.blastn_"$VIRUS".top.txt > $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".list 

#retrieve sequences
seqtk subseq $IN_DIR/4_filter/all_denovo.contigs.megahit.hicov.fasta $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".list > $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".fa

#align to reference
mafft --thread "$THREADS" --reorder --adjustdirection --maxiterate 10 --add $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".fa $db/fasta/"$VIRUS"_refs.fa > $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".ref_aligned.fa

#unalign draft genome alignment
sed '/^>/! s/\-//g' $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".ref_aligned.fa | sed 's/_R_//g' | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' > $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".ref_unaligned.fa

#retrieve draft denovo genomes
seqtk subseq $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".ref_unaligned.fa $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".list > $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".ref_unaligned.reoriented.fa

#create folder for mapping
mkdir $IN_DIR/5_remap

#retrieve individual draft sequences for mapping
while read i; do
grep -A1 ">Draft_${i}_" $IN_DIR/4_filter/"$VIRUS"_draft_genomes."$RUN".ref_unaligned.reoriented.fa > $IN_DIR/5_remap/"$i".draft.fa
done < $IN_DIR/IDs.list

#map trimmed (clipped) reads to draft genomes
while read i; do
    if [ "${CLIP_FLAG}" == "YES" ]; then 
        bbmap.sh -Xmx"$MEMORY"g threads="$THREADS" maxindel=200 minid=0.98 in=$IN_DIR/1_reads/"$i"_trim_clip_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_clip_2.fastq.gz outm=$IN_DIR/5_remap/"$i".remapped.bam ref=$IN_DIR/5_remap/"$i".draft.fa;
else bbmap.sh -Xmx"$MEMORY"g threads="$THREADS" maxindel=200 minid=0.98 in=$IN_DIR/1_reads/"$i"_trim_1.fastq.gz in2=$IN_DIR/1_reads/"$i"_trim_2.fastq.gz outm=$IN_DIR/5_remap/"$i".remapped.bam ref=$IN_DIR/5_remap/"$i".draft.fa
fi;
samtools sort -@ "$THREADS" -o $IN_DIR/5_remap/"$i".remapped.sorted.bam $IN_DIR/5_remap/"$i".remapped.bam
rm $IN_DIR/5_remap/"$i".remapped.bam
done < $IN_DIR/IDs.list

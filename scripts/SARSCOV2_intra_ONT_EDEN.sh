#!/bin/bash

### variables
# set working directory
	wrk=/path/to/output/location
# name of output folder
	prj=name_of_output_folder
# location of demultiplexed fastq ONT data
	fastq_loc=/path/to/fastq
# location of fast5 ONT data
	fast5_loc=/path/to/fast5
# location of sequencing sumary ONT data
	seq_sum_loc=/path/to/summary

# location of primer schemes
	prim=path/to/primer_schemes

# install location of varscan
	instal_varscan=/path/to/varscan/
# install location of ngmlr
	instal_ngmlr=/path/to/ngmlr-0.2.7/
# install location of sniffles
	instal_sniffles=/path/to/Sniffles-master



outp=$wrk/$prj

### Genome assembly of SARS-CoV-2 viruses: medaka & nanopolish

# gather fastq

mkdir $outp/fastq_files
cd $outp/fastq_files
for bc in {01..09};do mkdir barcode${bc};
cd barcode${bc};

artic guppyplex --skip-quality-check --min-length 1800 --max-length 2700 --directory $fastq_loc/barcode${bc} --prefix barcode${bc}_pass;
cd ../
done


# run medaka
mkdir $outp/medaka
cd $outp/medaka

for bc in {01..09};do mkdir barcode${bc};
cd barcode${bc};
artic minion --medaka --normalise 2000 --threads 4 --scheme-directory $prim --read-file $outp/fastq_files/barcode${bc}/barcode${bc}_pass_barcode${bc}.fastq nCoV-2019/Eden barcode${bc}_medaka;
cd ../
done


# run nanopolish
mkdir $outp/nanopolish
cd $outp/nanopolish


for bc in {01..09};do mkdir barcode${bc};
cd barcode${bc};
artic minion --normalise 2000 --threads 4 --scheme-directory $prim --read-file $outp/fastq_files/barcode${bc}/barcode${bc}_pass_barcode${bc}.fastq --sequencing-summary $seq_sum_loc/sequencing_summary.txt --fast5-directory $fast5_loc/fast5_pass/barcode${bc}/ nCoV-2019/Eden barcode${bc}_nanopolish;
cd ../
done
cd ../


### Variant calling of SARS-CoV-2 viruses on nanopolish data


cd $outp/nanopolish

for bc in {01..09};do cd barcode${bc}*;
samtools mpileup -f $prim/nCoV-2019/Eden/nCoV-2019.reference.fasta *_nanopolish.primertrimmed.rg.sorted.bam | java -jar $instal_varscan/VarScan.v2.4.3.jar pileup2snp > TSV_all.tsv
samtools mpileup -f $prim/nCoV-2019/Eden/nCoV-2019.reference.fasta *_nanopolish.primertrimmed.rg.sorted.bam | java -jar $instal_varscan/VarScan.v2.4.3.jar pileup2snp --min-var-freq 0.8 > TSV_80.tsv
cd ../
done

### Strucutural Variant calling of SARS-CoV-2 viruses on nanopolish data



for bc in {01..09};do cd barcode${bc}*;
$instal_ngmlr/ngmlr -t 4 -r $prim/nCoV-2019/Eden/nCoV-2019.reference.fasta -q $outp/fastq_files/barcode${bc}/barcode${bc}_pass_barcode${bc}.fastq -o barcode${bc}.sam -x ont
samtools view -S -b barcode${bc}.sam | samtools sort -o barcode${bc}.sorted.bam
$instal_sniffles/bin/sniffles-core-1.0.12/sniffles --threads 4 --min_support 20 --min_length 10 -m barcode${bc}.sorted.bam -v sniffles_SV_output.vcf
cd ../
done












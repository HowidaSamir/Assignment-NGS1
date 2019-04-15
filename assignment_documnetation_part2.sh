#5.2.HISAT alignment for shuffled samples:

#install Hisat
source activate ngs1
conda install -c bioconda hisat2 

mkdir -p ~mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/hisat_align/hisatIndex &&cd ~/home/ngs-01/nu-ngs01/Assignment_NGS1/hisat_align/hisatIndex

ln -s /home/ngs-01/nu-ngs01/Assignment-NGS1/RawData/ShuffledSamples/SX_1_agr_trimmed/gencode.v29.pc_transcripts.chr22.simplified.fa .

cd ~/workdir/Assignment/hisat_align```
for SAMPLE in 1 2 3 4 5;
    do
        R1=$HOME/workdir/Assignment/Aggressive_trimmed/shuff_SRR8797509_1.part_001.part_00${SAMPLE}.pe.trim.fastq
        R2=$HOME/workdir/Assignment/Aggressive_trimmed/shuff_SRR8797509_2.part_001.part_00${SAMPLE}.pe.trim.fastq
        hisat2 -p 1 -x hisatIndex/chr22_with_ERCC92 --dta --rna-strandness RF -1 $R1 -2 $R2 -S shuff_SRR8797509_part_00${SAMPLE}.sam
done


# Data Assembly 

# Assembly for the files resulted from BWA alignment on unshuffled data

```cd ~/workdir/Assignment/bwa_align```

## Prepare the SAM file for assembly
## install Samtools
## conda install samtools

## convert the SAM file into BAM file 
samtools view -bS SRR8797509_part_001.sam > SRR8797509_part_001.bam
samtools view -bS SRR8797509_part_002.sam > SRR8797509_part_002.bam
samtools view -bS SRR8797509_part_003.sam > SRR8797509_part_003.bam
samtools view -bS SRR8797509_part_004.sam > SRR8797509_part_004.bam
samtools view -bS SRR8797509_part_005.sam > SRR8797509_part_005.bam

## convert the BAM file to a sorted BAM file. 
samtools sort SRR8797509_part_001.bam -o SRR8797509_part_001.sorted.bam
samtools sort SRR8797509_part_002.bam -o SRR8797509_part_002.sorted.bam
samtools sort SRR8797509_part_003.bam -o SRR8797509_part_003.sorted.bam
samtools sort SRR8797509_part_004.bam -o SRR8797509_part_004.sorted.bam
samtools sort SRR8797509_part_005.bam -o SRR8797509_part_005.sorted.bam

## Export some useful statistics report for each sample indvidually
for f in 1 2 3 4 5;
    do
        samtools flagstat SRR8797509_part_00$f.sorted.bam > useful_stat_$f.txt;
done


## install required tool for assembly 
## install stringtie

## Assembly without known annotations
for SAMPLE in 1 2 3 4 5;
    do

        stringtie SRR8797509_part_00${SAMPLE}.sorted.bam --rf -l ref_free_${SAMPLE} -o ref_free_${SAMPLE}.gtf
done

# Assembly with known previous annotations
for SAMPLE in 1 2 3 4 5;
    do
        stringtie SRR8797509_part_00${SAMPLE}.sorted.bam --rf -l ref_sup_${SAMPLE} -G ~/workdir/sample_data/chr22_with_ERCC92.gtf -o ref_sup_${SAMPLE}.gtf 
done

# Assembly for the files resulted from Hisat alignment on shuffled data

```cd ~/workdir/Assignment/hisat_align```

## convert the SAM file into BAM file 
samtools view -bS shuff_SRR8797509_part_001.sam > shuff_SRR8797509_part_001.bam
samtools view -bS shuff_SRR8797509_part_002.sam > shuff_SRR8797509_part_002.bam
samtools view -bS shuff_SRR8797509_part_003.sam > shuff_SRR8797509_part_003.bam
samtools view -bS shuff_SRR8797509_part_004.sam > shuff_SRR8797509_part_004.bam
samtools view -bS shuff_SRR8797509_part_005.sam > shuff_SRR8797509_part_005.bam

## convert the BAM file to a sorted BAM file. 
samtools sort shuff_SRR8797509_part_001.bam -o shuff_SRR8797509_part_001.sorted.bam
samtools sort shuff_SRR8797509_part_002.bam -o shuff_SRR8797509_part_002.sorted.bam
samtools sort shuff_SRR8797509_part_003.bam -o shuff_SRR8797509_part_003.sorted.bam
samtools sort shuff_SRR8797509_part_004.bam -o shuff_SRR8797509_part_004.sorted.bam
samtools sort shuff_SRR8797509_part_005.bam -o shuff_SRR8797509_part_005.sorted.bam

## Export some useful statistics report for each sample indvidually
for f in 1 2 3 4 5;
    do
        samtools flagstat shuff_SRR8797509_part_00$f.sorted.bam > shuff_useful_stat_$f.txt;
done


## install required tool for assembly 
## install stringtie

## Assembly without known annotations
for SAMPLE in 1 2 3 4 5;
    do

        stringtie shuff_SRR8797509_part_00${SAMPLE}.sorted.bam --rf -l ref_free_${SAMPLE} -o ref_free_${SAMPLE}.gtf
done

# Assembly with known previous annotations
for SAMPLE in 1 2 3 4 5;
    do
        stringtie shuff_SRR8797509_part_00${SAMPLE}.sorted.bam --rf -l ref_sup_${SAMPLE} -G ~/workdir/sample_data/chr22_with_ERCC92.gtf -o ref_sup_${SAMPLE}.gtf 
done

















# Using GTF-Compare to Compare the Generated Annotation Files to a Reference Annotation.

#Create virtual evironment with conda


conda create -n ngs-gtf python=3.6 anaconda
source activate ngs-gtf
conda install -c conda-forge pypy3.5
wget https://bootstrap.pypa.io/get-pip.py
pypy3 get-pip.py


pypy3 -m pip install gffutils numpy tqdm 'intervaltree<3.0'


mkdir -p ~/workdir/Assignment/gtf-compare/gtfs && cd ~/workdir/Assignment/gtf-compare/gtfs
ln -s ~/workdir/Assignment/bwa_align/ref_sup_*.gtf .
ln -s ~/workdir/sample_data/chr22_with_ERCC92.gtf .

mkdir -p ~/workdir/Assignment/gtf-compare/method_one && cd ~/workdir/Assignment/gtf-compare/method_one
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/comp.py
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/stat.py

for f in 1 2 3 4 5;
      do
              pypy3 comp.py -r ../gtfs/ref_sup_*.gtf ../gtfs/chr22_with_ERCC92.gtf
              pypy3 stat.py



## Differential_expression
mkdir -p ~/workdir/assignment/diff_exp && cd ~/workdir/assignment/diff_exp/
mkdir ~/workdir/assignment/ngs1_project/out && cd ~/workdir/assignment/ngs1_project/out|mv ~/workdir/assignment/ngs1_project/main_reads/out/.fastq out| mv ~/workdir/assignment/ngs1_project/shuff_reads/out/.fastq out
#wget -c https://0x0.st/zK57.gz -O ref.tar.gz
#tar xvzf ref.tar.gz
#wget -c https://raw.githubusercontent.com/mr-eyes/nu-ngs01/master/Day-6/deseq1.r
#wget -c https://raw.githubusercontent.com/mr-eyes/nu-ngs01/master/Day-6/draw-heatmap.r#1 Setup enviornemnt
#conda activate ngs1
# conda install kallisto
# conda install samtools# Install subread, we will use featureCount : a software program developed for counting reads to genomic features such as genes, exons, promoters and genomic bins.
conda install subread# install r and dependicies
#conda install r
conda install -y bioconductor-deseq r-gplots
RUNLOG=runlog.txt#Step 2 (Quantification)Step
GTF=~/workdir/sample_data/chr22_with_ERCC92.gtf# Generate the counts.
featureCounts -a "$GTF" -g gene_name -o counts.txt  ~/workdir/Assignment/bwa_align/main_SRR8797509*.bam  ~/workdir/Assignment/hisat_align/shuff_SRR8797509*.bam
# Simplify the file to keep only the count columns.
cat counts.txt | cut -f 1,7-12 > simple_counts.txt
# Analyze the counts with DESeq1.
cat simple_counts.txt | Rscript deseq1.r 5x5 > results_deseq1.tsv#View only rows with pval < 0.05
cat results_deseq1.tsv | awk ' $8 < 0.05 { print $0 }' > filtered_results_deseq1.tsv
cat filtered_results_deseq1.tsv | Rscript draw-heatmap.r > hisat_output.pdf


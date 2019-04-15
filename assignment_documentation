
#!/bin/bash
# created a new reprository on github "NGS1-Assignment"
# created a new file "my_info" in the NGS1-Assignment using this code:
git clone https://github.com/HowidaSamir/Assignment-NGS1
cd Assignment-NGS1
echo "Name: Howida Samir Mohamed Fathi Hussein" > my_user.md
echo "Student ID: 181002" >> my_user.md 
echo "NU E-mail: H.Samir@nu.edu.eg" >> my_user.md 
git add my_user.md 
git commit -m "file created"
git push

#1.data download
mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1 && cd ~/home/ngs-01/nu-ngs01/Assignment_NGS1
mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/raw_data && cd ~/home/ngs-01/nu-ngs01/Assignment_NGS1/raw_data
wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra

#2.converting the SRA file to FASTQ file:
#SRA tools installation:
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.4.1/sratoolkit.2.4.1-ubuntu64.tar.gz
tar xzvf sratoolkit.2.4.1-ubuntu64.tar.gz


#conversion of 5M samples from SRA to FASTQ:
mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/raw_data/unshuffled5Msamples && cd ~/home/ngs-01/nu-ngs01/Assignment_NGS1/raw_data/unshuffled5Msamples

#extraction of the 5M fragments and splitting them into R1 & R2:
fastq-dump -I -X 5000000 --split-files SRR8797509 

#splitting the 5M extracted fragments into 5 samples(both R1 & R2) (unshuffled samples):
$ seqkit split2 -1 SRR8797509_1.fastq -2 SRR8797509_2.fastq -p 5 -O out -f

#shuffling R1 & R2 then splitting them into 5 shuffled samples:
mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/ShuffledSamples && cd ~mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/ShuffledSamples

seqkit shuffle -2 SRR8797509_1.fastq 

seqkit shuffle -2 SRR8797509_2.fastq

seqkit split2 -1 ShuffledR1.fastq -2 shuffledR2.fastq -p 5 -O ShuffledSamples -f

#3.FASTQC:
source activate ngs1
conda install -c bioconda fastqc 
conda install -c bioconda multiqc 

#copy all sample1 (R1 & R2, shuffled & unshuffled to one folder):
mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/sample_1_all && cd ~/home/ngs-01/nu-ngs01/Assignment_NGS1/sample_1_all

cp -R ShuffledSamples/shuffled_SRR8797509_1.part_001.part_001.fastq.gz  ShuffledSamples/shuffled_SRR8797509_2.part_001.part_001.fastq.gz  sample1all/

cp -R unshuffled5Msamples/SRR8797509_1.part_001.fastq  unshuffled5Msamples/SRR8797509_2.part_001.fastq  sample1all/

cd sample1all
for f in *fastq;
do
fastqc -t 1 -f fastq $f;
done

# merge all the fastq analysis results in one report:
multiqc -z -o . .

#4. The trimming:
source activate ngs1
conda install -c bioconda trimmomatic 

#4.1. Mild trimming for the unshuffled samples:
 
mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/unshuffled_mildtrimming && cd ~/home/ngs-01/nu-ngs01/Assignment_NGS1/unshuffled_mildtrimming

cd unshuffled5Msamples

for ((i=1;i<=5;i++)) ;
do
cp -R unshuffled5Msamples/SRR8797509_1.part_00$i.fastq  unshuffled5Msamples/SRR8797509_2.part_00$i.fastq  unshuffled_mildtrimming/
done

cd unshuffled_mildtrimming

mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/unshuffled_fastqc/SX_1_miled_trimmed
mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/unshuffled_fastqc/SX_1_miledtrimmed_unwanted

for ((i=1;i<=5;i++)) ; do
SX_1unsh="/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/unshuffled5Msamples/SRR8797509_1.part_00$i.fastq"
SX_2unsh="/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/unshuffled5Msamples/SRR8797509_2.part_00$i.fastq"
newSX_1t="/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/unshuffled5Msamples/SX_1_miled_trimmed/SRR8797509_1.part_00$i.fastq"
newSX_2t="/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/unshuffled5Msamples/SX_1_miled_trimmed/SRR8797509_2.part_00$i.fastq"
newSX_1up="/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/unshuffled5Msamples/SX_1_miledtrimmed_unwanted/SRR8797509_1.part_00$i.fastq"
newSX_2up="/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/unshuffled5Msamples/SX_1_miledtrimmed_unwanted/SRR8797509_.part_00$i.fastq"
adap="/home/ngs-01/miniconda3/envs/ngs1/share/trimmomatic-0.39-1/adapters"
trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile "$SX_1unsh" "$SX_2unsh" "$newSX_1t" "$newSX_1up" "$newSX_2t" "$newSX_2up" ILLUMINACLIP:"$adap"/TruSeq3-PE.fa:2:30:10:1 MINLEN:36
done

#4.2. Aggressive trimming for shuffled samples:

mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/shuffled_agrtrimming && cd ~/home/ngs-01/nu-ngs01/Assignment_NGS1/shuffled_agrtri$

cd ShuffledSamples

for ((i=1;i<=5;i++)) ;
do
cp -R ShuffledSamples/shuffled_SRR8797509_1.part_001.part_00$i.fastq  shuffled_SRR8797509_2.part_001.part_00$i.fastq  shuffled_agrtrimming/
done

cd shuffled_agrtrimming

mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/unshuffled_fastqc/SX_1_agr_trimmed
mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/unshuffled_fastqc/SX_1_agrtrimmed_unwanted

for ((i=1;i<=5;i++)) ;
do
SX_1sh="/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/ShuffledSamples/shuffled_SRR8797509_1.part_001.part_00$i.fastq"
SX_2sh="/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/ShuffledSamples/shuffled_SRR8797509_2.part_001.part_00$i.fastq"
newusSX_1t="/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/ShuffledSamples/SX_1_agr_trimmed/shuffled_SRR8797509_1.part_001.part_00$i.fastq"
newusSX_2t="/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/ShuffledSamples/SX_1_agr_trimmed/shuffled_SRR8797509_2.part_001.part_00$i.fastq"
newusSX_1up="/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/ShuffledSamples/SX_1_agrtrimmed_unwanted/shuffled_SRR8797509_1.part_001.part_00$i.fastq"
newusSX_2up="/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/ShuffledSamples/SX_1_agrtrimmed_unwanted/shuffled_SRR8797509_2.part_001.part_00$i.fastq"
adap="/home/ngs-01/miniconda3/envs/ngs1/share/trimmomatic-0.39-1/adapters"
trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile $SX_1sh "$SX_2sh" "$newusSX_1t" "$newusSX_1up" "$newusSX_2t" "$newusSX_2up" ILLUMINACLIP:"$adap"/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36
done

#5. the alignment:
#5.1.BWA alignment for the unshuffled samples:

#download the reference files:
 
mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/chr22data && cd ~/home/ngs-01/nu-ngs01/Assignment_NGS1/mkdir -p ~/home/ngs-01/nu-ngs01/Assignment_NGS1/shuffled_agrtrimming && cd ~/home/ngs-01/nu-ngs01/Assignment_NGS1/chr22data

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_transcripts.fa.gz

gunzip gencode.v29.pc_transcripts.fa.gz

wget https://transfer.sh/IbpI7/HBR_UHR_ERCC_ds_5pc.tar

tar -xvf HBR_UHR_ERCC_ds_5pc.tar


#OR download chr22 data:
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa

gunzip chr22_with_ERCC92.fa.gz

wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf

wget https://transfer.sh/IbpI7/HBR_UHR_ERCC_ds_5pc.tar

tar -xvf HBR_UHR_ERCC_ds_5pc.tar

#OR: 

cd ~/workdir/sample_data
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_transcripts.fa.gz
gunzip gencode.v29.pc_transcripts.fa.gz

wget https://transfer.sh/IbpI7/HBR_UHR_ERCC_ds_5pc.tar
tar -xvf HBR_UHR_ERCC_ds_5pc.tar

# Download the Transcriptome Annotation File
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
gunzip gencode.v29.annotation.gtf.gz

# Select the transcripts of Chr22
cd ~/workdir/sample_data/
source activate ngs1
READS=$(grep "^chr22" gencode.v29.annotation.gtf | awk -F'\t' '{print $9}' | awk -F';' '{print $1}' | awk -F' ' '{print $2}' | awk -F'"' '{print $2}' | sort | uniq)

for value in $READS
    do  
        echo "Processing: $value"
        seqkit grep -r -p ${value} gencode.v29.pc_transcripts.fa | awk -F'|' '{print $1}' >> gencode.v29.pc_transcripts.chr22.simplified.fa
    done
#honestly, i had the data already downloaded so i wasn't in need to do this step anyway :)

#install BWA:
source activate ngs1
conda install -c bioconda bwa

#indexing:
mkdir -p ~/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/chr22data/bwaIndex && cd ~/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/chr22data/bwaIndex

ln -s /home/ngs-01/nu-ngs01/Assignment-NGS1/RawData/unshuffled5Msamples/SX_1_miled_trimmed/gencode.v29.pc_transcripts.chr22.simplified.fa .

bwa index -a bwtsw /home/ngs-01/nu-ngs01/Assignment-NGS1/RawData/chr22data/gencode.v29.pc_transcripts.chr22.simplified.fa .

#the actual BWA alignment step:
mkdir -p ~/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/bwa_align && cd ~/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/bwa_align

for ((i=1;i<=5;i++)) ;

do

R1="$home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/unshuffled5Msamples/SX_1_miled_trimmed/SRR8797509_1.part_00$i.fastq"

R2="$/home/ngs-01/nu-ngs01/Assignment-NGS1/Raw Data/unshuffled5Msamples/SX_1_miled_trimmed/SRR8797509_2.part_00$i.fastq"

/usr/bin/time -v bwa mem bwaIndex/gencode.v29.pc_transcripts.chr22.simplified.fa $R1 $R2 > SX_1_SRR8797509_part_00$i.sam

done


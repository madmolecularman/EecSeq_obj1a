# This code repeats the analyses done in EAGER Obj1a 2022

This code assumes that there are two subdirectories in this working directory ./DNA which contain the raw sequencing files from the EAGER 2022 Obj1a sequenicng project and  ./Genome directory which contains files from the eastern oyster genome.

## Directory Structure

01_process_reads: start of the project and where all fastq files are placed
* barcodes: contains barcode files for demultiplexing
* clean:  clean fastq files
* demux:  demultiplexed fastq files
* raw:  raw fastqfiles

02_ddocent: a folder to store the output od a dDocent run. Additionally contains trimmed forward and reverse fq.gz, F.bam, and RGmd.bam files.

03_mapping: a folder to store file outputs from mapping reads to the oyster genome. Files include merged capture F.bam files, linked F.bam files, genome.depth profiles, and a lists of cds, gene, and exon coordinates from the oyster gff3 annotations. 

04_coverage_analysis
* 01_genome_region
* 02_exon_stats
* 03_target_interval
* 04_specificity

05_snp_calling: contains output from vcftools analysis

Genome:contains all reference files needed for the project

## Set working directory

```bash
WORKING_DIR = /home/jgreen/EAGERobj1a/
echo $WORKING_DIR
cd $WORKING_DIR
```

## Setup conda environment
```bash
conda update conda
mamba create eagerobj1a ddocent
conda activate eagerobj1a
```

## Copy full haplotype masked genome over to my Genome/ directory

```bash
#Copy genome
cp /RAID_STORAGE2/Shared_Data/Oyster_Genome/masked/masked.cvir.genome.fasta .
#Index genome
samtools faidx masked.cvir.genome.fasta
#Make genome.file
mawk -v OFS='\t' {'print $1,$2'} masked.cvir.genome.fasta > masked.genome.file
```

## Copy over top.bed files

```bash
cp /home/Genome_Resources/C_virginica/bed_and_GFF/ref_C_virginica-3.0_top.bed
```
## Copying files to 01_process_reads/raw directory and renaming to fit with ddocent naming convention
```bash
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/Capture1_B3/*L1_1.fq.gz Capture1_B3.R1.fastq.gz
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/Capture1_B3/*L1_2.fq.gz Capture1_B3.R2.fastq.gz
```

## Enmasse copying files and renaming
```bash
# Set string array with the prefix name of each read file.
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N1" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N1" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N1" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N1" "Capture4_N2")

## Check naming convention
for i in "${StringArray[@]}"
do
echo $i
done

## Copy Forward 1 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L1_1.fq.gz ${i}.R1.fastq.gz
done

## Copy Forward 2 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L2_1.fq.gz ${i}.R1.fastq.gz
done

## Copy Forward 3 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L3_1.fq.gz ${i}.reseq.R1.fastq.gz
done

## Copy Forward 4 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L4_1.fq.gz ${i}.reseq.R1.fastq.gz
done

## Copy Reverse 1 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L1_2.fq.gz ${i}.R2.fastq.gz
done

## Copy Reverse 2 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L2_2.fq.gz ${i}.R2.fastq.gz
done

## Copy Reverse 3 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L3_2.fq.gz ${i}.reseq.R2.fastq.gz
done

## Copy Reverse 4 reads
for i in "${StringArray[@]}"
do
cp /RAID_STORAGE2/Raw_Data/EAGER_2022/Sequencing_Run_Nov_2022/01.RawData/${i}/${i}_*_L4_2.fq.gz ${i}.reseq.R2.fastq.gz
done
```

## Run fastqc and multiqc to visualize files
```bash
fastqc -f fastq -t 20 *.fastq.gz
multiqc -o raw_multiqc_run .
```

## Results
Multiqc plots show uneven distribution of sequence counts specifically in samples Capture1_N1.R1, Capture1_N1.R2, Capture2_N1.R1, Capture2_N1.R2, Capture3_N1.reseq.R1, Capture3_N1.reseq.R2, Capture4_N1.reseq.R1, Capture4_N1.reseq.R2. See [multiqc report](https://github.com/madmolecularman/EecSeq_obj1a/blob/main/Bioinformatics/01_process_reads/multiqc_report_20230301_run1.html)

## Renaming reads

I need to be careful with some changes to the reads. 
1) Capture 2 = Capture 1
2) Capture 1 = Capture 2
3) Capture 3 = Capture 3
4) Capture 4 = Capture 4

Move capture 1 and capture 2 files into different directories
```bash
mkdir capture1 capture2

mv Capture2* capture1/
mv Capture1* capture2/
```

Rename Capture1 preffix to Capture2 preffix
```bash
rename Capture1 Capture2 Capture1*
rename Capture2 Capture1 Capture2*
```

Redo fastqc and multiqc
```bash
rm *fastqc*
fastqc -f fastq -t 20 *.fastq.gz
multiqc -o rename_multiqc .
```

## Fixing mixed files

Capture1_G3.R1, Capture1_G3.R2, Capture1_G3.reseq.R1, Capture1_G3.reseq.R2 are actually capture2

Capture2_G3.R1, Capture2_G3.R2, Capture2_G3.reseq.R1, Capture2_G3.reseq.R2 are actually capture1

Capture2_N2.R1, Capture2_N2.R2, Capture2_N2.reseq.R1, Capture2_N2.reseq.R1 are actually capture1

Capture1_N2.R1, Capture2_N2.R2, Capture2_N2.reseq.R1, Capture2_N2.reseq.R1 are actually capture2

```bash
mv Capture2_G3* capture1/
rename Capture2 Capture1 Capture2*
mv Capture1_G3* capture2/
rename Capture1 Capture2 Capture1*
mv Capture2_N2 capture1/
rename Capture2 Capture1 Capture2*
mv Capture1_N2 capture2/
rename Capture1 Capture2 Capture1*
```

Redo fastqc and multiqc
```bash
rm *fastqc*
fastqc -f fastq -t 20 *.fastq.gz
multiqc -o unmix_multiqc .
```

Now we have all the files organized which you can see in this [multiqc]() report. Most of the files are demultiplexed. Yet Capture[1-4]_N1 and paired reseq files have very messy 

## Processing short reads on files that are not demultiplexed

## Demultiplexing capture files

Use process_shortreads to demultiplex files from raw/ to the demux/samples/ folder. One of the major problems is that we are not sure if novogene was able to separate out the different capture pools specifically for Capture[1-4]_N1 samples and Capture[1-4]_N1 resequenced. We need to create a demultiplexing run that will filter the index and the inline barcode into the proper file. The question here is can I demultiplex just the index and then the adapter inline barcode?

Its important to understand what ![process_shortreads](https://catchenlab.life.illinois.edu/stacks/comp/process_shortreads.php) is doing and how to create the barcode files ![Stacks manual](https://catchenlab.life.illinois.edu/stacks/manual/index.php#clean)


## Using duplicated barcodes and inline_inline to invetigate ambigous barcode drop for R2 files

New directory for double inline barcodes. I am wanting to test out how we retain reads with double barcodes.

Make new barcode file

``bash
nano double_barcode.txt
```

```text
ATCGCG  ATCGCG	Capture1
CGATGT  CGATGT	Capture2
TTAGGC  TTAGGC	Capture3
TGGCCA  TGGCCA	Capture4
```

```bash
mkdir double_barcode_samples
cd double_barcode_samples
```

Making directories for Capture[1-4]_N1 and reseq files
```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N1" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N1" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N1" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N1" "Capture4_N2")

for i in "${StringArray[@]}"
do
mkdir ${i}
done
```

Also need to make reseq directories for the Capture[1-4]_N1

```bash
mkdir Capture1_N1_reseq Capture2_N1_reseq Capture3_N1_reseq Capture4_N1_reseq
```

## Decontaminating our Capture[1-4]_N1 reads

Use double barcode file located here: 

> /home/jgreen/EAGER_OBJ1a/01_process_reads/barcodes/double_barcode.txt

Code for process_shortreads on the Capture[1-4]_N1 and reseq files

```bash
declare -a StringArray=("Capture1_N1" "Capture2_N1" "Capture3_N1" "Capture4_N1")
for i in "${StringArray[@]}"
do
process_shortreads -P -1 $WORKING_DIR/01_process_reads/raw/${i}.R1.fastq.gz -2 $WORKING_DIR/01_process_reads/raw/${i}.R2.fastq.gz -o $WORKING_DIR/01_process_reads/demux/${i}/ -b $WORKING_DIR/01_process_reads/barcodes/double_barcode.txt --inline_inline -r
done

declare -a StringArray=("Capture1_N1.reseq" "Capture2_N1.reseq" "Capture3_N1.reseq" "Capture4_N1.reseq")
for i in "${StringArray[@]}"
do
process_shortreads -P -1 $WORKING_DIR/01_process_reads/raw/${i}.R1.fastq.gz -2 $WORKING_DIR/01_process_reads/raw/${i}.R2.fastq.gz -o $WORKING_DIR/01_process_reads/demux/${i}/ -b $WORKING_DIR/01_process_reads/barcodes/double_barcode.txt --inline_inline -r
done
```

## Fastqc and multiqc demux files

```bash
declare -a StringArray=("Capture1_N1" "Capture1_N1.reseq" "Capture2_N1" "Capture2_N1.reseq" "Capture3_N1" "Capture3_N1.reseq" "Capture4_N1" "Capture4_N1.reseq")

for i in "${StringArray[@]}"
do
fastqc -f fastq -t 20 ${i}/*.fq.gz
done

for i in "${StringArray[@]}"
do
multiqc ${i}/. -o ${i}_multiqc
done
```

Make process_shortreads file that is readable

```bash
for i in "${StringArray[@]}"
do
head -n 100 ${i}/process_shortreads.log > ${i}/${i}process_shortreads.short.log
done
```

## Cat reseq reads to original reads

```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")

for i in "${StringArray[@]}"
do
cp $WORKING_DIR/01_process_reads/raw/${i}.R1.fastq.gz $WORKING_DIR/01_process_reads/demux/${i}/
cp $WORKING_DIR/01_process_reads/raw/${i}.R2.fastq.gz $WORKING_DIR/01_process_reads/demux/${i}/
cat $WORKING_DIR/01_process_reads/raw/${i}.reseq.R1.fastq.gz >> $WORKING_DIR/01_process_reads/demux/${i}/${i}.R1.fastq.gz
cat $WORKING_DIR/01_process_reads/raw/${i}.reseq.R2.fastq.gz >> $WORKING_DIR/01_process_reads/demux/${i}/${i}.R2.fastq.gz
done
```

## Demultiplex all the other reads

```bash
for i in "${StringArray[@]}"
do
process_shortreads -P -1 $WORKING_DIR/01_process_reads/demux/${i}/${i}.R1.fastq.gz -2 $WORKING_DIR/01_process_reads/demux/${i}/${i}.R2.fastq.gz -o $WORKING_DIR/01_process_reads/demux/${i}/ -b $WORKING_DIR/01_process_reads/barcodes/double_barcode.txt --inline_inline -r
done
```

## Make all of the process_shortreads short logs for review

```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")
for i in "${StringArray[@]}"
do
head -n 100 ${i}/process_shortreads.log > ${i}/${i}process_shortreads.short.log
done
```

Add Capture[1-4]_N1 and reseq reads together (rework this 3/10)

```bash
```bash
#Capture1_N1
cat Capture2.1.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture2_N1.F.fq.gz
cat Capture2.2.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture2_N1.R.fq.gz
#Capture1_N1_reseq 
cat Capture2.1.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture2_N1.F.fq.gz
cat Capture2.2.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture2_N1.R.fq.gz
#Capture2_N1
cat Capture2.1.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture2_N1.F.fq.gz
cat Capture2.2.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture2_N1.R.fq.gz
#Capture2_N1_reseq
cat Capture2.1.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture2_N1.F.fq.gz
cat Capture2.2.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture2_N1.R.fq.gz
#Capture3_N1
cat Capture3.1.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture3_N1.F.fq.gz
cat Capture3.2.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture3_N1.R.fq.gz
#Capture3_N1_reseq
cat Capture3.1.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture3_N1.F.fq.gz
cat Capture3.2.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture3_N1.R.fq.gz
#Capture4_N1
cat Capture4.1.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture4_N1.F.fq.gz
cat Capture4.2.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture4_N1.R.fq.gz
#Capture4_N1_reseq
cat Capture4.1.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture4_N1.F.fq.gz
cat Capture4.2.fq.gz >> $WORKING_DIR/01_process_reads/clean/Capture4_N1.R.fq.gz

```

Need to pull Capture[1-4] reseq reads into the clean directory

```bash
#Capture1
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demux/${i}/Capture1.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demux/${i}/Capture1.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

#Capture2
declare -a StringArray=( "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demux/${i}/Capture2.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demux/${i}/Capture2.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

#Capture3
declare -a StringArray=("Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demux/${i}/Capture3.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demux/${i}/Capture3.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

#Capture4
declare -a StringArray=("Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demux/${i}/Capture4.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demux/${i}/Capture4.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done
```

Using chatGPT we can simpligy this code to the following

```bash
# List of prefixes for each capture
declare -a captures=("Capture1" "Capture2")

for capture in "${captures[@]}"
do
    # List of samples for each capture
    declare -a samples=("N1" "N2")
    
    for sample in "${samples[@]}"
    do
        # Construct the sample ID
        sample_id="${capture}_${sample}"
        
        # Copy the files
        cat "$WORKING_DIR/01_process_reads/demux/${sample_id}/${capture}.1.fq.gz" >> "$WORKING_DIR/01_process_reads/clean/Capture2_N2.F.fq.gz"
        cat "$WORKING_DIR/01_process_reads/demux/${sample_id}/${capture}.2.fq.gz" >> "$WORKING_DIR/01_process_reads/clean/Capture2_N2.R.fq.gz"
        cat "$WORKING_DIR/01_process_reads/demux/${sample_id}.reseq/${capture}.1.fq.gz" >> "$WORKING_DIR/01_process_reads/clean/Capture2_N2.F.fq.gz"
        cat "$WORKING_DIR/01_process_reads/demux/${sample_id}.reseq/${capture}.2.fq.gz" >> "$WORKING_DIR/01_process_reads/clean/Capture2_N2.R.fq.gz"
    done
done
```

Just need to redo capture2_N1 and capture2_n2, capture2_n1 takes reads from capture1_n1, capture1_n1.reseq, capture2_n1, capture2_n1.reseq

## Analyze clean files

## Run fastqc and multiqc to visualize files

```bash
fastqc -f fastq -t 20 *.fq.gz
multiqc . -o post_demux_multiqc_report
```

## Create a reference.file from reference.fasta that is originally the masked genome fasta

```bash
mv masked.genome.fasta $WORKING_DIR/02_ddocent/reference.fasta
samtools faidx reference.fasta
mawk -v OFS='\t' {'print $1,$2'} reference.fasta.fai > masked.genome.file
```

## Run ddocent on all individual files

Make config.file

```bash
nano config.file

Number of Processors
24
Maximum Memory
0
Trimming
yes
Assembly?
no
Type_of_Assembly
PE
Clustering_Similarity%
0.90
Minimum within individual coverage level to include a read for assembly (K1)
2
Minimum number of individuals a read must be present in to include for assembly (K2)
2
Mapping_Reads?
yes
Mapping_Match_Value
1
Mapping_MisMatch_Value
3
Mapping_GapOpen_Penalty
5
Calling_SNPs?
yes
Email
gree9242@uri.edu
```

## Download new dDocent script in scripts directory
```bash
wget "https://github.com/The-Eastern-Oyster-Genome-Project/2022_Eastern_Oyster_Haplotig_Masked_Genome/blob/main/scripts/dDocent_ngs"
```

I receive an error when downloading this as most of the file is in an html format. Create a file using nano. Copy and paste the whole script from the github into this new file

```bash 
nano dDocent_ngs.sh
```

change permissions on the files to allow it to be executable

```bash
chmod +x dDocent_ngs.sh
```

## Run the dDocent script in the directory containing our newly named files

```bash
/../scripts/dDocent_ngs.sh config.file
```

## Masking old bed files with haplotig bed file

Copy old bed files

```bash
cp /home/Genomic_Resources/C_virginica/bed_and_GFF/*bed .
```

Pull haplotig bed from github repo https://github.com/The-Eastern-Oyster-Genome-Project/2022_Eastern_Oyster_Haplotig_Masked_Genome/blob/main/Haplotig_Masking/Haplotig_Masking.md

```bash
wget https://raw.githubusercontent.com/The-Eastern-Oyster-Genome-Project/2022_Eastern_Oyster_Haplotig_Masked_Genome/main/Haplotig_Masking/Output/haplotigs.bed
```

Use bedtools subtract to remove haplotigs present inthe bed file from the other bed files

```bash
declare -a StringArray=("sorted.ref3.0.CDS" "sorted.ref3.0.exon" "sorted.ref3.0.gene" "sorted.ref3.0.UTR" )
for i in "${StringArray[@]}"
do
bedtools subtract -a ${i}.bed -b haplotigs.bed > ${i}.hmask.bed
done

declare -a StringArray=("sorted.ref3.0.CDS.sc" "sorted.ref3.0.exon.sc" "sorted.ref3.0.gene.sc" "sorted.ref3.0.UTR.sc" )
for i in "${StringArray[@]}"
do
bedtools subtract -a ${i}.bed -b haplotigs.bed > ${i}.hmask.bed
done
```

## Link all the .F.bam files from ddocent analysis to the 03_mapping folder

```bash
ln -s $WORKING_DIR/02_ddocent/*.F.bam .
```


# Merge capture pools into one file

```bash
#Capture1
samtools merge -@20 capture1_filter.merged.F.bam Capture1_B3.F.bam Capture1_B4.F.bam Capture1_G3.F.bam Capture1_G5.F.bam Capture1_K3.F.bam Capture1_K4.F.bam Capture1_M3.F.bam Capture1_M4.F.bam Capture1_N2.F.bam

#Capture2
samtools merge -@20 capture2_filter.merged.F.bam Capture2_B3.F.bam Capture2_B4.F.bam Capture2_G3.F.bam Capture2_G5.F.bam Capture2_K3.F.bam Capture2_K4.F.bam Capture2_M3.F.bam Capture2_M4.F.bam Capture2_N2.F.bam

#Capture3
samtools merge -@20 capture3_filter.merged.F.bam Capture3_B3.F.bam Capture3_B4.F.bam Capture3_G3.F.bam Capture3_G5.F.bam Capture3_K3.F.bam Capture3_K4.F.bam Capture3_M3.F.bam Capture3_M4.F.bam Capture3_N2.F.bam

#Capture4
samtools merge -@20 capture4_filter.merged.F.bam Capture4_B3.F.bam Capture4_B4.F.bam Capture4_G3.F.bam Capture4_G5.F.bam Capture4_K3.F.bam Capture4_K4.F.bam Capture4_M3.F.bam Capture4_M4.F.bam Capture4_N2.F.bam
#All capture
samtools merge -@20 all.filter.merged.bam capture1_filter.merged.F.bam capture2_filter.merged.F.bam capture3_filter.merged.F.bam capture4_filter.merged.F.bam
```

## Get total coverage counts per exon for each capture

sc = single copy gene bed file
hmask = haplotype masked bed file

We will be using the new haplotype masked bed file that Jon has made after the analysis finished up on the C. virginica genome

Capture 1
```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")

declare -a StringArray=("capture1_filter.merged" "capture2_filter.merged" "capture3_filter.merged" "capture4_filter.merged")

#CDS
#Test one one capture
bedtools coverage -hist -b $WORKING_DIR/03_mapping/Capture1_B3.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.CDS.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/Capture1_B3.hist.AllCDS.all.split.txt
#Looping for all samples CDS
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/${i}.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.CDS.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/01_genome_region/${i}.hist.AllCDS.all.split.txt
done
#Exon
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/${i}.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/01_genome_region/${i}.hist.AllExon.all.split.txt
done
#Gene
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/${i}.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.gene.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/01_genome_region/${i}.hist.AllGene.all.split.txt
done
#UTR
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/${i}.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.UTR.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -sorted -split | grep ^all > $WORKING_DIR/04_coverage_analysis/01_genome_region/${i}.hist.AllUTR.all.split.txt
done
#Intergenic
for i in "${StringArray[@]}"
do
bedtools coverage -hist -b $WORKING_DIR/03_mapping/${i}.F.bam -a $WORKING_DIR/Genome/sorted.ref3.intergenic.bed -g $WORKING_DIR/Genome/masked.genome.file -sorted  -split | grep ^all > $WORKING_DIR/04_coverage_analysis/01_genome_region/${i}.hist.AllIntergenic.all.split.txt
done
```

R plot for one file

```R
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)

setwd("/home/jgreen/EAGER_OBJ1a/04_coverage_analysis/01_genome_region/")
# initialize an empty list to store ggplot objects
plot_list <- list()

files <- c("Capture1_B3.hist", "Capture2_B3.hist", "Capture3_B3.hist", "Capture4_B3.hist")

for (file in files) {
  print(file)
  files_list <- list.files(pattern=file)
  print(files_list)
  
  files <- c(paste0(file, ".AllCDS.all.split.txt"), paste0(file, ".AllExon.all.split.txt"), paste0(file, ".AllGene.all.split.txt"), paste0(file, ".AllUTR.all.split.txt"))
  print(files)

  labs <- c("CDS","Exon","Gene","UTR")
  
  cov <- list()
  for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])[,c(2,5)]
    cov_cumul=1-cumsum(cov[[i]][,2])
    cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
    cov[[i]]$sample=labs[i]
  }

  cov_df=do.call("rbind",cov)
  names(cov_df)[1:2]=c("depth","fraction")

  pcbPalette <- c("#009E73" ,"#D55E00","#CC79A7","#56B4E9")

  p1 <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  + xlim(0,25)+
    scale_alpha(guide = 'none') +
    geom_line(size=1.5)+ 
    #geom_segment(aes(x=20, y=0, xend=20, yend=1, color="red"))+
    scale_color_manual(values=pcbPalette) +
    scale_fill_manual(values=pcbPalette) +
    ggtitle(file)+
    ylab("% of Bases > Depth")+
    xlab("Depth")+
    theme_bw() +
    theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
    theme(legend.title = element_blank()) + 
    theme(legend.position=c(0.50,0.75))
  
  plot_list[[length(plot_list)+1]] <- p1
  
  png(filename=paste0("Figure2_", file, ".png"), type="cairo",units="px", width=5600, 
      height=3000, res=600, bg="transparent")
  
  print(p1)
  
  dev.off()

}

# combine all the ggplot objects in the list using grid.arrange()
pcombined <- grid.arrange(grobs = plot_list, ncol = 2, top=textGrob("Capture B3 Read Depth across genome"))
#save the combined plot as PNG file
ggsave(file="Figure2_allcapture_B3.png", pcombined)
```

Lets take a look at one figure 

![Capture_N2](/home/jgreen/EAGER_OBJ1a/04_coverage_analysis/01_genome_region/Figure2_allcapture_N2.png)

These values are a lot lower than when I first did the analysis on just the single copy bed files. Additionally the 0 depth the % bases > 1 is 1. I don't know if this makes sense/ I am also missing a bam file for Capture4_n1.

# Generate data for Table with mapped read stats

```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")

for i in "${StringArray[@]}"
do 
nom=$(samtools view -@32 $WORKING_DIR/03_mapping/${i}.F.bam -c -L $WORKING_DIR/Genome/mtDNA.bed); denom=$(samtools view -@32 $WORKING_DIR/03_mapping/${i}.F.bam -c); dup=$(mawk '/Unknown/' $WORKING_DIR/02_ddocent/logfiles/${i}_dup_metrics.txt | cut -f9); paste <(echo $i) <(echo $(( `zcat $WORKING_DIR/01_process_reads/raw/${i}.R1.fastq.gz | wc -l` /4 ))) <(echo $(( `zcat $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz | wc -l` /4 ))) <(samtools view -@ 32 $WORKING_DIR/02_ddocent/${i}-RGmd.bam -c) <(python -c "print(round("$dup" * 100,2))") <(echo $denom) <(python -c "print(round("$nom"/"$denom" *100,2))") 
done > data.table2
echo -e "Pool\tRaw_Reads\tFiltered_Reads\tMapped_Reads\t%_Duplicate\tFiltered_Mapped_Reads\t%_mapping_to_mitochondrial_genome" > header
cat header data.table2 > table2.txt
```

## Generate data for Figure 3

Calculate Exon percentiles

Get total coverage counts for all merged

```bash
bedtools coverage -b $WORKING_DIR/03_mapping/all.filter.merged.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -mean -split > $WORKING_DIR/04_coverage_analysis/02_exon_stats/all.merged.cov.mean.exon.stats
```

Get total coverage counts per exon per capture
```bash
declare -a StringArray=("capture1_filter.merged" "capture2_filter.merged" "capture3_filter.merged" "capture4_filter.merged")

for i in "${StringArray[@]}"
do
bedtools coverage -b $WORKING_DIR/03_mapping/${i}.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed -g $WORKING_DIR/Genome/masked.genome.file -mean -split > $WORKING_DIR/04_coverage_analysis/02_exon_stats/${i}.cov.mean.filtered.exon.stats
done
```

Get total coverage counter per exon for each sample
```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")
for i in "${StringArray[@]}"
do
bedtools coverage -b $WORKING_DIR/03_mapping/${i}.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed -sorted -g $WORKING_DIR/Genome/masked.genome.file -mean > $WORKING_DIR/04_coverage_analysis/02_exon_stats/${i}.cov.mean.exon.stats
done
```

Get each capture name with .F.bam using string array
```bash
for i in "${StringArray[@]}"
do
echo $i.F.bam
done
```

Goal is to compare capture1 against all others

Remove mtDNA from stat files

```bash
# Individuals files
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")

for i in  "${StringArray[@]}"
do
mawk '!/NC_007175.2/' ${i}.cov.mean.exon.stats > ${i}.DNA.mean.exon.stats
done

#Merged files
declare -a StringArray=("capture1_filter.merged" "capture2_filter.merged" "capture3_filter.merged" "capture4_filter.merged")

for i in  "${StringArray[@]}"
do
mawk '!/NC_007175.2/' ${i}.cov.mean.filtered.exon.stats > ${i}.DNA.merged.mean.exon.stats
done

#All files
mawk '!/NC_007175.2/' all.merged.cov.mean.exon.stats > all.DNA.merged.mean.exon.stats
```

Calculate lower 10th percentile of exon sizes

```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")

for i in  "${StringArray[@]}"
do
mawk '{print $3 -$2}' ${i}.DNA.mean.exon.stats |sort -g | perl -e '$d=.1;@l=<>;print $l[int($d*@l)]'
done
```

Result: 58

```bash
declare -a StringArray=("capture1_filter.merged" "capture2_filter.merged" "capture3_filter.merged" "capture4_filter.merged")
for i in  "${StringArray[@]}"
do
mawk '{print $3 -$2}' ${i}.DNA.merged.mean.exon.stats |sort -g | perl -e '$d=.1;@l=<>;print $l[int($d*@l)]'
done
```

Result: 59

```bash
mawk '{print $3 -$2}' all.DNA.merged.mean.exon.stats |sort -g | perl -e '$d=.1;@l=<>;print $l[int($d*@l)]'
```

Result: 58


Calculate upper 10th percentile of exon sizes

```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")

for i in  "${StringArray[@]}"
do
mawk '{print $3 -$2}' ${i}.DNA.mean.exon.stats |sort -g | perl -e '$d=.9;@l=<>;print $l[int($d*@l)]'
done
```

Result: 523

```bash
declare -a StringArray=("capture1_filter.merged" "capture2_filter.merged" "capture3_filter.merged" "capture4_filter.merged")
for i in  "${StringArray[@]}"
do
mawk '{print $3 -$2}' ${i}.DNA.merged.mean.exon.stats |sort -g | perl -e '$d=.9;@l=<>;print $l[int($d*@l)]'
done
```

Result: 517

```bash
mawk '{print $3 -$2}' all.DNA.merged.mean.exon.stats |sort -g | perl -e '$d=.9;@l=<>;print $l[int($d*@l)]'
```

Result: 524

Mark exons into size classes based on size distribution and create data table
```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")

for i in  "${StringArray[@]}"
do
mawk '{if ( $3 -$2 > 524 ) print $0 "\tUpper"; else if ( $3 - $2 < 58 ) print $0 "\tLower"; else if ( $3 - $2 > 58 && $3 - $2 < 524) print $0 "\tMiddle" }' ${i}.DNA.mean.exon.stats > ${i}.mean.cov.exon.stats.class
done

echo -e "Chrom\tStart\tEnd\tDNA_Coverage\tExon_Size_Class" > header

for i in  "${StringArray[@]}"
do
cat header ${i}.mean.cov.exon.stats.class > ${i}.ExonMeanCoverage.txt
done
```

```bash
declare -a StringArray=("capture1_filter.merged" "capture2_filter.merged" "capture3_filter.merged" "capture4_filter.merged")

for i in  "${StringArray[@]}"
do
mawk '{if ( $3 -$2 > 524 ) print $0 "\tUpper"; else if ( $3 - $2 < 58 ) print $0 "\tLower"; else if ( $3 - $2 > 58 && $3 - $2 < 524) print $0 "\tMiddle" }' ${i}.DNA.merged.mean.exon.stats > ${i}.merged.mean.cov.exon.stats.class
done


for i in  "${StringArray[@]}"
do
cat header ${i}.merged.mean.cov.exon.stats.class > ${i}.merged.ExonMeanCoverage.txt
done
```

```bash
mawk '{if ( $3 -$2 > 524 ) print $0 "\tUpper"; else if ( $3 - $2 < 58 ) print $0 "\tLower"; else if ( $3 - $2 > 58 && $3 - $2 < 524) print $0 "\tMiddle" }' all.DNA.merged.mean.exon.stats > all.merged.mean.cov.exon.stats.class
cat header all.merged.mean.cov.exon.stats.class > all.merged.ExonMeanCoverage.txt
```

```r
library(MASS)
library(fields)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)
```

```r Make data frames
setwd("/home/jgreen/EAGER_OBJ1a/04_coverage_analysis/02_exon_stats")
df1 <- read.table("capture1_filter.merged.merged.ExonMeanCoverage.txt", header = TRUE)
df1 <-as.data.frame(df1)
df2 <- read.table("capture2_filter.merged.merged.ExonMeanCoverage.txt", header = TRUE)
df2 <-as.data.frame(df2)
df3 <- read.table("capture3_filter.merged.merged.ExonMeanCoverage.txt", header = TRUE)
df3 <-as.data.frame(df3)
df4 <- read.table("capture4_filter.merged.merged.ExonMeanCoverage.txt", header = TRUE)
df4 <-as.data.frame(df4)
df5 <- read.table("all.merged.ExonMeanCoverage.txt", header = TRUE)
df5 <-as.data.frame(df5)
```

# join the two data frames based on the id column
```r Merge12
setwd("/home/jgreen/EAGER_OBJ1a/04_coverage_analysis/02_exon_stats")
merged_cap12df <- merge(df1, df2[, c("Start","DNA_Coverage", "Exon_Size_Class")],by = "Start")
TotalExon <- merged_cap12df[merged_cap12df$DNA_Coverage.x != 0 & merged_cap12df$DNA_Coverage.y != 0,]
TotalExon <- TotalExon[, -5]

TotalExon$Exon_Size_Class <-factor(TotalExon$Exon_Size_Class.y, levels=c("Lower","Middle","Upper"))


TotalExon$Exon_Size_Class <- revalue(TotalExon$Exon_Size_Class.y, c("Lower"="Lower 10%", "Upper"="Upper 10%", "Middle"="Middle 80%"))

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
TotalExon$density <- get_density(TotalExon$DNA_Coverage.x,TotalExon$DNA_Coverage.y)

cbPalette <- c("#009E73","#D55E00","#56B4E9", "#0072B2","#E69F00", "#F0E442" , "#999999","#CC79A7")
t <- cbPalette[1]
cbPalette[1] <- cbPalette[2]
cbPalette[2] <- t

TotalExon$DNA_Coverage.y <- TotalExon$DNA_Coverage.y/6
TotalExon$DNA_Coverage.x <- TotalExon$DNA_Coverage.x/6


b1 <- bb<- ggplot(TotalExon, aes(x=DNA_Coverage.x+1,y=DNA_Coverage.y+1,alpha = 1/density),fill=Exon_Size_Class,color=Exon_Size_Class) + 
  geom_point(aes(color=TotalExon$Exon_Size_Class,fill=TotalExon$Exon_Size_Class,shape=TotalExon$Exon_Size_Class)) +
  geom_smooth(method="auto", alpha = 0.5, size = 0, se=TRUE)+
  stat_smooth(geom="line", alpha=0.75, size=0.5, linetype="dashed") +
  scale_alpha_continuous(guide = "none",range = c(.2, .95)) + 
  scale_shape_manual(values=c(15,16,17), name="Exon Size Percentile") +
  scale_fill_manual(values=cbPalette , name="Exon Size Percentile")+
  scale_color_manual(values=cbPalette, name="Exon Size Percentile") +
  scale_x_log10(limits=c(1,1300),expand=c(0.02,0), breaks = c(0,1,11,101,1001),
                labels = c("0","0","10","100","1,000"))+
  scale_y_log10(limits=c(1,1800),expand=c(0.02,0), breaks = c(0,1,11,101,1001),
                labels = c("0","0","10","100","1,000"))+
  xlab("Mean Capture1 DNA Reads per Exon Base Pair")+
  ylab("Mean Capture2 DNA Reads per Exon Base Pair") +
  theme_bw() +
  theme(legend.position = c(0.9,0.15))  

png(filename="Figure3_cap12.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
b1
dev.off()
```

Next step is to find exons with minimum thresholds of gDNA.  These will be our "target" sets along with confidence intervals.  Based on overal DNA coverage, we chose 35X as our "target" set and choose 15X boundaries around that number for confidence intervals.  We will create three `bed` files from our RNA exon coverage stats.

```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")
for i in  "${StringArray[@]}"
do
mawk 'BEGIN { FS = "\t" } ; $4 > 4' ${i}.cov.mean.exon.stats > ${i}.EiRc5.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 9' ${i}.cov.mean.exon.stats > ${i}.EiRc10.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 14' ${i}.cov.mean.exon.stats > ${i}.EiRc15.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 19' ${i}.cov.mean.exon.stats > ${i}.EiRc20.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 34' ${i}.cov.mean.exon.stats > ${i}.EiRc35.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 49' ${i}.cov.mean.exon.stats > ${i}.EiRc50.bed
done

declare -a StringArray=("capture1_filter" "capture2_filter" "capture3_filter" "capture4_filter")
for i in  "${StringArray[@]}"
do
mawk 'BEGIN { FS = "\t" } ; $4 > 4' ${i}.merged.cov.mean.filtered.exon.stats > ${i}.merged.EiRc5.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 9' ${i}.merged.cov.mean.filtered.exon.stats > ${i}.merged.EiRc10.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 14' ${i}.merged.cov.mean.filtered.exon.stats > ${i}.merged.EiRc15.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 19' ${i}.merged.cov.mean.filtered.exon.stats > ${i}.merged.EiRc20.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 34' ${i}.merged.cov.mean.filtered.exon.stats > ${i}.merged.EiRc35.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 49' ${i}.merged.cov.mean.filtered.exon.stats > ${i}.merged.EiRc50.bed
done

mawk 'BEGIN { FS = "\t" } ; $4 > 4' all.DNA.merged.mean.exon.stats > all.filter.merged.EiRc5.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 9' all.DNA.merged.mean.exon.stats > all.filter.merged.EiRc10.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 14' all.DNA.merged.mean.exon.stats > all.filter.merged.EiRc15.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 19' all.DNA.merged.mean.exon.stats > all.filter.merged.EiRc20.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 49' all.DNA.merged.mean.exon.stats > all.filter.merged.EiRc50.bed
mawk 'BEGIN { FS = "\t" } ; $4 > 34' all.DNA.merged.mean.exon.stats > all.filter.merged.EiRc35.bed
```

### Calculating data for table 3

We will use a BASH function to automate this for us:

```bash
counts_per_target(){

#Calculate number of exons with more than 1X coverage
EXONC=$(bedtools coverage -b $WORKING_DIR/03_mapping/$1.F.bam -a $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 5X targets with more than 1X coverage
X5XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/$1.F.bam -a all.filter.merged.EiRc5.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 10X targets with more than 1X coverage
X10XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/$1.F.bam -a all.filter.merged.EiRc10.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 15X targets with more than 1X coverage
X15XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/$1.F.bam -a all.filter.merged.EiRc15.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 20X targets with more than 1X coverage
X20XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/$1.F.bam -a all.filter.merged.EiRc20.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 35X targets with more than 1X coverage
X35XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/$1.F.bam -a all.filter.merged.EiRc35.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 50X targets with more than 1X coverage
X50XC=$(bedtools coverage -b $WORKING_DIR/03_mapping/$1.F.bam -a all.filter.merged.EiRc50.bed -counts -sorted -g $WORKING_DIR/Genome/masked.genome.file | mawk '$4 > 0' | wc -l) 

#Calculate the total number of targets for each set
EXON=$(cat $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed | wc -l )
X5X=$(cat all.filter.merged.EiRc5.bed | wc -l)
X10X=$(cat all.filter.merged.EiRc10.bed | wc -l)
X15X=$(cat all.filter.merged.EiRc15.bed | wc -l)
X20X=$(cat all.filter.merged.EiRc20.bed | wc -l)
X35X=$(cat all.filter.merged.EiRc35.bed | wc -l)
X50X=$(cat all.filter.merged.EiRc50.bed | wc -l)


#Print results in pretty percentages
echo $1
echo `python -c "print(round("$EXONC"/"$EXON" * 100,1))"`"%"
echo `python -c "print(round("$X5XC"/"$X5X" * 100,1))"`"%"
echo `python -c "print(round("$X10XC"/"$X10X" * 100,1))"`"%"
echo `python -c "print(round("$X15XC"/"$X15X" * 100,1))"`"%"
echo `python -c "print(round("$X20XC"/"$X20X" * 100,1))"`"%"
echo `python -c "print(round("$X35XC"/"$X35X" * 100,1))"`"%"
echo `python -c "print(round("$X50XC"/"$X50X" * 100,1))"`"%"
  
}

export -f counts_per_target
```

Now with this function we can use `paste` and subshells to produce the table
```bash
paste <(echo -e "Targets\nAll Exons\n5XR Exons\n10XR Exons\n15XR Exons\n20XR Exons\n35XR Exons\n50XR Exons") <(counts_per_target capture1_filter.merged) <(counts_per_target capture2_filter.merged) <(counts_per_target capture3_filter.merged) <(counts_per_target capture4_filter.merged) > Table3_merged.txt

declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N1" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N1" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N1" "Capture4_N2")
for i in  "${StringArray[@]}"
do
echo "<(counts_per_target $i)"
done

paste <(echo -e "Targets\nAll Exons\n5XR Exons\n10XR Exons\n15XR Exons\n20XR Exons\n35XR Exons\n50XR Exons") <(counts_per_target Capture1_B3) <(counts_per_target Capture1_B4) <(counts_per_target Capture1_G3) <(counts_per_target Capture1_G5) <(counts_per_target Capture1_K3) <(counts_per_target Capture1_K4) <(counts_per_target Capture1_M3) <(counts_per_target Capture1_M4) <(counts_per_target Capture1_N2) <(counts_per_target Capture2_B3) <(counts_per_target Capture2_B4) <(counts_per_target Capture2_G3) <(counts_per_target Capture2_G5) <(counts_per_target Capture2_K3)  <(counts_per_target Capture2_K4) <(counts_per_target Capture2_M3) <(counts_per_target Capture2_M4) <(counts_per_target Capture2_N1) <(counts_per_target Capture2_N2) <(counts_per_target Capture3_B3) <(counts_per_target Capture3_B4) <(counts_per_target Capture3_G3) <(counts_per_target Capture3_G5) <(counts_per_target Capture3_K3) <(counts_per_target Capture3_K4) <(counts_per_target Capture3_M3) <(counts_per_target Capture3_M4) <(counts_per_target Capture3_N1) <(counts_per_target Capture3_N2) <(counts_per_target Capture4_B3) <(counts_per_target Capture4_B4) <(counts_per_target Capture4_G3) <(counts_per_target Capture4_G5) <(counts_per_target Capture4_K3) <(counts_per_target Capture4_K4) <(counts_per_target Capture4_M3) <(counts_per_target Capture4_M4) <(counts_per_target Capture4_N1) > Table3_individual.txt
```

## Generate data for figure 4

Figure 4 is per bp sensitivity looking at coverage across our target sets, near targets (definied as 150 bp around the edge of targets, and off target (everything that is not near or on target).

First steps involve creating our different interval sets using bedtools.


```bash
bedtools flank -i all.filter.merged.EiRc5.bed -b 150 -g $WORKING_DIR/Genome/masked.genome.file | bedtools sort -faidx $WORKING_DIR/Genome/masked.genome.file >  all.filter.merged.EiRc5.150.neartarget.bed 
bedtools slop -i all.filter.merged.EiRc5.bed -b 150 -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc5.150.slop.bed 
bedtools complement -i all.filter.merged.EiRc5.150.slop.bed -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc5.150.offtarget.bed 

bedtools flank -i all.filter.merged.EiRc10.bed -b 150 -g $WORKING_DIR/Genome/masked.genome.file | bedtools sort -faidx $WORKING_DIR/Genome/masked.genome.file >  all.filter.merged.EiRc10.150.neartarget.bed 
bedtools slop -i all.filter.merged.EiRc10.bed -b 150 -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc10.150.slop.bed 
bedtools complement -i all.filter.merged.EiRc10.150.slop.bed -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc10.150.offtarget.bed 

bedtools flank -i all.filter.merged.EiRc15.bed -b 150 -g $WORKING_DIR/Genome/masked.genome.file | bedtools sort -faidx $WORKING_DIR/Genome/masked.genome.file >  all.filter.merged.EiRc15.150.neartarget.bed
bedtools slop -i all.filter.merged.EiRc15.bed -b 150 -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc15.150.slop.bed
bedtools complement -i all.filter.merged.EiRc15.150.slop.bed -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc15.150.offtarget.bed

bedtools flank -i all.filter.merged.EiRc5.bed -b 300 -g $WORKING_DIR/Genome/masked.genome.file | bedtools sort -faidx $WORKING_DIR/Genome/masked.genome.file >  all.filter.merged.EiRc5.300.neartarget.bed 
bedtools slop -i all.filter.merged.EiRc5.bed -b 300 -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc5.300.slop.bed 
bedtools complement -i all.filter.merged.EiRc5.300.slop.bed -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc5.300.offtarget.bed 

bedtools flank -i all.filter.merged.EiRc10.bed -b 300 -g $WORKING_DIR/Genome/masked.genome.file | bedtools sort -faidx $WORKING_DIR/Genome/masked.genome.file >  all.filter.merged.EiRc10.300.neartarget.bed 
bedtools slop -i all.filter.merged.EiRc10.bed -b 300 -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc10.300.slop.bed 
bedtools complement -i all.filter.merged.EiRc10.300.slop.bed -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc10.300.offtarget.bed 

bedtools flank -i all.filter.merged.EiRc15.bed -b 300 -g $WORKING_DIR/Genome/masked.genome.file | bedtools sort -faidx $WORKING_DIR/Genome/masked.genome.file >  all.filter.merged.EiRc15.300.neartarget.bed
bedtools slop -i all.filter.merged.EiRc15.bed -b 300 -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc15.300.slop.bed
bedtools complement -i all.filter.merged.EiRc15.300.slop.bed -g $WORKING_DIR/Genome/masked.genome.file > all.filter.merged.EiRc15.300.offtarget.bed

```

With the target sets defined we again use bedtools to calculate coverage levels across these various genomic regions, and below we use GNU-parallel to speed things up.

```bash
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc5.150.neartarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc5.150.neartarget.all.txt' 
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc5.150.offtarget.bed  -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc5.150.offtarget.all.txt' 
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc5.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc5.150.all.txt' 

ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc10.150.neartarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc10.150.neartarget.all.txt' 
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc10.150.offtarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc10.150.offtarget.all.txt' 
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc10.bed -sorted -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc10.150.all.txt' 

ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc15.150.neartarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc15.150.neartarget.all.txt' 
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc15.150.offtarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc15.150.offtarget.all.txt' 
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc15.bed -sorted -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc15.150.all.txt' 
```

Do the same for the 300 bp segments

```bash
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc5.300.neartarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc5.300.neartarget.all.txt'
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc5.300.offtarget.bed  -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc5.300.offtarget.all.txt'
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc5.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc5.300.all.txt'

ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc10.300.neartarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc10.300.neartarget.all.txt'
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc10.300.offtarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc10.300.offtarget.all.txt'
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc10.bed -sorted -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc10.300.all.txt'

ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc15.300.neartarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc15.300.neartarget.all.txt' 
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc15.300.offtarget.bed -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc15.300.offtarget.all.txt' 
ls $WORKING_DIR/03_mapping/*F.bam | sed 's/.F.bam//g' | parallel 'bedtools coverage -hist -b {}.F.bam -a all.filter.merged.EiRc15.bed -sorted -sorted -g masked.genome.file  | grep ^all > {}.hist.EiRc15.300.all.txt'
```

```r
lsetwd("/home/jgreen/EAGER_OBJ1a/04_coverage_analysis/03_target_interval")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

make_graph <- function(j){
  print(files <- list.files(pattern=paste(j,"_filter.merged.hist.*EiRc5.150*", sep = "")))
  labs <- c("On target","Near target", "Off target")
  
  cov <- list()
  for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])[,c(2,5)]
    cov_cumul=1-cumsum(cov[[i]][,2])
    cov[[i]]$cov_cumul <- c(1,cov_cumul[-length(cov_cumul)])
    cov[[i]]$sample=labs[i]
  }
  
  cov_df=do.call("rbind",cov)
  names(cov_df)[1:2]=c("depth","fraction")
  
  print(files <- list.files(pattern=paste(j,"_filter.merged.hist.*EiRc10.150*", sep = "")))
  cov2 <- list()
  for (i in 1:length(files)) {
    cov2[[i]] <- read.table(files[i])[,c(2,5)]
    cov_cumul2=1-cumsum(cov2[[i]][,2])
    cov2[[i]]$cov_cumul <- c(1,cov_cumul2[-length(cov_cumul2)])
    cov2[[i]]$sample=labs[i]
  }
  cov2_df=do.call("rbind",cov2)
  names(cov2_df)[1:2]=c("depth","fraction")
  
  print(files <- list.files(pattern=paste(j,"_filter.merged.hist.*EiRc15.150*", sep = "")))
  cov3 <- list()
  for (i in 1:length(files)) {
    cov3[[i]] <- read.table(files[i])[,c(2,5)]
    cov_cumul3=1-cumsum(cov3[[i]][,2])
    cov3[[i]]$cov_cumul <- c(1,cov_cumul3[-length(cov_cumul3)])
    cov3[[i]]$sample=labs[i]
  }
  cov3_df=do.call("rbind",cov3)
  names(cov3_df)[1:2]=c("depth","fraction")
  
  cov_df <- subset(cov_df, depth <51)
  cov2_df <- subset(cov2_df, depth <51)
  cov3_df <- subset(cov3_df, depth <51)
  
  cov2_df$high <-cov3_df$cov_cumul
  
  cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7")
  cbPalettep <- c("#0072B2", "#009E73","#D55E00" )
  
  cov_df$sample <- factor(cov_df$sample,levels=c("On target","Near target", "Off target"))
  cov2_df$sample <- factor(cov2_df$sample,levels=c("On target","Near target", "Off target"))
  
  p <- ggplot(cov_df, aes(x= depth, y=cov_cumul, color=sample))  +
    #xlim(0,200)+
    geom_ribbon(data=cov2_df,aes(ymin=cov_cumul,ymax=high, color=sample, fill=sample, alpha=0.4)) +
    scale_alpha(guide = 'none') +
    geom_line(size=1.5)+ 
    scale_color_manual(values=cbPalettep) +
    scale_fill_manual(values=cbPalettep) +
    ylab("% of Target Bases > Depth")+
    #scale_x_log10()+
    xlab("Depth")+
    ggtitle(eval(j)) +
    theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold",hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.position="none")
  #theme(legend.position=c(0.92,0.88))
  
  return(p)
}

sample_names=c("capture1","capture2","capture3","capture4")


Capture1 <- make_graph(sample_names[1])

Capture2 <- make_graph(sample_names[2])
Capture2 <- Capture2 + theme(axis.title.y=element_text(color="transparent"))


Capture3 <- make_graph(sample_names[3])
Capture3 <- Capture3 + theme(axis.title.y=element_text(color="transparent"))


Capture4 <- make_graph(sample_names[4])
Capture4 <- Capture4 + theme(axis.title.y=element_text(color="transparent"))


pdf(file="Figure4.pdf",width=14, height=6.5, bg="transparent")
multiplot(Capture1,Capture3,Capture2,Capture4, cols=2)

dev.off()

pdf(file="Figure4Legend.pdf",width=14, height=6.5, bg="transparent")

Capture4 <- Capture4 + theme(legend.position="bottom")
Capture4
dev.off()
```

## Calculating Specificity

First let's calculate near and off-target intervals for all exons

```bash
bedtools flank -i sorted.ref3.0.exon.sc.hmask.bed -b 150 -g masked.genome.file | bedtools sort -faidx masked.genome.file >  sorted.ref3.0.exon.sc.hmask.neartarget.bed
bedtools slop -i sorted.ref3.0.exon.sc.hmask.bed -b 150 -g masked.genome.file > sorted.ref3.0.exon.sc.hmask.slop.bed
bedtools complement -i sorted.ref3.0.exon.sc.hmask.slop.bed -g masked.genome.file > sorted.ref3.0.exon.sc.hmask.offtarget.bed
```

Now we can create a specificity table for all exons and for expressed targets using a few more BASH functions

```bash
specExon(){

exon_reads=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c -L $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed)
exon_nearr=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c -L $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.neartarget.bed)
exon_nearo=$(samtools view $WORKING_DIR/03_mapping/$1.F.bam  -h -@32 -L $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.bed | samtools view - -@32 -c -L $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.neartarget.bed)
exon_offtr=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c -L $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.offtarget.bed)
exon_nearO=$(samtools view $WORKING_DIR/03_mapping/$1.F.bam  -h -@32 -L $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.slop.bed | samtools view - -@32 -c -L $WORKING_DIR/Genome/sorted.ref3.0.exon.sc.hmask.offtarget.bed)
total=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c)


echo -e $1"\t"`python -c "print(round("$exon_reads"/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"%\t"
  
}

export -f specExon
```

```bash
spec5X(){

exon_reads=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c -L all.filter.merged.EiRc5.bed)
exon_nearr=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c -L all.filter.merged.EiRc5.150.neartarget.bed )
exon_nearo=$(samtools view $WORKING_DIR/03_mapping/$1.F.bam  -h -@32 -L all.filter.merged.EiRc5.bed | samtools view - -@32 -c -L all.filter.merged.EiRc5.150.neartarget.bed )
exon_offtr=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c -L all.filter.merged.EiRc5.150.offtarget.bed)
exon_nearO=$(samtools view $WORKING_DIR/03_mapping/$1.F.bam  -h -@32 -L all.filter.merged.EiRc5.150.neartarget.bed  | samtools view - -@32 -c -L all.filter.merged.EiRc5.150.offtarget.bed)
total=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c)


echo -e `python -c "print(round("$exon_reads"/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"%\t"
  
}

export -f spec5X
```

```bash
spec10X(){

exon_reads=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c -L all.filter.merged.EiRc10.bed)
exon_nearr=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c -L all.filter.merged.EiRc10.150.neartarget.bed )
exon_nearo=$(samtools view $WORKING_DIR/03_mapping/$1.F.bam  -h -@32 -L all.filter.merged.EiRc10.bed | samtools view - -@32 -c -L all.filter.merged.EiRc10.150.neartarget.bed )
exon_offtr=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c -L all.filter.merged.EiRc10.150.offtarget.bed)
exon_nearO=$(samtools view $WORKING_DIR/03_mapping/$1.F.bam  -h -@32 -L all.filter.merged.EiRc10.150.neartarget.bed  | samtools view - -@32 -c -L all.filter.merged.EiRc10.150.offtarget.bed)
total=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c)


echo -e `python -c "print(round("$exon_reads"/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"%\t"
  
}

export -f spec10X
```

```bash
spec15X(){

exon_reads=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c -L all.filter.merged.EiRc15.bed)
exon_nearr=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c -L all.filter.merged.EiRc15.150.neartarget.bed )
exon_nearo=$(samtools view $WORKING_DIR/03_mapping/$1.F.bam  -h -@32 -L all.filter.merged.EiRc15.bed | samtools view - -@32 -c -L all.filter.merged.EiRc15.150.neartarget.bed )
exon_offtr=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c -L all.filter.merged.EiRc15.150.offtarget.bed)
exon_nearO=$(samtools view $WORKING_DIR/03_mapping/$1.F.bam  -h -@32 -L all.filter.merged.EiRc15.150.neartarget.bed  | samtools view - -@32 -c -L all.filter.merged.EiRc15.150.offtarget.bed)
total=$(samtools view -@32 $WORKING_DIR/03_mapping/$1.F.bam -c)


echo -e `python -c "print(round("$exon_reads"/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_nearr" - "$exon_nearo")/"$total" * 100,1))"`"%\t"`python -c "print(round(("$exon_offtr" - "$exon_nearO")/"$total" * 100,1))"`"%\t"
  
}
export -f spec15X
```

Now we use all the functions to create a table
```bash
echo -e "Pool\t%_in_Exons\t%_Near_Exons\t%Off_Target_Exons\t\t%_on_Target5X\t%_Near_Target5X\t%Off_Target5X\t\t%_on_Target10X\t%_Near_Target10X\t%Off_Target10X\t\t%_on_Target15X\t%_Near_Target15X\t%Off_Target15X" > Spec.Table

paste <(specExon capture1_filter.merged) <(spec5X capture1_filter.merged) <(spec10X capture1_filter.merged) <(spec15X capture1_filter.merged) >> Spec.Table
paste <(specExon capture2_filter.merged) <(spec5X capture2_filter.merged) <(spec10X capture2_filter.merged) <(spec15X capture2_filter.merged) >> Spec.Table
paste <(specExon capture3_filter.merged) <(spec5X capture3_filter.merged) <(spec10X capture3_filter.merged) <(spec15X capture3_filter.merged) >> Spec.Table
paste <(specExon capture4_filter.merged) <(spec5X capture4_filter.merged) <(spec10X capture4_filter.merged) <(spec15X capture4_filter.merged) >> Spec.Table
```

|Pool     |%_in_Exons |%_Near_Exons |%Off_Target_Exons |%_on_Target5X |%_Near_Target5X |%Off_Target5X |%_on_Target10X |%_Near_Target10X |%Off_Target10X |%_on_Target15X |%_Near_Target15X |%Off_Target15X |
|:--------|:----------|:------------|:-----------------|:-------------|:---------------|:-------------|:--------------|:----------------|:--------------|:--------------|:----------------|:--------------|
|capture1 |44.1%      |9.1%         |46.8%             |43.5%         |8.7%            |47.8%         |43.3%          |8.6%             |48.1%          |43.0%          |8.5%             |48.5%          |
|capture2 |41.4%      |8.5%         |50.1%             |40.8%         |8.1%            |51.1%         |40.7%          |8.0%             |51.3%          |40.4%          |7.8%             |51.8%          |
|capture3 |41.8%      |6.8%         |51.4%             |41.2%         |6.4%            |52.4%         |41.1%          |6.4%             |52.6%          |40.8%          |6.2%             |52.9%          |
|capture4 |39.6%      |9.3%         |51.1%             |39.1%         |9.0%            |52.0%         |38.9%          |8.8%             |52.3%          |38.5%          |8.7%             |52.8%          |

## Generating data for figure 7

Change back to DNA directory

```bash
cd $WORKING_DIR/Genome
```

The first step is to generate depth per bp data for each capture pool.

```bash
ls *.F.bam | sed 's/.F.bam//g' | parallel "samtools depth -aa {}.F.bam > {}.genome.depth"
```

Next, we extract the region of interest.

```bash
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' capture1_filter.merged.genome.depth | mawk '$2 < 32866205' > capture1_sub1.graph.depth
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' capture2_filter.merged.genome.depth | mawk '$2 < 32866205' > capture2_sub1.graph.depth
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' capture3_filter.merged.genome.depth | mawk '$2 < 32866205' > capture3_sub1.graph.depth
mawk '$1 ~ /NC_035780.1/ && $2 > 32736205' capture4_filter.merged.genome.depth | mawk '$2 < 32866205' > capture4_sub1.graph.depth
```

Next, we add a column for the pool identifier and concatenate into single data table.

```bash
sed -i 's/$/\tcapture1/g' capture1_sub1.graph.depth &
sed -i 's/$/\tcapture2/g' capture2_sub1.graph.depth &
sed -i 's/$/\tcapture3/g' capture3_sub1.graph.depth &
sed -i 's/$/\tcapture4/g' capture4_sub1.graph.depth & 

echo -e "Contig\tbp\tDepth\tSample" > header
cat header capture1_sub1.graph.depth > TotalCovCap1.txt
cat header capture2_sub1.graph.depth > TotalCovCap2.txt
cat header capture3_sub1.graph.depth > TotalCovCap3.txt
cat header capture4_sub1.graph.depth > TotalCovCap4.txt
```

For the graph, we all need the annotations for the gene regions, the exons, and CDS

```bash
mawk '$1 ~ /NC_035780.1/ && $4 > 32736205' ref_C_virginica-3.0_top_level.gff3 | mawk '$5 < 32866205' | mawk '$3 == "exon"' | cut -f1,4,5,9 | uniq -w 30 | sed 's/ID=.*product=//g' | sed 's/;trans.*//g' | sed 's/%.*//g' > exons
cat <(echo -e "Contig\tStart\tEnd\tTreatment") exons > exon.list

mawk '$1 ~ /NC_035780.1/ && $4 > 32736205' ref_C_virginica-3.0_top_level.gff3 | mawk '$5 < 32866205' | mawk '$3 == "mRNA"' | cut -f1,4,5,9 | uniq -w 30 | sed 's/ID=.*product=//g' | sed 's/;trans.*//g' | sed 's/%.*//g' > genes
cat <(echo -e "Contig\tStart\tEnd\tTreatment") genes > genes.list


mawk '$1 ~ /NC_035780.1/ && $4 > 32736205' ref_C_virginica-3.0_top_level.gff3 | mawk '$5 < 32866205' | mawk '$3 == "CDS"' | cut -f1,4,5,9 | uniq -w 30 | sed 's/ID=.*product=//g' | sed 's/;trans.*//g' | sed 's/%.*//g' > CDS
cat <(echo -e "Contig\tStart\tEnd\tTreatment") CDS > CDS.list
```

We also need to perform similar steps for the RNA data


R code for Figure 7

```r
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(scales)
library(zoo)

cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")
DepCap1 <- read.table("TotalCovCap1.txt", header = TRUE)
DepCap1 <- read.table("TotalCovCap2.txt", header = TRUE)
DepCap1 <- read.table("TotalCovCap3.txt", header = TRUE)
DepCap1 <- read.table("TotalCovCap4.txt", header = TRUE)


DepC <- as.data.frame(DepCap1)
DepC$Sample <- factor(DepC$Sample,levels=c("capture1"))
DepR <- as.data.frame(DepCap2)
DepR$Sample <- factor(DepR$Sample,levels=c("capture2"))

exons <- read.table("exon.list", header = TRUE, sep = "\t")
exons <- as.data.frame(exons)

genes <- read.table("genes.list", header = TRUE, sep = "\t")
genes <- as.data.frame(genes)

cds <- read.table("cds.list", header = TRUE, sep = "\t")
cds <- as.data.frame(cds)

subDepC <-subset(DepC, bp <32755000 & bp > 32739000)
subDepR <-subset(DepR, bp <32755000 & bp > 32739000)
subexons <-subset(exons, End <32755205 & End > 32740205)
subgenes <-subset(genes, End <32800757 & Start < 32754201)
subcds <-subset(cds, End <32800757 & Start < 32755000)
subDepR$Depth <- subDepR$Depth / -1
submean.cov <- ddply(subDepC, .(Contig,bp), summarize,  Depth=mean(Depth))
submeanR.cov <- ddply(subDepR, .(Contig,bp), summarize,  Depth=mean(Depth))
subgenes$End[4] <- 32755000


pie(rep(1, length(cbPalette)), labels = sprintf("%d (%s)", seq_along(cbPalette), 
                                                cbPalette), col = cbPalette)
redcol <-"#940000"

cbPalette <- c("#D55E00", "#009E73", "#56B4E9" ,"#0072B2" ,"#E69F00" ,"#F0E442" ,"#999999" ,"#CC79A7","#7570B3")
cbPalettedd <- c( "#009E73","#D55E00", "#E69F00")


dd <- ggplot(subDepC, aes(x= bp, y=Depth)) +
  geom_area(aes(group=Sample),position = "identity",color=alpha("grey30",0.25),fill=cbPalette[4], alpha=0.1, linetype="dotted")+  
  geom_line(data=submean.cov,aes(y=rollmean(Depth, 100, na.pad=TRUE)),colour=cbPalette[4], size =1.0, alpha=0.9)  +
  geom_line(data=submeanR.cov,aes(y=rollmean(Depth, 100, na.pad=TRUE)),colour=redcol, size =1.0, alpha=0.9)  +
  geom_area(data=subDepR, aes(group=Sample),position = "identity",color=alpha("grey30",0.25),fill=redcol, alpha=0.1, linetype="dotted")+
  scale_color_manual(values=cbPalettedd) +
  geom_segment(data=subgenes, aes(x = Start, y = 715, xend = End, yend = 715), size = 6,color=cbPalette[9], alpha=1)+
  geom_segment(data=subexons,aes(x = Start, y = 715, xend = End, yend = 715, color=Treatment),size = 4, alpha=1) +
  geom_segment(data=subcds,aes(x = Start, y = 715, xend = End, yend = 715),size = 1, color="grey90", alpha=1) +
  theme_bw()+
  coord_cartesian(xlim = c(32740000,32755000))+
  xlim(32740000,32755000) +
  scale_y_continuous(limits=c(-415,735),labels=c("250","0","500"), breaks=c(-250,0,500),expand=c(0.01,0)) +
  theme(legend.position="none")

png(filename="Figure7.png", type="cairo",units="px", width=5600, 
    height=3000, res=600, bg="transparent")
dd
dev.off()
```

## Code for SNP calling and statistics

First, use dDocent_ngs to call raw variants

#wget https://raw.githubusercontent.com/jpuritz/EecSeq/master/Bioinformatics/SNPConfig
#bash ./dDocent_ngs.sh SNPConfig

Next decompose raw variants into SNPs and InDel calls

```bash
vcfallelicprimitives TotalRawSNPs.vcf.gz --keep-info --keep-geno > decomp.raw.vcf
```

Next, remove InDel calls

```bash
mawk '$0 ~ /#/ || $0 ~ /TYPE\=snp/' decomp.raw.vcf > snp.raw.vcf
```
Next, filter using VCFtools to generate stats

```bash
#Total SNPs: 
vcftools --vcf snp.raw.vcf 

#Total SNPs with quality score higher than 20: 
vcftools --vcf snp.raw.vcf --minQ 20

#Total Exome SNPS: 
vcftools --vcf snp.raw.vcf --bed sorted.ref3.0.exon.sc.bed --minQ 20

#Exome SNPs with minimum mean 16X coverage: 
vcftools --vcf snp.raw.vcf --minQ 20 --bed sorted.ref3.0.exon.sc.bed --min-meanDP 16

#Exome SNPs with minimum mean 32X coverage: 
vcftools --vcf snp.raw.vcf --minQ 20  --bed sorted.ref3.0.exon.sc.bed --min-meanDP 32

#Exome SNPs with minimum mean 48X coverage: 
vcftools --vcf snp.raw.vcf --minQ 20  --bed sorted.ref3.0.exon.sc.bed --min-meanDP 48

#Exome SNPs with minimum mean 80X coverage: 
vcftools --vcf snp.raw.vcf --minQ 20  --bed sorted.ref3.0.exon.sc.bed --min-meanDP 80

#SNPs with minimum mean 80X coverage outside of exome: 
vcftools --vcf snp.raw.vcf --minQ 20  --exclude-bed sorted.ref3.0.exon.sc.bed --min-meanDP 80
```

## Graveyard

# Used for orginally renaming files for demultiplexing

Rename fq suffix to fastq
```bash
for file in *.fq.gz; do mv -- "$file" "${file%fq.gz}fastq.gz"; done
```

Rename R1 files to R2 to reflect reverse read naming convention
```bash
for file in *R1.fastq.gz; do mv -- "$file" "${file%R1.fastq.gz}R2.fastq.gz"; done
for file in *R1_fastqc.zip; do mv -- "$file" "${file%R1_fastqc.zip}R2_fastqc.zip"; done
for file in *R1_fastqc.html; do mv -- "$file" "${file%R1_fastqc.html}R2_fastqc.html"; done
```

Rename F1 and F2 files to R1 to reflect forward read naming convention

```bash
# For F1 file -> R1 file naming
for file in *F1.fastq.gz; do mv -- "$file" "${file%F1.fastq.gz}R1.fastq.gz"; done
for file in *F1_fastqc.zip; do mv -- "$file" "${file%F1_fastqc.zip}R1_fastqc.zip"; done
for file in *F1_fastqc.html; do mv -- "$file" "${file%F1_fastqc.html}R1_fastqc.html"; done

# For F2 file -> R1 file naming
for file in *F2.fastq.gz; do mv -- "$file" "${file%F2.fastq.gz}R1.fastq.gz"; done
for file in *F2_fastqc.zip; do mv -- "$file" "${file%F2_fastqc.zip}R1_fastqc.zip"; done
for file in *F2_fastqc.html; do mv -- "$file" "${file%F2_fastqc.html}R1_fastqc.html"; done

```


Making directories for all original samples
```bash
#declare -a StringArray=("Capture1_B3" "Capture1_G5" "Capture1_M3" "Capture1_N2" "Capture2_G3" "Capture2_K4" "Capture2_N1" "Capture3_B4" "Capture3_K3" "Capture3_M4" "Capture4_B3" "Capture4_G5" "Capture4_M3" "Capture4_N2" "Capture1_B4" "Capture1_K3" "Capture1_M4" "Capture2_B3" "Capture2_G5" "Capture2_M3" "Capture2_N2" "Capture3_G3" "Capture3_K4" "Capture3_N1" "Capture4_B4" "Capture4_K3" "Capture4_M4" "Capture1_G3" "Capture1_K4" "Capture1_N1" "Capture2_B4" "Capture2_K3" "Capture2_M4" "Capture3_B3" "Capture3_G5" "Capture3_M3" "Capture3_N2" "Capture4_G3" "Capture4_K4" "Capture4_N1")

#for i in "${StringArray[@]}"; do; mkdir ${i}; done
```

## Continuing with global capture files

For right now I am going to keep the capture files as one files. Meaning all individuals will be group into one category based on their size selection of the probe and gDNA fragment length. The Multiqc shows great stats across the board and I will be able to do an initial pass for coverage differences across all the treatment irrespective of individual differences or comparisons. 

## Going to concatenate all the files from each capture into one file
```bash
for i in 1:4
do
cat ../raw/Capture${i}_*.F*.fq.gz >> capture${i}_F.fq.gz
done

for i in 1:4
do
cat ../raw/Capture${i}_*.R*.fq.gz >> capture${i}_R.fq.gz
done

## Run Fastqq and multiqc
fastqc -f fastq -t 20 *.fq.gz
multiqc .
```

## Rename files

```bash
mv capture1_F.fq.gz C1_Sample1.F.fq.gz
mv capture1_R.fq.gz C1_Sample1.R.fq.gz
mv capture2_F.fq.gz C2_Sample1.F.fq.gz
mv capture2_R.fq.gz C2_Sample1.R.fq.gz
mv capture3_F.fq.gz C3_Sample1.F.fq.gz
mv capture3_R.fq.gz C3_Sample1.R.fq.gz
mv capture4_F.fq.gz C4_Sample1.F.fq.gz
mv capture4_R.fq.gz C4_Sample1.R.fq.gz
```

```bash
#CDS
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.0.CDS.sc.hmask.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllCDS.all.split.txt
#Exon
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.0.exon.sc.hmask.bed -g reference.file -sorted -split | grep ^all > C1.hist.AllExon.all.split.txt
#Gene
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.0.gene.sc.hmask.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllGene.all.split.txt
#UTR
bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.0.UTR.sc.hmask.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllUTR.all.split.txt
#Intergenic
#bedtools coverage -hist -b C1_Sample1.F.bam -a sorted.ref3.intergenic.bed -g reference.file -sorted  -split | grep ^all > C1.hist.AllIntergenic.all.split.txt
```

Capture 2
```bash
#CDS
bedtools coverage -hist -b C2_Sample1.F.bam -a sorted.ref3.0.CDS.sc.bed -g reference.file -sorted  -split | grep ^all > C2.hist.AllCDS.all.split.txt
#Exon
bedtools coverage -hist -b C2_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -sorted -split | grep ^all > C2.hist.AllExon.all.split.txt
#Gene
bedtools coverage -hist -b C2_Sample1.F.bam -a sorted.ref3.0.gene.sc.bed -g reference.file -sorted  -split | grep ^all > C2.hist.AllGene.all.split.txt
#UTR
bedtools coverage -hist -b C2_Sample1.F.bam -a sorted.ref3.0.UTR.sc.bed -g reference.file -sorted  -split | grep ^all > C2.hist.allUTR.all.split.txt
#Intergenic
bedtools coverage -hist -b C2_Sample1.F.bam -a sorted.ref3.intergenic.bed -g reference.file -sorted  -split | grep ^all > C2.hist.AllIntergenic.all.split.txt
```

Capture 3
```bash
#CDS
bedtools coverage -hist -b C3_Sample1.F.bam -a sorted.ref3.0.CDS.sc.bed -g reference.file -sorted  -split | grep ^all > C3.hist.AllCDS.all.split.txt
#Exon
bedtools coverage -hist -b C3_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -sorted -split | grep ^all > C3.hist.AllExon.all.split.txt
#Gene
bedtools coverage -hist -b C3_Sample1.F.bam -a sorted.ref3.0.gene.sc.bed -g reference.file -sorted  -split | grep ^all > C3.hist.AllGene.all.split.txt
#UTR
bedtools coverage -hist -b C3_Sample1.F.bam -a sorted.ref3.0.UTR.sc.bed -g reference.file -sorted  -split | grep ^all > C3.hist.allUTR.all.split.txt
#Intergenic
bedtools coverage -hist -b C3_Sample1.F.bam -a sorted.ref3.intergenic.bed -g reference.file -sorted  -split | grep ^all > C3.hist.AllIntergenic.all.split.txt
```

Capture 4
```bash
#CDS
bedtools coverage -hist -b C4_Sample1.F.bam -a sorted.ref3.0.CDS.sc.bed -g reference.file -sorted  -split | grep ^all > C4.hist.AllCDS.all.split.txt
#Exon
bedtools coverage -hist -b C4_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -sorted -split | grep ^all > C4.hist.AllExon.all.split.txt
#Gene
bedtools coverage -hist -b C4_Sample1.F.bam -a sorted.ref3.0.gene.sc.bed -g reference.file -sorted  -split | grep ^all > C4.hist.AllGene.all.split.txt
#UTR
bedtools coverage -hist -b C4_Sample1.F.bam -a sorted.ref3.0.UTR.sc.bed -g reference.file -sorted  -split | grep ^all > C4.hist.allUTR.all.split.txt
#Intergenic
bedtools coverage -hist -b C4_Sample1.F.bam -a sorted.ref3.intergenic.bed -g reference.file -sorted  -split | grep ^all > C4.hist.AllIntergenic.all.split.txt
```

## Generate data for Figure 3

Calculate Exon percentiles

Get total coverage counts per exon

```bash
#Capture 1
bedtools coverage -b C1_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -mean -split > cov.mean.filtered.C1.exon.stats
#Capture 2
bedtools coverage -b C2_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -mean --split> cov.mean.filtered.C2.exon.stats
#Capture 3
bedtools coverage -b C3_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -mean > cov.mean.filtered.C3.exon.stats
#Capture 4
bedtools coverage -b C4_Sample1.F.bam -a sorted.ref3.0.exon.sc.bed -g reference.file -mean > cov.mean.filtered.C4.exon.stats
```

```bash
samtools merge -@20 filter.merged.bam C1_Sample1.F.bam C2_Sample1.F.bam C3_Sample1.F.bam C4_Sample1.F.bam
```

Get total coverage counts per exon
```bash
bedtools coverage -b filter.merged.bam -a sorted.ref3.0.exon.sc.bed -sorted -g reference.file -mean > capture.cov.mean.filtered.merged.exon.stats
```

Not sure what to do with these stats files since we dont have RNA coverage. Think this through. Maybe use percentiles as x - axis

Remove mtDNA from stat files

```bash
mawk '!/NC_007175.2/' capture.cov.mean.filtered.merged.exon.stats > capture.DNA.mean.stats
```

Calculate lower 10th percentile of exon sizes
```bash
mawk '{print $3 -$2}' cov.mean.filtered.C1.exon.stats |sort -g | perl -e '$d=.1;@l=<>;print $l[int($d*@l)]'
```
Resut: `59`

Calculate upper 10th percentile of exon sizes
```bash
mawk '{print $3 -$2}' cov.mean.filtered.C1.exon.stats |sort -g | perl -e '$d=.9;@l=<>;print $l[int($d*@l)]'
```
Result: `517`

Mark exons into size classes based on size distribution and create data table
```bash
mawk '{if ( $3 -$2 > 517 ) print $0 "\tUpper"; else if ( $3 - $2 < 59 ) print $0 "\tLower"; else if ( $3 - $2 > 58 && $3 - $2 < 518) print $0 "\tMiddle" }' cov.mean.filtered.C1.exon.stats > C1.rnd.mean.cov.stats.class
echo -e "Chrom\tStart\tEnd\tDNA_Coverage\tExon_Size_Class" > header
cat header C1.rnd.mean.cov.stats.class > C1.ExonMeanCoverage.txt
```

Next step is to find exons with minimum thresholds of gDNA.  These will be our "target" sets along with confidence intervals.  Based on overal DNA coverage, we chose 35X as our "target" set and choose 15X boundaries around that number for confidence intervals.  We will create three `bed` files from our RNA exon coverage stats.

```bash
mawk 'BEGIN { FS = "\t" } ; $11 > 19' cov.mean.filtered.C1.exon.stats > C1.EiRc20.bed
mawk 'BEGIN { FS = "\t" } ; $11 > 49' cov.mean.filtered.C1.exon.stats > C1.EiRc50.bed
mawk 'BEGIN { FS = "\t" } ; $11 > 34' cov.mean.filtered.C1.exon.stats > C1.EiRc35.bed
```
### Calculating data for table 3

We will use a BASH function to automate this for us:

```bash
counts_per_target(){

#Calculate number of exons with more than 1X coverage
EXONC=$(bedtools coverage -b $1.F.bam -a sorted.ref3.0.exon.sc.bed -counts -sorted -g reference.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 35X targets with more than 1X coverage
X35XC=$(bedtools coverage -b $1.F.bam -a C1.EiRc20.bed -counts -sorted -g reference.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 20X targets with more than 1X coverage
X20XC=$(bedtools coverage -b $1.F.bam -a C1.EiRc50.bed -counts -sorted -g reference.file | mawk '$4 > 0' | wc -l) 
#Calculate number of 50X targets with more than 1X coverage
X50XC=$(bedtools coverage -b $1.F.bam -a C1.EiRc35.bed -counts -sorted -g reference.file | mawk '$4 > 0' | wc -l) 

#Calculate the total number of targets for each set
EXON=$(cat sorted.ref3.0.exon.sc.bed | wc -l )
X35X=$(cat C1.EiRc20.bed | wc -l)
X20X=$(cat C1.EiRc50.bed | wc -l)
X50X=$(cat C1.EiRc35.bed | wc -l)

#Print results in pretty percentages
echo $1
echo `python -c "print(round("$EXONC"/"$EXON" * 100,1))"`"%"
echo `python -c "print(round("$X20XC"/"$X20X" * 100,1))"`"%"
echo `python -c "print(round("$X35XC"/"$X35X" * 100,1))"`"%"
echo `python -c "print(round("$X50XC"/"$X50X" * 100,1))"`"%"
  
}

export -f counts_per_target
```

Now with this function we can use `paste` and subshells to produce the table
```bash
paste <(echo -e "Targets\nAll Exons\n20XR Exons\n35XR Exons\n50XR Exons") <(counts_per_target C1_Sample1) <(counts_per_target C2_Sample1) <(counts_per_target C3_Sample1) <(counts_per_target C4_Sample1) > Table3.txt
```

First I have to create a barcode file that is tab delimited

```bash
nano combined_barcodes.txt
```

```text
ATCGCG	CGTGAT	Capture1
CGATGT	ACATCG	Capture2
TTAGGC	GCCTAA	Capture3
TGGCCA	TGGCCA	Capture4
```

Be careful to create this files in nane and don't edit it in the VScode text editor. VScode will change the endline with a "M^" character that is incompatible with process_shortreads program. To view if there are there use the following

```bash
cat -A <barcodes.tsv>
```

Make directories for all of the files

```bash
mkdir samples
cd samples
```

Making directories for Capture[1-4]_N1 and reseq files
```bash
declare -a StringArray=("Capture1_N1" "Capture1_N1_reseq" "Capture2_N1" "Capture2_N1_reseq" "Capture3_N1" "Capture3_N1_reseq" "Capture4_N1" "Capture4_N1_reseq")

for i in "${StringArray[@]}"
do
mkdir ${i}
done
```

Renaming samples directory to combined barcode samples directory

```bash
mv samples/ combined_barcode_samples
```

Try process_shortreads on one pair of files (ran from 01_process_reads/raw/ directory)

```bash
process_shortreads -P -1 Capture1_N1_R1.fastq.gz -2 Capture1_N1_R2.fastq.gz -o ../demux/combined_barcode_samples/Capture1_N1/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture1_N1.reseq.R1.fastq.gz -2 Capture1_N1.reseq.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture1_N1_reseq/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture2_N1.R1.fastq.gz -2 Capture2_N1.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture2_N1/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture2_N1.reseq.R1.fastq.gz -2 Capture2_N1.reseq.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture2_N1_reseq/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture3_N1.R1.fastq.gz -2 Capture3_N1.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture3_N1/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture3_N1.reseq.R1.fastq.gz -2 Capture3_N1.reseq.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture3_N1_reseq/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture4_N1.R1.fastq.gz -2 Capture4_N1.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture4_N1/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
process_shortreads -P -1 Capture4_N1.reseq.R1.fastq.gz -2 Capture4_N1.reseq.R2.fastq.gz -o ../demux/combined_barcode_samples/Capture4_N1_reseq/ -b ../barcodes/combined_barcode.txt --inline_inline -r -c -q
```

Make smaller process_shortread log files for easier viewing

```bash
declare -a StringArray=("Capture1_N1" "Capture1_N1_reseq" "Capture2_N1" "Capture2_N1_reseq" "Capture3_N1" "Capture3_N1_reseq" "Capture4_N1" "Capture4_N1_reseq")

for i in "${StringArray[@]}"
do
echo ${i}
done

for i in "${StringArray[@]}"
do
head -n 100 ${i}/process_shortreads.log > ${i}/${i}process_shortreads.short.log
done
```

See 01_process_reads/demux/combined_barcode_samples with each directories process_shortreads.short.log file. What we can see are large number of ambiguous reads being discarded amny of which are in the sequences not recorded. 

| File                    | Retained Reads | Low Quality | Ambiguous Barcodes | Trimmed Reads | Orphaned paired-end reads | Total    |
|-------------------------|----------------|-------------|--------------------|---------------|---------------------------|----------|
| Capture1_N1_R1.fastq.gz | 648            | 0           | 29671776           | 0             | 0                         | 29672424 |


| Total Sequences      | 29672424 |
|----------------------|----------|
| Ambiguous Barcodes   | 29671776 |
| Low Quality          | 0        |
| Trimmed Reads        | 0        |
| Orphaned Paired-ends | 0        |
| Retained Reads       | 648      |

| Barcode       | Filename | Total | Retained |
|---------------|----------|-------|----------|
| ATCGCG-CGTGAT | Capture1 | 34    | 34       |
| CGATGT-ACATCG | Capture2 | 310   | 310      |
| TTAGGC-GCCTAA | Capture3 | 122   | 122      |
| TGGCCA-TGGCCA | Capture4 | 182   | 182      |

**Sequences not recorded**

| Barcode       | Total    |
|---------------|----------|
| CGATGT-CGATGT | 14620908 |
| CGATGT-GGGGGG | 32992    |
| CGATGT-AGATGT | 22422    |

This is the case for Capture[1-3]_N1 and their reseq files. Capture4_N1 does not have this issue as its adapter barcode is a palindrome. It seems like the forward and reverse reads contain the barcode not its compliment. Going to test a few things

## Using index_null process short reads to investigate issue with ambiguous reads
New directory for single inline barcodes. I am wanting to test out how we retain reads with single barcodes.

Make new barcode file

``bash
nano single_barcode.txt
```

```text
ATCGCG	Capture1
CGATGT	Capture2
TTAGGC	Capture3
TGGCCA	Capture4
```

```bash
mkdir single_barcode_samples
cd single_barcode_samples
```

Making directories for Capture[1-4]_N1 and reseq files
```bash
declare -a StringArray=("Capture1_N1" "Capture1_N1_reseq" "Capture2_N1" "Capture2_N1_reseq" "Capture3_N1" "Capture3_N1_reseq" "Capture4_N1" "Capture4_N1_reseq")

for i in "${StringArray[@]}"
do
mkdir ${i}
done
```

Try process_shortreads on one pair of files (ran from 01_process_reads/raw/ directory)

```bash
process_shortreads -P -1 Capture1_N1_R1.fastq.gz -2 Capture1_N1_R2.fastq.gz -o ../demux/single_barcode_samples/Capture1_N1/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture1_N1.reseq.R1.fastq.gz -2 Capture1_N1.reseq.R2.fastq.gz -o ../demux/single_barcode_samples/Capture1_N1_reseq/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture2_N1.R1.fastq.gz -2 Capture2_N1.R2.fastq.gz -o ../demux/single_barcode_samples/Capture2_N1/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture2_N1.reseq.R1.fastq.gz -2 Capture2_N1.reseq.R2.fastq.gz -o ../demux/single_barcode_samples/Capture2_N1_reseq/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture3_N1.R1.fastq.gz -2 Capture3_N1.R2.fastq.gz -o ../demux/single_barcode_samples/Capture3_N1/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture3_N1.reseq.R1.fastq.gz -2 Capture3_N1.reseq.R2.fastq.gz -o ../demux/single_barcode_samples/Capture3_N1_reseq/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture4_N1.R1.fastq.gz -2 Capture4_N1.R2.fastq.gz -o ../demux/single_barcode_samples/Capture4_N1/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
#process_shortreads -P -1 Capture4_N1.reseq.R1.fastq.gz -2 Capture4_N1.reseq.R2.fastq.gz -o ../demux/single_barcode_samples/Capture4_N1_reseq/ -b ../barcodes/single_barcode.txt --inline_null -r -c -q
```

Make smaller process_shortread log files for easier viewing (run from single_barcode_sample)

```bash
declare -a StringArray=("Capture1_N1" "Capture1_N1_reseq" "Capture2_N1" "Capture2_N1_reseq" "Capture3_N1" "Capture3_N1_reseq" "Capture4_N1" "Capture4_N1_reseq")

for i in "${StringArray[@]}"
do
echo ${i}
done

for i in "${StringArray[@]}"
do
head -n 100 ${i}/process_shortreads.log > ${i}/${i}process_shortreads.short.log
done
```

See 01_process_reads/demux/single_barcode_samples with each Capture1_N1 process_shortreads.short.log file. We retain more reads! Yay! But we have 14 million ambiguous barcodes? Is it not reading on of the paired files? We can also see the the ambiguous barcodes are not recorded.

| File                    | Retained Reads | Low Quality | Ambiguous Barcodes | Trimmed Reads | Orphaned paired-end reads | Total    |
|-------------------------|----------------|-------------|--------------------|---------------|---------------------------|----------|
| Capture1_N1_R1.fastq.gz | 15129309       | 32341       | 14510774           | 0             | 0                         | 29672424 |


| Total Sequences      | 29672424 |
|----------------------|----------|
| Ambiguous Barcodes   | 14510774 |
| Low Quality          | 32341    |
| Trimmed Reads        | 0        |
| Orphaned Paired-ends | 0        |
| Retained Reads       | 15129309 |

| Barcode | Filename | Total    | Retained |
|---------|----------|----------|----------|
| ATCGCG  | Capture1 | 60438    | 60311    |
| CGATGT  | Capture2 | 14938002 | 14906109 |
| TTAGGC  | Capture3 | 92194    | 92011    |
| TGGCCA  | Capture4 | 71016    | 70878    |


Sequences not recorded
| Barcode | Total |
|---------|-------|
| GTTTGG  | 15032 |
| GTGTTC  | 14364 |


Now lets place exact copies of the barcodes doubling them in the file and try to pull out both the forward and reverse files with that.

```bash
process_shortreads -P -1 Capture1_N1_R1.fastq.gz -2 Capture1_N1_R2.fastq.gz -o ../demux/double_barcode_samples/Capture1_N1/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture1_N1.reseq.R1.fastq.gz -2 Capture1_N1.reseq.R2.fastq.gz -o ../demux/double_barcode_samples/Capture1_N1_reseq/ -b ../barcodes/double_barcode.txt --inline_inline -r
```

Make smaller process_shortread log files for easier viewing (run from single_barcode_sample)

I am now more confused. So here we see a drop in about half of our reads in the ambiguous reads. We have lost a lot of the retain barcodes from the other capture compaired to the single approach.

| File                    | Retained Reads | Low Quality | Ambiguous Barcodes | Trimmed Reads | Orphaned paired-end reads | Total    |
|-------------------------|----------------|-------------|--------------------|---------------|---------------------------|----------|
| Capture1_N1_R1.fastq.gz | 14658692       | 31070       | 14982662           | 0             | 0                         | 29672424 |

| Total Sequences      | 29672424 |
|----------------------|----------|
| Ambiguous Barcodes   | 14982662 |
| Low Quality          | 31070    |
| Trimmed Reads        | 0        |
| Orphaned Paired-ends | 0        |
| Retained Reads       | 14658692 |

| Barcode       | Filename | Total    | Retained |
|---------------|----------|----------|----------|
| ATCGCG-ATCGCG | Capture1 | 34       | 34       |
| CGATGT-CGATAT | Capture2 | 14687394 | 14656325 |
| TTAGGC-TTAGGC | Capture3 | 2152     | 2151     |
| TGGCCA-TGGCCA | Capture4 | 182      | 182      |

Sequences not recorded
| Barcode       | Total |
|---------------|-------|
| CGATGT-GGGGGG | 32992 |
| CGATGT-AGATGT | 22422 |

Looking at files barcodes

```bash
Capture1_N1_R1.fastq.gz
zcat Capture1_N1_R1.fastq.gz | head -100000 | mawk 'NR%4==2' | cut -c1-7 | sort | uniq -c | mawk '$1 > 10' | sort -h
zcat Capture1_N1_R2.fastq.gz | head -100000 | mawk 'NR%4==2' | cut -c1-7 | sort | uniq -c | mawk '$1 > 10' | sort -h
zcat Capture1_N1.reseq.R1.fastq.gz | head -100000 | mawk 'NR%4==2' | cut -c1-7 | sort | uniq -c | mawk '$1 > 10' | sort -h
zcat Capture1_N1.reseq.R2.fastq.gz | head -100000 | mawk 'NR%4==2' | cut -c1-7 | sort | uniq -c | mawk '$1 > 10' | sort -h
```

After meeting with Jon on 3/2/2023 we have a few ideas of what the issue are regarding these samples.

1. Probe contamination. The elongated barcodes that are visualized in the multiqc report show a GT overthang compared to the other "clean" samples.
2. Adapter ligation. The adapters could have been ligated weirdly. But since we found clean samples in those pools this was discarded.
3. Sample contamination. Since we shared a lane with other samples our could've gotten contaminated. Novogene sent us back the samples demultipelxed and may have mistakenly sent back these reads or there was an issue where they were localized to on the sequencer. 

To move forward in processing the reads we are going to do the following:

1. Concatenate all the demultiplexed samples into one read. 
2. Process_shortread on each of the contaminated reads.
  * Then concatenate them into one file.
3. Explore other contamination programs
4. Normalization

I will keep 3 and 4 on the back burner in order to begin ddocent mapping, snp calling, and coverage.

## Decontaminating our Capture[1-4]_N1 reads

Use double barcode file located here: 

> /home/jgreen/EAGER_OBJ1a/01_process_reads/barcodes/double_barcode.txt

Code for process_shortreads on the Capture[1-4]_N1 and reseq files

```bash
process_shortreads -P -1 Capture1_N1_R1.fastq.gz -2 Capture1_N1_R2.fastq.gz -o ../demux/double_barcode_samples/Capture1_N1/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture1_N1.reseq.R1.fastq.gz -2 Capture1_N1.reseq.R2.fastq.gz -o ../demux/double_barcode_samples/Capture1_N1_reseq/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture2_N1.R1.fastq.gz -2 Capture2_N1.R2.fastq.gz -o ../demux/double_barcode_samples/Capture2_N1/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture2_N1.reseq.R1.fastq.gz -2 Capture2_N1.reseq.R2.fastq.gz -o ../demux/double_barcode_samples/Capture2_N1_reseq/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture3_N1.R1.fastq.gz -2 Capture3_N1.R2.fastq.gz -o ../demux/double_barcode_samples/Capture3_N1/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture3_N1.reseq.R1.fastq.gz -2 Capture3_N1.reseq.R2.fastq.gz -o ../demux/double_barcode_samples/Capture3_N1_reseq/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture4_N1.R1.fastq.gz -2 Capture4_N1.R2.fastq.gz -o ../demux/double_barcode_samples/Capture4_N1/ -b ../barcodes/double_barcode.txt --inline_inline -r
process_shortreads -P -1 Capture4_N1.reseq.R1.fastq.gz -2 Capture4_N1.reseq.R2.fastq.gz -o ../demux/double_barcode_samples/Capture4_N1_reseq/ -b ../barcodes/double_barcode.txt --inline_inline -r
```

Names for the files: Capture1_N1.R1.fastq.gz Capture1_N1.R2.fastq.gz Capture2_N1.R1.fastq.gz Capture2_N1.R2.fastq.gz Capture3_N1.R1.fastq.gz Capture3_N1.R2.fastq.gz Capture4_N1.R1.fastq.gz Capture4_N1.R2.fastq.gz

```bash
#Capture1_N1
cat Capture2.1.fq.gz >> ../../../clean/Capture2_N1.F.fq.gz
cat Capture2.2.fq.gz >> ../../../clean/Capture2_N1.R.fq.gz
#Capture1_N1_reseq 
cat Capture2.1.fq.gz >> ../../../clean/Capture2_N1.F.fq.gz
cat Capture2.2.fq.gz >> ../../../clean/Capture2_N1.R.fq.gz
#Capture2_N1
cat Capture2.1.fq.gz >> ../../../clean/Capture2_N1.F.fq.gz
cat Capture2.2.fq.gz >> ../../../clean/Capture2_N1.R.fq.gz
#Capture2_N1_reseq
cat Capture2.1.fq.gz >> ../../../clean/Capture2_N1.F.fq.gz
cat Capture2.2.fq.gz >> ../../../clean/Capture2_N1.R.fq.gz
#Capture3_N1
cat Capture3.1.fq.gz >> ../../../clean/Capture3_N1.F.fq.gz
cat Capture3.2.fq.gz >> ../../../clean/Capture3_N1.R.fq.gz
#Capture3_N1_reseq
cat Capture3.1.fq.gz >> ../../../clean/Capture3_N1.F.fq.gz
cat Capture3.2.fq.gz >> ../../../clean/Capture3_N1.R.fq.gz
#Capture4_N1
cat Capture4.1.fq.gz >> ../../../clean/Capture4_N1.F.fq.gz
cat Capture4.2.fq.gz >> ../../../clean/Capture4_N1.R.fq.gz
#Capture4_N1_reseq
cat Capture4.1.fq.gz >> ../../../clean/Capture4_N1.F.fq.gz
cat Capture4.2.fq.gz >> ../../../clean/Capture4_N1.R.fq.gz
```

## Rename files for ddocent analysis

```bash
rename R1 F *R1.fastq.gz
rename R2 R *R2.fastq.gz
rename fastq fq *.fastq.gz
```

[def]: /home/jgreen/EAGER_OBJ1a/04_coverage_analysis/01_genome_region/Figure2_allcapture_N2.png
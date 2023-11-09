# This file analyzes the raw sequenced files and processes them for analysis

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

## Merging reseq files with original files

```bash

declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N1" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N1" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N1" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N1" "Capture4_N2")

for i in "${StringArray[@]}"
do
cat ${i}.reseq.R1.fastq >> ${i}.R1.fastq
cat ${i}.reseq.R2.fastq >> ${i}.R2.fastq
done

```

## Renaming reads

There ae some issues with the way that CAptures have been labeled and some interindividual issues as well. Globally it looks like Capture 1 should be Capture 2 and Capture 2 should be Capture 1. From out work with the insert sizes we know that when I swapped these the insert sizes incorrectly identified each library, meaning that Capture 1 looked like a 350 bp library and Capture 2 looked like a 150 bp library. It should be the opposite. One way of knowing this was that samples G3 and N2 look like they are in the wrong directories, so I switched them out. In fact those were the only samples to be correctly binned, seeing as they were double switched.

I am going to keep the following code as evidence of what I did to lead to this erroneous result, but in this case I am not carrying out the code until the demultiplexing session.

#################################################################### DO NOT DO #######################################################################

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

########################################################### DO NOT DO ####################################################################

## Processing short reads on files that are not demultiplexed

## Demultiplexing capture files

Use process_shortreads to demultiplex files from raw/ to the demux/ folder. One of the major problems is that we are not sure if novogene was able to separate out the different capture pools specifically for Capture[1-4]_N1 samples and Capture[1-4]_N1 resequenced. We need to create a demultiplexing run that will filter the index and the inline barcode into the proper file. The question here is can I demultiplex just the index and then the adapter inline barcode?

Its important to understand what ![process_shortreads](https://catchenlab.life.illinois.edu/stacks/comp/process_shortreads.php) is doing and how to create the barcode files ![Stacks manual](https://catchenlab.life.illinois.edu/stacks/manual/index.php#clean)


## Using duplicated barcodes and inline_inline to investigate ambigous barcode drop all files

New directory for double inline barcodes. I am wanting to test out how we retain reads with double barcodes.

Make new barcode file

``bash
nano single_barcode.txt
```

```text
ATCACG	Capture1
CGATGT	Capture2
TTAGGC	Capture3
TGGCCA	Capture4
```

```bash
mkdir demultiplex
cd demultiplex
```

Making directories for all Capture files
```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N1" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N1" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N1" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N1" "Capture4_N2")

for i in "${StringArray[@]}"
do
mkdir ${i}
done
```

## Running process_shortreads on all files

Use single barcode file located here: 

> /home/jgreen/EAGER_OBJ1a/01_process_reads/barcodes/double_barcode.txt

Code for process_shortreads on all Captures

```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N1" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N1" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N1" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N1" "Capture4_N2")

for i in "${StringArray[@]}"
do
process_shortreads -P -1 $WORKING_DIR/01_process_reads/demultiplex/${i}.R1.fastq.gz -2 $WORKING_DIR/01_process_reads/demultiplex/${i}.R2.fastq.gz -o $WORKING_DIR/01_process_reads/demultiplex/${i}/ -b $WORKING_DIR/01_process_reads/barcodes/double_barcode.txt --inline_inline -r -c -q -D
done
```

## Make all of the process_shortreads short logs for review

```bash
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")
for i in "${StringArray[@]}"
do
head -n 100 ${i}/process_shortreads.log > ${i}/${i}.process_shortreads.short.log
done
```

Need to pull Capture[1-4] reads into the clean directory

```bash
#Capture1
#Note All Capture 1 are read as Capture 2 in the process_shortreads output. Need to pull Capture2.1.fq.gz and Capture2.2.fq.gz and rename them to the proper Capture category. The only ones to omit from this are G3 and N2.
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture2.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture2.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

declare -a StringArray=("Capture1_G3"  "Capture1_N2")
for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture1.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture1.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

#Capture2
#Note All Capture 2 are read as Capture 2 in the process_shortreads output. Need to pull Capture1.1.fq.gz and Capture1.2.fq.gz and rename them to the proper Capture category. The only ones to omit from this are G3 and N2.
declare -a StringArray=( "Capture2_B3" "Capture2_B4" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture1.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture1.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

declare -a StringArray=("Capture2_G3"  "Capture2_N2")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture2.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture2.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

#Capture3
declare -a StringArray=("Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture3.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture3.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done

#Capture4
declare -a StringArray=("Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")

for i in "${StringArray[@]}"
do
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture4.1.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.F.fq.gz
cat $WORKING_DIR/01_process_reads/demultiplex/${i}/Capture4.2.fq.gz > $WORKING_DIR/01_process_reads/clean/${i}.R.fq.gz
done
```

Using chatGPT we can simplify this code to the following

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
    done
done
```

## Analyze clean files

## Run fastqc and multiqc to visualize files

```bash
fastqc -f fastq -t 20 *.fq.gz
multiqc . -o post_demux_multiqc_report
```

## Merging PE files that have close insert sizes

We will be merging PE FASTQ files that are overlapping into SE files. Using Geneious we found that Capture 2,3, and 4 have read intervals less than 2x their read length which can lead to more overlaps between PE read. 

![Alt text](B3_subsample_insertsize-1.png)

The result will be a PE and SE fastq file which will be mapped to the reference. From here we will merged the PE and SE bam files together and proceed with normalization.

## Merging PE reads from a single individual sample

To merge PE reads we will be using the program [FASTP](https://github.com/OpenGene/fastp)

```bash
#Single file
fastp -m --in1 $WORKING_DIR/01_process_reads/clean/Capture1_B3.F.fq.gz --out1 Capture1_B3.F.fastp.fq.gz --in2 $WORKING_DIR/01_process_reads/clean/Capture1_B3.R.fq.gz --out2 Capture1_B3.R.fastp.fq.gz --merged_out Capture1_B3_merged.fq.gz -D -p

#All files
declare -a StringArray=("Capture1_B3" "Capture1_B4" "Capture1_G3" "Capture1_G5" "Capture1_K3" "Capture1_K4" "Capture1_M3" "Capture1_M4" "Capture1_N2" "Capture2_B3" "Capture2_B4" "Capture2_G3" "Capture2_G5" "Capture2_K3" "Capture2_K4" "Capture2_M3" "Capture2_M4" "Capture2_N2" "Capture3_B3" "Capture3_B4" "Capture3_G3" "Capture3_G5" "Capture3_K3" "Capture3_K4" "Capture3_M3" "Capture3_M4" "Capture3_N2" "Capture4_B3" "Capture4_B4" "Capture4_G3" "Capture4_G5" "Capture4_K3" "Capture4_K4" "Capture4_M3" "Capture4_M4" "Capture4_N2")

for i in "${StringArray[@]}"
do
fastp -m --in1 ../clean/${i}.F.fq.gz --out1 ${i}.R1.fq.gz --in2 ../clean/${i}.R.fq.gz --out2 ${i}.R2.fq.gz --merged_out ${i}.R1.fq.gz -D -p -R "${i}" -j ${i}.json -h ${i}.html
done
```

## Run fastqc and multiqc to visualize files

```bash
fastqc -f fastq -t 20 *.fq.gz
multiqc . -o post_merge_multiqc_report
```

## Create directory in 02_ddocent for single and PE fastp filtered reads
```bash
mkdir 01_PE
mkdir 02_SE
mkdir 04_probes
```

Note: Naming convention for dDocent files is very important
* PE must contain F.fq.gz and R.fq.gz for their identifier with correct prefix names listed in the popmap file
* SE must contain F.fq.gz for their identifier with correct prefix names listed in popmap file.

If you need to rename the filenames

```bash
find . -type f -name "*.R1*" -exec bash -c 'mv "$0" "${0/R1/F}"' {} \;
find . -type f -name "*.R2*" -exec bash -c 'mv "$0" "${0/R2/R}"' {} \;
find . -type f -name "*R1*" -exec bash -c 'mv "$0" "${0/R1/F}"' {} \;
```

## Create a reference.file from reference.fasta that is originally the masked genome fasta

```bash
#Copy genome
cp /RAID_STORAGE2/Shared_Data/Oyster_Genome/masked/masked.cvir.genome.fasta .
#Index genome
samtools faidx masked.cvir.genome.fasta
#Make genome.file
mawk -v OFS='\t' {'print $1,$2'} masked.cvir.genome.fasta.fai > masked.genome.file
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

Use bedtools subtract to remove haplotigs present in the bed file from the other bed files

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
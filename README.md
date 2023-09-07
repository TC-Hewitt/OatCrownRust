# OatCrownRust

## Read mapping statistics and coverage plots

After FASTQ files are mapped to a reference FASTA with an aligner such as [BWA](https://github.com/lh3/bwa) or [HISAT2](https://github.com/DaehwanKimLab/hisat2), simple count and coverage statistics are obtained from the subsequent SAM/BAM files with the following scripts using [SAMtools](https://github.com/samtools/samtools) output piped into AWK operations

**1. for a single BAM, count total number of mapped reads and store in var _map_count_**
```
module load samtools/1.12 #your system might have different environment tool command
map_count=$(samtools view -c -F 260 sample01.bam)
```
>_-c_, counts total n alignments instead of printing each alignment<br />
>_-F 260_, exclude unmapped reads or non-primary alignments ("260" set in _FLAG_ field)
<br />

**2. Get read depth at each ref position in BAM and pipe to AWK to generate coverage stats**
```
samtools depth -a sample01.bam | awk -v OFS='\t' -v sname="sample01" -v nmap="$map_count" \
'{sum+=$3} ; $3>10{c++} END{print sname, nmap, sum/NR, c+0, c*100/NR}' > sample01.covstats.out
```
>_-a_, output all positions (including those with zero depth)<br />
>_-v OFS='\t'_, sets awk output to be tab-delimited<br />
>_-v sname="sample01"_, sets awk var for sample name<br />
>_-v nmap="$map_count"_, sets awk var for _map_count_ determined in previous step<br />
>_{sum+=$3}_, cosecutively add depth value in 3rd field to _sum_ var<br />
>_$3>10{c++}_, adds 1 count to var _c_ if depth value in 3rd field > 10 (i.e. 10X coverage)<br />
>_END{print ...}_, after final line of _samtools depth_ output reached, writes the following to file _sample01.covstats.out_:<br />
>    sample name (_sname_)<br />
>    total mapped reads (_nmap_)<br />
>    mean coverage (_sum/NR_, where _NR_=total n of rows processed =n ref positions)<br />
>    total n positions with coverage over 10X (_c+0_)<br />
>    percent of all positions with coverage over 10X (_c*100/NR_)<br />
<br />

multiple BAM files can be processed in parrallel to save time. For example, below is an array job script for 100 samples using the [Slurm Workload Manager](https://github.com/SchedMD/slurm) on a shared cluster
```
#!/usr/bin/env bash

#SBATCH --job-name=covstats
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mem=2GB
#SBATCH --array=0-99

SAMPLES=( $(cat samplenames.txt) ) #file with 1 sample name per line (sans file exts)

module load samtools/1.12

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    map_count=$(samtools view -c -F 260 ${SAMPLES[$i]}.bam)
    samtools depth -a ${SAMPLES[$i]}.bam | awk -v OFS='\t' -v sname=${SAMPLES[$i]} -v nmap="$map_count" \
    '{sum+=$3} ; $3>10{c++} END{print sname, nmap, sum/NR, c+0, c*100/NR}' > temp${i}
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi
```
>the above can be wrapped in file _covstats.sh_ and run as _>sbatch covstats.sh_ <br />
>the output line of each sample is written to seperate temp files assigned _temp0_ to _temp99_ <br />
<br />

**3. the output stats from each sample can now be collated into a single tab-delim table by creating a header and concatenating with all _temp_ files**
```
echo -e "sample\tmapped_reads\tmean_cov\tbp_over_10X\tpcnt_ref_over_10X" > header
cat header temp{0..99} > combined.covstats.txt
rm header temp{0..99}
```
>header names describe the output fields<br />
>final matrix is written to _combined.covstats.txt_<br />
>temporary _header_ and _temp_ files are deleted<br />
<br />

**4. various stats can then be plotted from the table using a simple R script like below**
```
#!/usr/bin/env Rscript
infile = commandArgs(trailingOnly=TRUE)
stats <- read.table(infile, header = T, sep = "\t")
outname <- strsplit(infile, '.txt') 
png(file=paste0(outname, "_plots.png"), width=2000, height=600)
par(mfrow=c(2,1), mar=c(6,4,4,4))
barplot(stats$mean_cov, names=stats$sample, main="mean cov", ylim=c(0,100), las=2)
barplot(stats$pcnt_ref_over_10X, names=stats$sample, main="perc bp over 10X", ylim=c(0,100), las=2)
dev.off()
```
>the above can be wrapped in file _plot_covstats.R_ and run as _>Rscript plot_covstats.R_ <br />
>this creates two basic barplots, one for mean coverage and another for percent of ref bases over 10X coverage, in a single PNG image<br />
<br />

## Assessing variation with population SNP counts

preamble/description.

1. step 1:

`single line command/script`

2. step 2:
```
multi-line command/script
multi-line command/script
```
>further comments/description.

## Generating chromosome maps with karyoploteR

preamble/description.

1. step 1:

`single line command/script`

2. step 2:
```
multi-line command/script
multi-line command/script
```
>further comments/description.

## Generating synteny plots with gggenomes

preamble/description.

1. step 1:

`single line command/script`

2. step 2:
```
multi-line command/script
multi-line command/script
```
>further comments/description.

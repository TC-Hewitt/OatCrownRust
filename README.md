# OatCrownRust

## Read mapping statistics and coverage plots

After FASTQ files are mapped to a reference FASTA with an aligner such as BWA or HISAT2, simple count and coverage statistics are obtained from the subsequent SAM/BAM files with the following scripts using samtools output piped into AWK operations

1. for a single BAM, count total number of mapped reads and store in var _map_count_
```
module load samtools/1.12 \# your system might have different environment tool command
map_count=$(samtools view -c -F 260 sample01.bam)
```
>-c, counts total n alignments instead of printing each alignment
>-F 260, exclude unmapped reads or non-primary alignments ("260" set in FLAG field)

2. Get read depth at each ref position in BAM and pipe to AWK to generate coverage stats
```
samtools depth -a sample01.bam | awk -v OFS='\t' -v sname="sample01" -v nmap="$map_count" \
'{sum+=$3} ; $3>10{c++} END{print sname, nmap, sum/NR, c+0, c*100/NR}' > sample01.covstats.out
```
>-a, output all positions (including those with zero depth)
>-v OFS='\t', sets awk output to be tab-delimited
>-v sname="sample01", sets awk var for sample name
>-v nmap="$map_count", sets awk var for _map_count_ determined in previous step
>{sum+=$3}, cosecutively add depth value in 3rd field to _sum_ var
>$3>10{c++}, adds 1 count to var _c_ if depth value in 3rd field > 10 (i.e. 10X coverage)
>END{print ...}, after final line of samtools depth output reached, prints: 
>    sample name (sname)
>    total mapped reads (nmap)
>    mean coverage (sum/NR, where NR=total n of rows processed =n ref positions)
>    total n positions with coverage over 10X (c+0)
>    percent of all positions with coverage over 10X (c*100/NR)
>writes to out file _sample01.covstats.out_

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

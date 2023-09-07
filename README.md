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

**2. get read depth at each ref position in BAM and pipe to AWK to generate coverage stats**
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

Here we want to obtain the number of variable SNP sites (VSS) in a populations of crown rust fungi. These are genomic positions at which SNP genotypes are heterogeneous within a population (i.e. not all individuals possess the same genotype) and can therefore be used as a relative measure of overall genotypic diversity when comparing populations (provided each population is called against the same reference genome). Variant calling can be performed independently for each population (represented by collections of fungal isolates) using the same reference with appropriate filters applied. 

In this example, we are comparing a population of North American (US) and South African (SA) isolates. However, variant calling was performed simultaneously on all isolates using [FreeBayes](https://github.com/freebayes/freebayes) to produce a single, multisample VCF containing variant data from 20 US isolates and 29 SA isolates. Extra steps are required to subset the US and SA populations and apply additional filtering before both the total SNP sites (TSS) and VSS of each population are counted.

**1. initial filtering** - the raw VCF file _cr_collection.vcf_ is 1st filtered on quality using [vcffilter](https://github.com/vcflib/vcflib) with various standard params
```
module load vcflib/1.0.1
vcffilter -f "QUAL > 20 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & AC > 0 & DP > 500" \
cr_collection.vcf > cr_collection.QA.vcf
```
>_-f_, specifies a filter to apply to the info fields of records, removes sites which do not pass the filter<br />
>_QUAL > 20_, requires quality score no less than 20<br />
>_QUAL / AO > 10_, requires additional contribution of each observation at least 10 log units (~Q10 per read)<br />
>_SAF > 0 & SAR > 0_, requires reads present on both strands<br />
>_RPR > 1 & RPL > 1_, requires at least two reads “balanced” to each side of the site<br />
>_AC > 0_, requires total number of alt alleles in called genotypes more than 0<br />
>_DP > 500_, requires greater than 500 total read depth (with 49 samples this means a min ~10X per sample on average)<br />
>_> cr_collection.QA.vcf_, write filtered records to new "Quality Assured" file<br />
<br />

**2. feature filtering** - we are only interested in biallelic SNPs likely to be informative, so we will apply additional filtering to the VCF using [BCFtools](https://github.com/samtools/bcftools)
```
module load bcftools/1.15.1
bcftools view -m2 -M2 -v snps -i 'F_MISSING < 0.1' -e 'MAF < 0.05' \
--output-type z -o cr_collection.QA.biaSNPs.vcf.gz cr_collection.QA.vcf
```
>_-m2 -M2_, keep sites with min 2 alleles and max 2 alleles<br />
>_-v snps_, keep SNPs only (omits indels and MNPs)<br />
>_-i 'F_MISSING < 0.1'_, include only sites wih less than 10% missing data<br />
>_-e 'MAF < 0.05'_, exclude sites with a minor allele frequency (MAF) of less than 5% (note: skip MAF filtering if you dont want to discard potential rare/novel genotypes)<br />
>_--output-type z_, output gz compressed vcf<br />
>_-o cr_collection.QA.biaSNPs.vcf.gz_, output file name<br />
>_cr_collection.QA.vcf_, input file name<br />
<br />

**3. subset populations** - we will now subset the US and SA populations from this VCF based on sample names
```
bcftools view -S samplenames_US.txt cr_collection.QA.biaSNPs.vcf.gz > cr_pop_US.vcf
bcftools view -S samplenames_SA.txt cr_collection.QA.biaSNPs.vcf.gz > cr_pop_SA.vcf
```
>_-S samplenames_US.txt_, subset US pop using sample names listed in samplenames_US.txt (1 per line)<br />
>_-S samplenames_SA.txt_, subset SA pop using sample names listed in samplenames_SA.txt<br />

Names must match exactly as present in the VCF file. In case they are not known, it is possible to retrieve them from the VCF with a bash script
`zcat cr_collection.QA.biaSNPs.vcf.gz | egrep "^#" | tail -1`
>_zcat_ will print gz compressed files, otherwise use _cat_ if VCF uncompressed<br />
>_egrep_ will find only header lines (starts with "#")<br />
>_tail -1_ will print out the final header line which should contain the names of the first 8 or 9 default fields followed by the sample fields<br />
<br />

**4. get TSS counts** - we now have seperate US and SA files to get their total SNP count. However, since they were subsetted from a single VCF file, each may contain genotypes that are homozygous reference (ref/ref, RR) in every sample. These can be excluded from the final count
```
bcftools view -i 'GT[*]="alt"' cr_pop_US.vcf > cr_pop_US.noRR.vcf
bcftools view -i 'GT[*]="alt"' cr_pop_SA.vcf > cr_pop_SA.noRR.vcf
```
>_-i ..._, this includes only sites where at least one sample has an alt allele, written to a new VCF

TSS and other stats can now be obtained from these VCFs, contained in _bcftools stats_ output
```
bcftools stats cr_pop_US.noRR.vcf > cr_pop_US.noRR.stats.out
bcftools stats cr_pop_US.noRR.vcf > cr_pop_US.noRR.stats.out
```
alternatively, TSS can be counted directly without writing new files
```
bcftools view -i 'GT[*]="alt"' cr_pop_US.vcf | egrep -v "^#" | wc -l
# TSS for US pop =756,980

bcftools view -i 'GT[*]="alt"' cr_pop_US.vcf | egrep -v "^#" | wc -l
# TSS for SA pop =548,091
```
>RR sites are filtered out as before<br />
>_egrep -v_, exclude lines that start with "#" (header lines)<br />
>_wc -l_, count and print number of output lines (= total n SNP records)<br />
<br />

**5. get VSS counts** - different filter params can be used to obtain the number of variable SNP sites. As before, we can write to new files and then run _bcftools stats_, but here we will print the counts directly like above
```
nUS=20
nSA=29

bcftools view -e "COUNT(GT='RR')=$nUS | COUNT(GT='RA')=$nUS | COUNT(GT='AA')=$nUS | COUNT(GT='Aa')=$nUS" \
cr_pop_US.vcf | egrep -v "^#" | wc -l
# VSS for US pop =652,504

bcftools view -e "COUNT(GT='RR')=$nSA | COUNT(GT='RA')=$nSA | COUNT(GT='AA')=$nSA | COUNT(GT='Aa')=$nSA" \
cr_pop_SA.vcf | egrep -v "^#" | wc -l
# VSS for SA pop =47,369
```
>exact sample sizes of the US and SA populations (as present in VCF subsets) are first assigned to to vars _nUS_ and _nSA_ <br />
>_-e ..._, this excludes sites where a count condition is equal to the sample size (_$nUS_ or _$nSA_), i.e. identical genotype (RR=ref/ref hom, RA=ref/alt het, AA=alt/alt hom, Aa=alt/alt het) in every sample. Thus, what remains are only the variable sites<br />
<br />

To save space, uncompresed leftover VCF files can be bgzipped with _samtools_ then indexed with _bcftools_ <br />
`bgzip -@ 4 file.vcf && bcftools index file.vcf.gz`
>_-@ 4_, num threads to use to perform compression<br />
<br />

For context, the SA population is entirely asexual, lacking any recombination, whereas the US population reproduces both sexually and asexually. Therefore, the US population is expected to be more diverse. This is highly evident in the VSS count, which is ~14 times higher than that of the SA population. TSS counts only indicate SNPs against the reference and are heavily influenced by reference choice. Subtracting the VSS from the TSS in the SA population (548,091-47,369) gives ~500K genotypes that were likely already present at the founding of this population, and diversity within this population is represented by only ~47K SNP sites.<br />
<br />

## Generating chromosome maps with [karyoploteR](https://github.com/bernatgel/karyoploteR)

preamble/description.


## Generating synteny plots with [gggenomes](https://github.com/thackl/gggenomes)

preamble/description.

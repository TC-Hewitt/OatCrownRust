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
>final matrix is written to _combined.covstats.txt_ <br />
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
>the above can be wrapped in file _plot_covstats.R_ and run as _>Rscript plot_covstats.R combined.covstats.txt_ <br />
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

Names must match exactly as present in the VCF file. In case they are not known, it is possible to retrieve them from the VCF with a bash script<br />
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

These are instructions to generate simple chromosome maps annotating regions of interest, similar to what appears in the paper figure. More detailed maps incorporating features such as density plots and coverage graphs are covered in the [karyoploteR tutorial](https://bernatgel.github.io/karyoploter_tutorial/). Here we will make basic plots based on a diploid chromosome (chr) set where A and B haplotypes are separated.

**1. Input for chr lengths** - this is simply a tab delim table of chr names, start and end coords which can be obtaind from the _.fasta.fai_ index of the genome FASTA file<br />

`samtools faidx mygenome.fasta`

the following bash scripts can be used to create a regions table from _mygenome.fasta.fai_ compatible with [GRanges](https://bioconductor.org/packages/3.18/bioc/html/GenomicRanges.html) formatting
```
printf "chr\tstart\tend\n" > mygenome.granges.txt
TAB='\t'
cut -f1,2 mygenome.fasta.fai | sed "s/${TAB}/${TAB}0${TAB}/g" >> mygenome.granges.txt
```
this will look like the following in this example, where haplotypes are inidcated by __A_ and __B_ suffixes.
```
chr start   end
chr1_A  0   8018478
chr2_A  0   6810466
chr3_A  0   6979198
chr4_A  0   7574082
chr5_A  0   6841886
chr6_A  0   5186056
chr1_B  0   7419164
chr2_B  0   6665088
chr3_B  0   6167249
chr4_B  0   8040728
chr5_B  0   7117848
chr6_B  0   5207541
```
<br />

**2. Input for regions of interest** - similarly, we will have a table _loci.granges.txt_ with start and end coords for specific loci, with extra fields for locus name and chosen colours as HEX values
```
chr start   end name    colour
chr1_A  2140000 2200000 loc1A '#61d04f'
chr1_B  1980000 2020000 loc1B '#61d04f'
chr1_A  4010000 4030000 loc2A '#2297e6'
chr1_B  3760000 3770000 loc2B '#2297e6'
chr2_A  4590000 4910000 loc3A '#542a18'
chr2_B  4680000 4990000 loc3B '#542a18'
chr6_A  480000  700000  loc4A '#cd0bbc'
chr6_B  570000  770000  loc4B '#cd0bbc'
chr6_A  2740000 2840000 loc5A '#ffa500'
chr6_B  2860000 2930000 loc5B '#ffa500'
```
<br />

**3. Input for text annotation** - we also need a table _loci.labels.txt_ indicating the position and label for our loci. The positions can be anywhere, but here we are just using the start coordinate of each locus with the colour values matching those of the regions
```
chr pos label   colour
chr1_A  2140000 loc1A    '#61d04f'
chr1_B  1980000 loc1B    '#61d04f'
chr1_A  4010000 loc2A    '#2297e6'
chr1_B  3760000 loc2B    '#2297e6'
chr2_A  4590000 loc3A    '#542a18'
chr2_B  4680000 loc3B    '#542a18'
chr6_A  480000  loc4A    '#cd0bbc'
chr6_B  570000  loc4B    '#cd0bbc'
chr6_A  2740000 loc5A    '#ffa500'
chr6_B  2860000 loc5B    '#ffa500'
```
<br />

**4. Input for other local features** - a previous analysis was performed with [SweeD](https://github.com/alachins/sweed) to identify potential seletive sweep regions. It generated the output _selective_sweep.granges.bed_ which we will also plot
```
#chr start   end
chr1_A  112428  172428
chr1_A  890413  950413
chr1_B  3223380 3292289
chr1_B  4567533 4628523
chr2_A  422444  484441
chr2_B  1038753 1098753
chr2_A  1899788 1961785
chr2_B  3122417 3182417
chr2_A  3752711 3812711
chr2_B  4863466 4923466
chr2_A  5693536 5755533
chr2_B  6070114 6169070
chr5_B  217674  277674
```
<br />

**5. Prepare inputs in R** - now in R, load the requisite libraries and import the data. We will also specify the different haplotypes by chr name.
```
library(karyoplotR)
library(svglite)

mygenome <- toGRanges("mygenome.granges.txt")
loci <- toGRanges("loci.granges.txt")
sweeps <- toGRanges("selective_sweep.granges.bed")
labels <- read.table("loci.labels.txt", header=TRUE)
hapA <- c("chr1_A", "chr2_A", "chr3_A", "chr4_A", "chr5_A", "chr6_A")
hapB <- c("chr1_B", "chr2_B", "chr3_B", "chr4_B", "chr5_B", "chr6_B")
```
<br />

**6. Plot A chrs** - we will create separate plots for the A and B chrs so they can later be arranged side-by-side (without subsetting, all chrs will be plotted as a single stack).
First initialise an output image _chrsA.svg_ using [svglite](https://github.com/r-lib/svglite) to save the plot as a scalable vector graphic (SVG) that can be easily resized and modified in a free SVG editor such as [Inkscape](https://inkscape.org/), then run the plot script
```
svglite("chrsA.svg", width=8, height=6)

kp_hapA <- plotKaryotype(genome=mygenome, chromosomes=hapA)
kpPlotRegions(
    kp_hapA, loci, 
    data.panel="ideogram", 
    r0=NULL, r1=NULL, 
    col=loci$colour, 
    border='#000000', 
    avoid.overlapping = FALSE
    )
kpPlotRegions(
    kp_hapA, sweeps, 
    data.panel=1, 
    col="#c00000", 
    avoid.overlapping = FALSE, 
    r0=-0.6, r1=-0.4)
kpAddBaseNumbers(
    kp_hapA, 
    tick.dist = 1000001, 
    tick.len = 20, 
    tick.col="#000000", 
    cex = 0.01
    )
kpPlotMarkers(
    kp_hapA, 
    chr=labels$chr, 
    x=labels$pos, 
    labels=labels$label, 
    text.orientation = "horizontal", 
    line.color = "#000000", 
    label.color = labels$colour, 
    marker.parts = c(0.1, 0.8, 0.1), 
    r1=0.4, cex=0.8, 
    ignore.chromosome.ends=FALSE)

dev.off()
```
>_kp_hapA <- plotKaryotype(...)_, this creates the plot object _kp_hapA_ based on the chr lengths in _mygenome_ and the subset of chrs listed in _hapA_ <br />
>_kpPlotRegions(kp_hapA, loci, ...)_, this draws the loci regions overlayed on the chrs<br />
>_kpPlotRegions(kp_hapA, sweeps, ...)_, this draws the selective sweep regions below the chrs (coloured red)<br />
>_kpAddBaseNumbers(kp_hapA, ...)_, this adds distance ticks along each chr set 1Mb apart<br />
>_kpPlotMarkers(kp_hapA, chr=labels$chr, x=labels$pos, labels=labels$label, ...)_, this adds text labels to each locus specified by names and coordinates in _labels_ <br />

![chrsA_small](https://github.com/TC-Hewitt/OatCrownRust/assets/33470968/b4bc15d7-0b03-45d6-b128-94a61dfd4e42)

**7. Plot B chrs** - we will now do the same for the B chrs
```
svglite("chrsB.svg", width=8, height=6)

kp_hapB <- plotKaryotype(genome=mygenome, chromosomes=hapB)
kpPlotRegions(
    kp_hapB, loci, 
    data.panel="ideogram", 
    r0=NULL, r1=NULL, 
    col=loci$colour, 
    border='#000000', 
    avoid.overlapping = FALSE)
kpPlotRegions(
    kp_hapB, 
    sweeps, 
    data.panel=1, 
    col="#c00000", 
    avoid.overlapping = FALSE, 
    r0=-0.6, r1=-0.4
    )
kpAddBaseNumbers(
    kp_hapB, 
    tick.dist = 1000001, 
    tick.len = 20, 
    tick.col="#000000", 
    cex = 0.01
    )
kpPlotMarkers(kp_hapB, 
    chr=labels$chr, 
    x=labels$pos, 
    labels=labels$label, 
    text.orientation = "horizontal", 
    line.color = "#000000", 
    label.color = labels$colour, 
    marker.parts = c(0.1, 0.8, 0.1), 
    r1=0.4, cex=0.8, 
    ignore.chromosome.ends=FALSE)

dev.off()
```
![chrsB_small](https://github.com/TC-Hewitt/OatCrownRust/assets/33470968/19912205-27df-4b14-9097-49b358176bd2)

this illustrates that one of the loci (_loc3B_) overlaps a region that has a positive signal for selection<br />
<br />

## Generating synteny plots with [gggenomes](https://github.com/thackl/gggenomes)

Here we are creating a synteny plot of genomic regions from three different strains of the same species. One of the assemblies is fully phased into full-length chromosomes, while the other two are partially phased into contigs only. Through GWAS, we have identified significant marker-trait association regions (MTAR) on chromosomes 6 (A and B). We want to visualise annotated gene content and synteny between orthologous regions from all three assemblies. Here, an MTAR is defined as a genomic interval given by the span of SNP markers above a _p_-value threshold in an association peak, which can be retrieved from the output matrix of GWAS programs such as [TASSEL](https://www.maizegenetics.net/tassel).

![ROI_from_GWAS](https://github.com/TC-Hewitt/OatCrownRust/assets/33470968/857a35bf-c785-4a87-af56-3ca3bff481e2)


In our example, the MTARs on chrs 6 in _mygenome.fasta_ are defined by the following coordinates:
```
chr6_A:480000-700000
chr6_B:570000-770000
```
Contigs from the other two assemblies _12SD80_refassem.fa_ and _12NC29_refassem.fa_ that are orthologous to these MTARs were established based on whole assembly alignments between _mygenome.fasta_, _12SD80_refassem.fa_ and _12NC29_refassem.fa_ using a tool such as [Minimap2](https://github.com/lh3/minimap2) or [D-Genies](https://dgenies.toulouse.inra.fr/) to visualise dotplots. In the below illustration, a contig-level assembly was aligned to a chr-level assembly with the dotplot zoomed-in to region of interest on _chr9_. Here we see that _contig_1_ and _contig_8_ map to this region but wih some gaps remaining.

![dotplot_synteny_juxta](https://github.com/TC-Hewitt/OatCrownRust/assets/33470968/165c0ee1-cd72-40e5-88a8-6ff39d3c8adc)


In our example, the orthologous contigs (orthotigs) in _12SD80_refassem.fa_ and _12NC29_refassem.fa_ have the following sequence IDs:
```
# from 12SD80
000065F_004
000065F
000065F_003

# from 12NC29
000020F
000020F_001
```

We also have previously generated annotations for each assembly in the form of GFF files containing predicted gene models, and for _mygenome.fasta_, an additional GFF file for sequence repeats. We will first prepare all the inputs in UNIX/bash before moving to RStudio to make the synteny plot.<br />
<br />

**1. retrieve MTAR and orthotig sequences**<br />
extract MTAR seqs from _mygenome.fasta_ for each haplotype (A and B)
```
module load samtools/1.12
samtools faidx mygenome.fasta chr6_A:480000-700000 | sed -r 's/:[0-9]+-[0-9]+//g' > MTAR_6A.fa
samtools faidx mygenome.fasta chr6_B:570000-770000 | sed -r 's/:[0-9]+-[0-9]+//g' > MTAR_6B.fa
```
>_sed -r ..._, this removes the coordinates which samtools automatically adds to the output fasta headers (otherwise seq IDs won\'t match GFFs)

now pull out the orthotigs from _12SD80_refassem.fa_ and _12NC29_refassem.fa_. This can be done with a tool such as [_get_contigs.py_](https://github.com/TC-Hewitt/Misc_NGS/blob/master/get_contigs.py) from my other repo [Misc_NGS](https://github.com/TC-Hewitt/Misc_NGS)
```
~/Misc_NGS/get_contigs.py -i 12SD80_refassem.fa -l 000065F_004 000065F 000065F_003 -o MTAR_6_12SD80_orthotigs.fa
~/Misc_NGS/get_contigs.py -i 12NC29_refassem.fa -l 000020F 000020F_001 -o MTAR_6_12NC29_orthotigs.fa
```
>_-i_, input fasta file<br />
>_-l_, space-separated list of seq IDs<br />
>_-o_, output fasta file<br />
<br />

**2. retrieve corresponding annotations for MTAR seqs and orthotigs**<br />

using an AWK script, we will use the same coordinates as previously to pull out annotations for each MTAR from _mygenome.models.gff_ then use _grep_ to keep only _mRNA_ features
```
awk '$1 == "chr6_A" && $4 >= 480000 && $5 <= 700000' mygenome.models.gff | grep mRNA > MTAR_6A.mrna.gff
awk '$1 == "chr6_B" && $4 >= 570000 && $5 <= 770000' mygenome.models.gff | grep mRNA > MTAR_6B.mrna.gff
```

we will do the same for the repeat annotations
```
awk '$1 == "chr6_A" && $4 >= 480000 && $5 <= 700000' mygenome.repeats.gff > MTAR_6A.repeats.gff
awk '$1 == "chr6_B" && $4 >= 570000 && $5 <= 770000' mygenome.repeats.gff > MTAR_6B.repeats.gff
```

now we have annotations for just the MTARs; however, the _start_ and _end_ coordinates of each feature are still relative to the whole chrs 6A and 6B; to modify their positions to be relative to _MTAR_6A.fa_ and _MTAR_6B.fa_ we can use AWK to subtract just the start coordinate of the respective MTAR from both the _start_ and _end_ values of each feature
```
awk -F $'\t' ' { $4 = $4 - 480000; $5 = $5 - 480000; print; } ' OFS=$'\t' MTAR_6A.mrna.gff > MTAR_6.mrna.modpos.gff
awk -F $'\t' ' { $4 = $4 - 570000; $5 = $5 - 570000; print; } ' OFS=$'\t' MTAR_6B.mrna.gff >> MTAR_6.mrna.modpos.gff
```
>here we write the outputs from both haplotypes to the same file _MTAR_6.mrna.modpos.gff_ <br />

```
awk -F $'\t' ' { $4 = $4 - 480000; $5 = $5 - 480000; print; } ' OFS=$'\t' MTAR_6A.repeats.gff > MTAR_6.repeats.modpos.gff
awk -F $'\t' ' { $4 = $4 - 570000; $5 = $5 - 570000; print; } ' OFS=$'\t' MTAR_6B.repeats.gff >> MTAR_6.repeats.modpos.gff
```
>here we write the outputs from both haplotypes to the same file _MTAR_6.repeats.modpos.gff_ <br />
this isn\'t always necessary since the R package _gggenomes_ allows intervals to be custom defined in generating synteny plots; however, to save overhead of importing the entire chr or assembly and GFFs files into RStudio, we have opted to prep the inputs beforehand<br />

likewise, we will retrieve the _mRNA_ annotations for the orthotigs from _12SD80_ and _12NC29_ using a _bash_ script; we do not need to modify the feature positions since they are already relative to the individual contigs
```
# 12SD80
for i in 000065F_004 000065F 000065F_003
do
    grep "^$i\s.*mRNA" 12SD80_refassem.models.gff3 >> MTAR_6_12SD80_orthotigs.mrna.gff
done

#12NC29
for i in 000020F 000020F_001
do
    grep "^$i\s.*mRNA" 12NC29_refassem.models.gff3 >> MTAR_6_12NC29_orthotigs.mrna.gff
done
```
<br />

**3. compile fasta files and annotation files of MTARs and orthotigs for later export to RStudio**<br />
combine MTAR and orthotig fastas to _MTAR_6_allseqs.fa_
```
cat \
MTAR_6A.fa \
MTAR_6B.fa \
MTAR_6_12SD80_orthotigs.fa \
MTAR_6_12NC29_orthotigs.fa \
> MTAR_6_allseqs.fa
```

combine MTAR and orthotig gffs to _MTAR_6_allseqs.mrna.gff_
```
cat \
MTAR_6.mrna.modpos.gff \
MTAR_6_12SD80_orthotigs.mrna.gff \
MTAR_6_12NC29_orthotigs.mrna.gff \
> MTAR_6_allseqs.mrna.gff
```

since we only have repeat annotations for the MTARs, we do not need to combine them with those of the orthotigs (unavailable)<br />
<br />

**4. perform an all-vs-all alignment of _MTAR_6_allseqs.fa_ with [Minimap2](https://github.com/lh3/minimap2)** <br />
the PAF output will be later used for drawing links between regions in the synteny plot
```
module load minimap2/2.24
minimap2 -X MTAR_6_allseqs.fa MTAR_6_allseqs.fa > MTAR_6_allseqs.paf
```
>_-X_, skip self and dual mappings (for the all-vs-all mode)<br />
<br />


**5. generate ortholog clusters of gene models using [OrthoFinder](https://github.com/davidemms/OrthoFinder)** <br />
we will use protein fasta files previously created along with the GFF annotations of the three assemblies; a package such as [AGAT](https://github.com/NBISweden/AGAT) can be used to generate translated protein sequences from a gff and nucleotide fasta<br />

the orthogroups will be later used for colouring genes and drawing links between genes in the synteny plot<br />

first we will retrieve the matching sequences from each protein fasta with [_get_contigs.py_](https://github.com/TC-Hewitt/Misc_NGS/blob/master/get_contigs.py) using _MTAR_6_allseqs.mrna.gff_ as a reference
```
~/Misc_NGS/get_contigs.py -i mygenome.proteins.faa -t MTAR_6_allseqs.mrna.gff -s gene_ -o MTAR_6.proteins.faa
~/Misc_NGS/get_contigs.py -i 12SD80_refassem.proteins.faa -t MTAR_6_allseqs.mrna.gff -s gene_ -o MTAR_6_12SD80_orthotigs.proteins.faa
~/Misc_NGS/get_contigs.py -i 12NC29_refassem.proteins.faa -t MTAR_6_allseqs.mrna.gff -s gene_ -o MTAR_6_12NC29_orthotigs.proteins.faa
```
>_-i_, input fasta file<br />
>_-t_, reference file containing IDs of seqs to retrieve<br />
>_-s gene__, prefix of seq IDs to look for in reference file (i.e. starts with "gene_")<br />
>_-o_, output fasta file<br />
<br />

now we will create a separate directory _translated_ to move our retrieved protein seqs and run _orthofinder_ on
```
mkdir translated
mv *.faa translated/

module load orthofinder/2.5.4
orthofinder -a 4 -f translated/
```
>_-a_, number of parallel analysis threads<br />
>_-f_, dir containing fasta proteomes<br />
<br />

the output we want for our purposes is found in *translated/OrthoFinder/Results_?????/Orthogroups/Orthogroups.txt* where the wildcard "?????" substitutes for the date the orthofinder output was created (e.g. "Feb28"); _Orthogroups.txt_ is not in a tractable format, so we can convert it to long format using [_Orthogroups_txtConvert.py_](https://github.com/TC-Hewitt/OatCrownRust/blob/main/Orthogroups_txtConvert.py) for later use in R
```
./Orthogroups_txtConvert.py -i translated/OrthoFinder/Results_?????/Orthogroups/Orthogroups.txt \
> MTAR_6_allseqs.orthogroups.txt
```
<br />

**6. pull out only secreted subset from orthogroups**<br />
in this example, we are particularly interested in genes that encode secreted proteins and want to highlight these in our synteny plot later; here they are denoted by the string "secreted" in the attribute field (column 9) of the gff files; for example, column 9 of a single gene record in _MTAR_6_allseqs.mrna.gff_ looks like the following:
```
id=gene_1-T1;parent=gene_1;product=secreted protein;
```

of course, this depends on the specific annotation workflow used in assigning attributes to gene models, but for our example, our genes of interest can be easily retrieved with _grep_; as we only need the transcript name (_id=..._), we can strip other columns and attributes to output a simple list of IDs
```
grep "secreted" MTAR_6_allseqs.mrna.gff \
| cut -f 9 | sed -r "s/id=([^;]+).*$/\1/g" \
> MTAR_6_allseqs.secreted.txt
```
>_cut -f 9_, limit to col 9<br />
>_sed -r ..._, regex operation discards all chars except those matching group 1 (any consecutive chars after "id=" except ";" in the parentheses)<br />
<br />

we will then use the IDs to extract the orthogroups from _Orthogroups.txt_ and write to a temp file; this will pull out all the genes of an orthogroup, as long as the group contains at least one secreted member in _MTAR_6_allseqs.secreted.txt_ (this is desirable since we may want to know if the secretion state varies between orthologs of different assemblies/strains)
```
for i in $(cat MTAR_allseqs.secreted.txt) 
do
    grep $i translated/OrthoFinder/Results_?????/Orthogroups/Orthogroups.txt >> temp1.txt 
done
```
given the above, _temp1.txt_ is likely to contain duplicate lines which we can remove with _sort_ and _uniq_ and write to a new temp file
```
sort temp1.txt | uniq > temp2.txt
```
we can now convert this to a desirable format using _Orthogroups_txtConvert.py_ then delete the temp files
```
./Orthogroups_txtConvert.py -i temp2.txt > MTAR_6_allseqs.orthogroups.secreted.txt
rm temp1.txt temp2.txt
```
_MTAR_6_allseqs.orthogroups.secreted.txt_ now contains only proteins that are either secreted or belong to the same orthogroup as a secreted protein<br />
<br />

**7. create zip folder of all input files for export to Rstudio**<br />
now all input files for our synteny plot are ready, we will zip them in _exportR_ for easy transfer to RStudio
```
zip exportR \
MTAR_6_allseqs.fa \
MTAR_6_allseqs.gff \
MTAR_6.repeats.modpos.gff \
MTAR_6_allseqs.paf \
MTAR_6_allseqs.orthogroups.txt \
MTAR_6_allseqs.orthogroups.secreted.txt
```

now in RStudio, we want to draw a synteny plot using _gggenomes_ with the following features:
- sequences arranged in bins based on haplotype
- links bewteen adjacent sequences showing blocks of homology
- directional gene annotations
- stranded repeat blocks
- genes coloured by orthogroup (secreted only)
- lines between adjacent sequences linking genes of the same orthogroup
- gene labels (secreted only) 

the following was conducted in RStudio using R 4.0.2 and gggenomes 0.9.5.9000<br />
<br />

**8. load libraries and read in input files**
```
library(gggenomes)
library(svglite)

#assign inputs
my_seqs <- read_seqs("MTAR_6_allseqs.fa")
my_genes <- read_gff3("MTAR_6_allseqs.gff")
my_repeats <- read_gff("MTAR_6.repeats.modpos.gff")
my_links <- read_paf("MTAR_6_allseqs.paf")
my_clusts <- read.table("MTAR_6_allseqs.orthogroups.txt", header=T, sep="\t")
my_clusts_sec <- read.table("MTAR_6_allseqs.orthogroups.secreted.txt", header=T, sep="\t")

#remove 3 leading zeroes from 12SD80 and 12NC29 derived seq_ids for cleaner labels later on
my_seqs$seq_id <- sub("^000", "", my_seqs$seq_id)
my_genes$seq_id <- sub("^000", "", my_genes$seq_id)
my_links$seq_id <- sub("^000", "", my_links$seq_id)
my_links$seq_id2 <- sub("^000", "", my_links$seq_id2)
```
>we wont end up using _my_clusts_ in our final plot and will only use _my_clusts_sec_ as we want to highlight only secreted genes<br />
<br />

**9. create preliminary plot for assessing sequence arrangements**
```
synplot1 <- gggenomes(
     genes = my_genes,
     seqs = my_seqs,
     feats = NULL,
     links = my_links,
     .id = "file_id",
     spacing = 0.05,
     wrap = NULL,
     adjacent_only = FALSE, #will also plot links b/w non-adjacent seqs
     infer_bin_id = seq_id,
     infer_start = min(start, end),
     infer_end = max(start, end),
     infer_length = max(start, end),
     theme = c("clean", NULL),
     .layout = NULL) + geom_link(alpha=0.1) + geom_seq() + geom_bin_label()
```

we can select a few sequences to plot at a time to help figure out how they should be ordered and binned
```
synplot1 %>% pick("chr6_A", "chr6_B", "065F", "065F_004", "065F_003")
```
for example, this plots the A and B MTARs along with orthotigs from _12SD80_. Homology links are shown between all seqs.

![synplot_exp1](https://github.com/TC-Hewitt/OatCrownRust/assets/33470968/f7c40d1e-dd1d-4da0-bd57-bab18ca2080b)


on first look we can see the plot may benefit from alternative ordering and orientation of seqs. For instance, contigs _065F_003_ and _065F_004_ appear to align side-by-side against the larger _065F_ contig, meaning they should be placed in the same bin (each bin representing an inferred haplotype). By default, all seqs are plotted as a single stack unless they are assigned to bins allowing seqs to be positioned in series<br />
<br />

**10. assign seqs to bins and create a new plot**<br />
first we need to know the lengths of our seqs to enter into the table. These can be found in the _my_seqs_ object created earlier

```
print(my_seqs)
#   seq_id  seq_desc    length
    <chr>   <chr>   <int>
    chr6_A  NA  220001
    chr6_B  NA  200001
    000065F_004 NA  103286
    000065F NA  373524
    000065F_003 NA  237372
    000020F NA  481180
    000020F_001 NA  241894
```
we can then define our bin information in a new table object
```
my_bins <- tribble(
  ~bin_id, ~seq_id, ~length,
  "mygenome_1", "chr6_A", 220001,
  "mygenome_2", "chr6_B", 200001,
  "12SD80_1", "065F", 373524,
  "12SD80_2", "065F_004", 103286,
  "12SD80_2", "065F_003", 237372,  
  "12NC29_1", "020F", 481180,
  "12NC29_2", "020F_001", 241894,
)
```
we may also want to filter out smaller links from our _my_links_ object so they do not clutter the plot. Here we will only keep links of thickness 5kb or more and assign to _my_links_5kb_
```
my_links_5kb <- filter(my_links, map_length >= 5000)
```
we can now create our initial plot substituting in _my_bins_ and _my_links_5kb_
```
synplot2 <- gggenomes(
     genes = my_genes,
     seqs = my_bins, # new seq table
     feats = my_repeats,
     links = my_links_5kb, # new links obj
     .id = "file_id",
     spacing = 0.05,
     wrap = NULL,
     adjacent_only = TRUE, # links will only be shown between adjacent seqs
     infer_bin_id = bin_id, # make sure changed to bin_id
     infer_start = min(start, end),
     infer_end = max(start, end),
     infer_length = max(start, end),
     theme = c("clean", NULL),
     .layout = NULL) +
     geom_bin_label() + 
     geom_seq() + 
     geom_link(fill="lightgrey", color="lightgrey") + 
     geom_seq_label()
```

![synplot_exp2](https://github.com/TC-Hewitt/OatCrownRust/assets/33470968/7bb20abd-d5b1-472f-bafc-8d03fb238211)


our seqs are now ordered the way we want, except their orientations are still mismatched<br />
<br />

**11. flip seqs to match orientations**<br />
we can see the alignment can be fixed by flipping seqs _065F_, _065F_004_ and _065F_003_, which are seq numbers 3, 4 and 5, corresponding to the order they are listed in _my_bins_. We assign the untwisted alignment to a new plot object
```
synplot2_flip <- synplot2 %>% flip_seqs(3, 4, 5)
```

![synplot_exp3](https://github.com/TC-Hewitt/OatCrownRust/assets/33470968/2611676c-0c5f-43e8-82b3-9a712155256b)


it is now looking much better, but seq _020F_ is showing some unaligned flanking regions which we might want to trim<br />
<br />

**12. assign seq boundaries to plot only desired regions**<br />
we can create a table defining the start and end coordinates of each seq region we want to show. Note that coordinates are always based on the original seqs (unflipped), so the previous _synplot2_ can be used a reference when defining region coordinates. The default length bar at the bottom can be used as a guide.
```
my_loci <- tribble(
  ~seq_id, ~start, ~end,
  "chr6_A", 1, 220001,
  "chr6_B", 1, 200001,
  "065F", 1, 373524,
  "065F_004", 1, 103286,
  "065F_003", 1, 237372,  
  "020F", 70000, 450000,
  "020F_001", 1, 241894,
)
```
we can then use the focus function with _my_loci_ to create a new plot object showing only desired regions
```
syplot2_flip_focus <- synplot2_flip %>% focus(.track_id=seqs, .loci=my_loci, .locus_id=seq_id)
```

![synplot_exp4](https://github.com/TC-Hewitt/OatCrownRust/assets/33470968/4304ba3e-0f8c-472d-b8da-7f7e1ce0dff7)


now that we have a bare plot combining our desired ordering, orientation and boundaries, we can start adding features and aesthetics<br />
<br />

**13. modify package functions for custom fill and links**<br />
by default, _gggenomes_ function _geom_gene()_ draws features of type "mRNA" from a gff with a lightened color scheme. However, we want our genes to have full solid colour. To temporarily disable automatic lighten,  open and edit the _geom_gene()_ code in the console using `trace(gggenomes::geom_gene, edit=TRUE)`. Now in the text editor window, change the line `rna_def <- aes(fill = colorspace::lighten(fill, .5), color = colorspace::lighten(colour, .5))` to simply `rna_def <- aes()`<br />

our plot already shows links between regions of homology as connected polygons. We also want to draw links between genes in the the same cluster but as single lines. However, the default behaviour of _geom_link()_ is to draw links as polygons with the same thickness of the gene, which can look cluttered. To instead draw links as single lines, we can define a _geom_link_line()_ function using base _ggplot::geom_segment_
```
geom_link_line <- function(mapping = NULL, data = links(), stat = "identity",
          position = "identity", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
          ...){
  default_aes <- aes(y=y, yend=yend, x=(x+xend)/2, xend=(xmin+xmax)/2)
  mapping <- gggenomes:::aes_intersect(mapping, default_aes)

  layer(geom = GeomSegment, mapping = mapping, data = data, stat = stat, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(na.rm = na.rm, ...))
}
```
we can now include this function when adding features to the plot<br />
<br />

**14. add features and custom aesthetics**<br />
we want to colour only secreted genes according to orthgroup (OG) cluster so first setup a colour palette to use. From [_RColorBrewer_](https://r-graph-gallery.com/38-rcolorbrewers-palettes.html) we'll use _Set1_, which has 9 discrete colours but we can expand it with a colour ramp as _my_clusts_sec_ has 17 OGs
```
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
```
>only using the first 8 colours of _Set1_ since the 9th colour is grey, which we want for _na.values_, i.e. all other genes not in _my_clusts_sec_<br />
<br />

now we will put everything together with previous _synplot2_flip_focus_ in new plot _synplot3_
```
synplot3 <- synplot2_flip_focus %>% add_clusters(my_clusts_sec) +
geom_feat(position=position_strand(), colour="black") + 
geom_gene(aes(fill=cluster_id, colour=cluster_id), size=7) + 
geom_gene_tag(aes(label=parent_ids), vjust = -1, hjust = -0.15, size = 2) + 
geom_link_line(data=links(.track_id=2), color="black", linetype=2) +
scale_fill_manual("cluster", values=getPalette(17), na.value="darkgrey") + 
scale_color_manual("cluster", values=getPalette(17), na.value="darkgrey")
```
>_geom_feat_, these are the repeat regions positioned by +/- strand and coloured black<br />
>_geom_gene_, these are the genes with border and fill according to _cluster_id_ from _my_clusts_sec_ <br />
>_geom_gene_tag_, gene labels from _my_seqs_ using _parent_id_ attribute so trancript number isn't shown (see step 6)<br />
>_geom_link_line_, our custom function (step 13) to link genes with lines instead of polygons<br />
>_scale_fill_manual_, assign fill colour to clusters<br />
>_scale_color_manual_, assign border colour to clusters (matching fill)<br />
<br />

**15. shift bins for optimal display and save to final plot _synplot4_** <br />
by default, the seq bins are left-justified but this can make the links look stretched; we can automatically centre the bins if desired
```
synplot4 <- synplot3 %>% shift(center=TRUE)
```
alternatively, it is easy to manually shift bins horizontally (this may take trial-and-error but plot scale bar can be used as an aid)
```
synplot4 <- synplot3 %>% shift(bins=c(1,2,4,6), by=c(50000, 50000, 25000, 50000))
```

![synplot_exp5_small](https://github.com/TC-Hewitt/OatCrownRust/assets/33470968/a8bb91a2-ff44-4775-8a22-a59f2867fcef)


now we are happy with our synteny plot, we can save it as an SVG by invoking _svglite_ within _ggsave_
```
ggsave(file="MTAR_6_synteny_plot.svg", device=svg, plot=synplot4, width=12, height=6)
```
>in the above plot, only gene labels for the first bin are shown for simplicity, coordinates were manually added to the labels for _chr6_A_, _chr6_B_ and _020F_, some superfluous minor homolgy links were removed for cleaner visuals, and scale bar was reduced to show just one increment. The benefit of saving as SVG is that small adjustments and fixes like these are easy to do in a free SVG editor like [Inkscape](https://inkscape.org/) <br />
<br />

more plotting capabilities can be found in the _gggenomes_ [tutorial](https://thackl.github.io/gggenomes/articles/emales.html) pages

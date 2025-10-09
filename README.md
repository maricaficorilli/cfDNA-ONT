# Fragmenting the Future: Development of a Comprehensive Fragmentomics Pipeline Powered by Nanopore Sequencing Data
After sequencing, the raw signal data were used for basecalling with the Dorado 
basecaller, and the resulting reads were aligned to the human reference genome 
(hg38) using minimap2. Tumor fraction and copy number variations (CNVs) were 
estimated using the HMMcopy and ichorCNA packages from Bioconductor/R. 

For motif analysis, we extracted the first four nucleotides at the 5′ end of each 
aligned read using a method based on that described in Katsman et al., 2022. To 
reduce complexity and identify recurring patterns among samples, Non-negative 
Matrix Factorization (NMF) was applied to the matrix of motif frequencies. This 
approach allowed the extraction of three main profiles, each characterized by a 
specific set of dominant motifs, representing underlying patterns that explain 
the variability observed among the samples.

## Basecalling and Demultiplexing

### Command Basacaling

```bash
/path/to/dorado/bin/dorado basecaller hac \
  --modified-bases 5mCG_5hmCG \
  --emit-moves \
  --device cuda:all \
  --no-trim \
  --reference /path/to/reference/genome.fa \
  /path/to/input_fast5_directory \
  | samtools view -b -o /path/to/output/output_pass.bam


```
### Command Demultiplexing
```bash
dorado demux \
  --output-dir ./bam_demuxed \
  --no-classify \
  output_pass.bam

```

### Reference
- Dorado: [https://github.com/nanoporetech/dorado](https://github.com/nanoporetech/dorado)

### Alignment and Filtering with minimap2 

```bash
# Convert BAM to FASTQ
samtools fastq -@ 4 -O -n /path/to/input.bam | gzip > /path/to/output.fastq.gz

# Align FASTQ to reference genome with minimap2 and filter reads
/path/to/minimap2 -ax map-ont --MD -L /path/to/hg38.p13.mmi /path/to/output.fastq.gz \
  | samtools view -h -q 20 -F 0x4 -F 0x100 -F 0x800 \
  | awk '( $9 < 700 || $1 ~ /^@/ )' \
  | samtools view -bS -o /path/to/output.filtered.bam
```

### Reference
- minimap2: https://github.com/lh3/minimap2
- Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18, September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191  





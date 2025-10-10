#!/usr/bin/env python3
import pysam
from pyfaidx import Fasta
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json

# === INTRO MESSAGE ===
print("🧬 Motif analysis script starting...")
print("This script will extract 4-mer end-motifs from a Nanopore BAM file,")
print("generate a top20 barplot, a pie chart of base composition,")
print("save the top20 motifs with counts and frequencies, and create a log file.\n")

# === INPUT ===
# File paths (customizable)
bam_file = "path/to/sample.bam"        # Path to the BAM file
ref_fasta = "path/to/reference.fasta"  # Path to the reference FASTA
chr_list_file = "path/to/chr_list.txt" # Chromosome list file

# Analysis parameters
sample_name = "sample1"    # Sample name
min_mapq = 20              # Minimum mapping quality
max_read_length = 700      # Maximum read length

# === LOAD CHROMOSOME LIST ===
with open(chr_list_file) as f:
    allowed_chroms = [line.strip() for line in f if line.strip()]

# === OPEN FILES ===
bam = pysam.AlignmentFile(bam_file, "rb")
fasta = Fasta(ref_fasta)

motif_counts = Counter()
total_reads = 0
used_reads = 0
mapq_filtered = 0
clipped_filtered = 0
chrom_filtered = 0
n_filtered = 0
length_filtered = 0

# === EXTRACT 4-MER END MOTIFS ===
for read in bam.fetch():
    total_reads += 1

    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        continue
    if read.mapping_quality < min_mapq:
        mapq_filtered += 1
        continue
    if read.query_length > max_read_length:
        length_filtered += 1
        continue

    chrom = bam.get_reference_name(read.reference_id)
    if chrom not in allowed_chroms:
        chrom_filtered += 1
        continue

    # Strand-specific clipping check
    if read.is_reverse:
        if read.cigartuples and read.cigartuples[-1][0] in (4, 5):
            clipped_filtered += 1
            continue
        start = read.reference_end - 4
        seq = fasta[chrom][start:read.reference_end].seq.upper()
        seq = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
    else:
        if read.cigartuples and read.cigartuples[0][0] in (4, 5):
            clipped_filtered += 1
            continue
        start = read.reference_start
        seq = fasta[chrom][start:start+4].seq.upper()

    if len(seq) != 4 or "N" in seq:
        n_filtered += 1
        continue

    motif_counts[seq] += 1
    used_reads += 1

bam.close()

print(f"✅ Total reads processed: {total_reads}")
print(f"✅ Reads used: {used_reads}\n")

# === FREQUENCIES ===
motif_freqs = {m: (c / used_reads) * 100 for m, c in motif_counts.items()}
bases = ["A", "C", "G", "T"]
ordered_motifs = [m for b in bases for m in sorted(k for k in motif_freqs if k.startswith(b))]
freqs = [motif_freqs.get(m, 0) for m in ordered_motifs]

colors_map = {"A": "#1f77b4", "C": "#d62728", "G": "#2ca02c", "T": "#ff7f0e"}
colors = [colors_map[m[0]] for m in ordered_motifs]

# === TOP20 SELECTION ===
top20 = sorted(motif_counts, key=motif_counts.get, reverse=True)[:20]

# === BARPLOT COMPLETO ===
fig, ax = plt.subplots(figsize=(18, 6))
ax.bar(range(len(ordered_motifs)), freqs, color=colors, width=0.5)
ax.yaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

ymax = max(freqs) if freqs else 1
yticks = np.arange(0, ymax + 0.25, 0.25)
ax.set_yticks(yticks)
ax.set_yticklabels([f"{y:.2f}" for y in yticks], fontsize=9)
ax.set_ylim(0, ymax * 1.08)

ax.set_xticks(range(len(ordered_motifs)))
ax.set_xticklabels(ordered_motifs, rotation=90, fontsize=6)
ax.set_ylabel("Frequency (%)", fontsize=12)

for idx, motif in enumerate(ordered_motifs):
    if motif in top20:
        ax.text(idx, freqs[idx] + 0.02, motif,
                ha="center", va="bottom", fontsize=7, rotation=90)

block_size = len(ordered_motifs) // 4 if ordered_motifs else 1
for i, base in enumerate(bases):
    mid = i * block_size + block_size / 2
    ax.text(mid, ymax * 1.08, f"{base}-end",
            ha="center", va="bottom", fontsize=16,
            fontweight="bold", color=colors_map[base])

plt.tight_layout()
plt.savefig(f"{sample_name}_end_motif_plot.png", dpi=300, bbox_inches="tight")
plt.close()

# === PIE CHART ===
block_perc = {b: sum(motif_freqs[m] for m in motif_freqs if m.startswith(b)) for b in bases}
fig2, ax2 = plt.subplots(figsize=(6, 6))
ax2.pie(block_perc.values(),
        labels=[f"{b}-end" for b in bases],
        autopct="%.1f%%", startangle=90,
        colors=[colors_map[b] for b in bases],
        wedgeprops=dict(edgecolor="white"))
ax2.set_title("End-motif composition by block")
plt.savefig(f"{sample_name}_end_motif_pie.png", dpi=300, bbox_inches="tight")
plt.close()

# === TOP20 BARPLOT ===
df_top20 = pd.DataFrame({
    "Motif": top20,
    "Frequency(%)": [motif_freqs[m] for m in top20],
    "Count": [motif_counts[m] for m in top20]
})
colors_top20 = [colors_map[m[0]] for m in top20]

fig3, ax3 = plt.subplots(figsize=(12, 6))
ax3.bar(range(len(top20)), df_top20["Frequency(%)"], color=colors_top20, width=0.6)
ax3.set_xticks(range(len(top20)))
ax3.set_xticklabels(top20, rotation=90, fontsize=9)
ax3.set_ylabel("Frequency (%)", fontsize=12)
ax3.set_title("Top20 End-motifs")
ax3.yaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
plt.tight_layout()
plt.savefig(f"{sample_name}_top20_barplot.png", dpi=300, bbox_inches="tight")
plt.close()

# === OUTPUT TSV TOP20 ===
df_top20.to_csv(f"{sample_name}_top20.tsv", sep="\t", index=False)

# === LOG FILE ===
log = {
    "total_reads": total_reads,
    "used_reads": used_reads,
    "mapq_filtered": mapq_filtered,
    "clipped_filtered": clipped_filtered,
    "chrom_filtered": chrom_filtered,
    "n_filtered": n_filtered,
    "length_filtered": length_filtered
}
with open(f"{sample_name}.motif.log", "w") as f:
    json.dump(log, f, indent=2)

print("\n✅ Outputs saved:")
print(f"- {sample_name}_end_motif_plot.png")
print(f"- {sample_name}_end_motif_pie.png")
print(f"- {sample_name}_top20_barplot.png")
print(f"- {sample_name}_top20.tsv")
print(f"- {sample_name}.motif.log")

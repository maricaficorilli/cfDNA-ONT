#!/usr/bin/env python3
"""
CancerDetector_for_nanopore.py
------------------------------------------------
End-to-end pipeline for cfDNA tumor fraction estimation from Nanopore sequencing data:

BAM (MM/ML) → marker overlap → readProbFile → CancerDetector.py → tumor_burden.txt → pie chart

Features:
✔ Supports Nanopore BAM files with MM/ML methylation tags
✔ Uses a fixed cancer-specific marker panel (markers.hg38.bed by default)
✔ Robust parsing of malformed or non-numeric MM tag offsets
✔ Optional --verbose flag for detailed diagnostic logging
"""

import re
import sys
import os
import gzip
import pysam
import numpy as np
import pandas as pd
from pybedtools import BedTool
from tqdm import tqdm
import subprocess
import matplotlib.pyplot as plt
import tempfile

# ----------------- ARGUMENT CONFIGURATION -----------------
if len(sys.argv) < 2:
    print("❌ Usage: python CancerDetector_for_nanopore.py <bam_file> [markers_file] [CancerDetector_path] [--verbose]")
    sys.exit(1)

bam_file = sys.argv[1]
markers_file = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].startswith("--") else "markers.hg38.bed"
cancerdetector_path = sys.argv[3] if len(sys.argv) > 3 and not sys.argv[3].startswith("--") else "CancerDetector.py"
verbose = "--verbose" in sys.argv

lambda_value = float(os.getenv("CFTOOLS_LAMBDA", 0.5))
sample_id = os.path.basename(bam_file).replace(".sorted", "").replace(".bam", "")
output_file = f"{sample_id}.readProbFile.txt"
autosomes = {f"chr{i}" for i in range(1, 23)}

print("="*50)
print(f"📘 Sample detected: {sample_id}")
print(f"📂 BAM file: {bam_file}")
print(f"🧬 Marker file: {markers_file}")
print(f"💾 Output: {output_file}")
print(f"🧩 CancerDetector path: {cancerdetector_path}")
print("="*50)

# ----------------- FUNCTION DEFINITIONS -----------------
def parse_MM_field(mm, verbose=False):
    """Parse MM tag into base, strand, modification type, and offsets (robust Nanopore version)."""
    mm = mm.rstrip(';')
    blocks = mm.split(';')
    result = []
    bad_offset_warned = False
    for b in blocks:
        if not b:
            continue
        m = re.match(r"^([ACGT])([+\-])([A-Za-z0-9]+)(?:,)?(.*)$", b)
        if not m:
            continue
        base, strand, modcode, offsets_str = m.groups()
        if modcode.lower() not in {"m", "h", "c", "5mc", "5hmc"}:
            continue
        offsets = []
        if offsets_str:
            for x in offsets_str.split(','):
                try:
                    offsets.append(int(x))
                except ValueError:
                    if not bad_offset_warned:
                        print("⚠️ Invalid offset detected in MM tag (e.g., '?', '*') — skipping such entries.")
                        bad_offset_warned = True
        if verbose:
            print(f"Parsed {len(offsets)} offsets for base {base} ({modcode}).")
        result.append({'base': base, 'strand': strand, 'mod': modcode, 'offsets': offsets})
    return result

def offsets_to_read_positions(offsets):
    if not offsets:
        return [], [], []
    cum = np.cumsum(offsets).tolist()
    return cum, [v + i for i, v in enumerate(cum)], [v + 1 for v in cum]

def map_readpos_to_refpos(read):
    mapping = {}
    for read_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
        if read_pos is not None and ref_pos is not None:
            mapping[read_pos] = ref_pos
    return mapping

# ----------------- LOAD MARKERS -----------------
print(f"📂 Loading markers from {markers_file} ...")
if not os.path.exists(markers_file):
    print(f"❌ Error: Marker file {markers_file} not found. Please provide markers.hg38.bed in the correct path.")
    sys.exit(1)

opener = gzip.open if markers_file.endswith(".gz") else open
with opener(markers_file, "rt") as f:
    first_line = f.readline().strip()
    if first_line.lower().startswith("chr"):
        print("🧾 Header detected in marker file, skipping first line.")
    else:
        f.seek(0)
    markers = BedTool(f).filter(lambda x: x.chrom in autosomes).saveas()

# ----------------- SCAN BAM FILE -----------------
print(f"📘 Opening BAM: {bam_file} ...")
bam = pysam.AlignmentFile(bam_file, "rb")

site_records = []
total_mod_calls = 0
reads_with_mod = 0

print("🔍 Scanning BAM and parsing MM/ML tags...")
for read in tqdm(bam.fetch(until_eof=True)):
    if read.is_unmapped or read.reference_name not in autosomes:
        continue
    if not (read.has_tag("MM") and read.has_tag("ML")):
        continue

    mm_parsed = parse_MM_field(read.get_tag("MM"), verbose=verbose)
    if not mm_parsed:
        continue

    ml_raw = read.get_tag("ML")
    ml_list = list(ml_raw) if isinstance(ml_raw, (list, tuple)) else list(ml_raw)
    readpos_to_ref = map_readpos_to_refpos(read)
    read_seq = read.query_sequence
    ml_index = 0
    read_had_mod = False

    for block in mm_parsed:
        base, offsets = block["base"], block["offsets"]
        if not offsets:
            continue
        cum, cum_plus1, cum_plus1b = offsets_to_read_positions(offsets)
        for candidate in (cum, cum_plus1, cum_plus1b):
            for adjust in (0, -1):
                read_positions = [rp + adjust for rp in candidate]
                ref_positions = [readpos_to_ref.get(rp) for rp in read_positions if rp >= 0]
                if len(ref_positions) >= int(0.6 * len(read_positions)):
                    chosen_refpos = ref_positions
                    break
            else:
                continue
            break
        else:
            pos_list = [i for i, ch in enumerate(read_seq) if ch.upper() == base.upper()]
            chosen_refpos = [readpos_to_ref.get(rp) for rp in pos_list[:len(offsets)]]

        for refpos in chosen_refpos:
            if refpos is None or ml_index >= len(ml_list):
                ml_index += 1
                continue
            p = float(ml_list[ml_index]) / 255.0
            site_records.append([read.reference_name, int(refpos), int(refpos)+1, p])
            total_mod_calls += 1
            ml_index += 1
            read_had_mod = True

    if read_had_mod:
        reads_with_mod += 1

bam.close()
print(f"📊 Reads with modifications: {reads_with_mod:,}")
print(f"🧬 Total mod sites extracted: {total_mod_calls:,}")

if total_mod_calls == 0:
    print("⚠️ No usable ML/MM parsed from BAM. Exiting.")
    sys.exit(1)

# ----------------- INTERSECT WITH MARKERS -----------------
print("📎 Intersecting modification sites with marker panel...")
with tempfile.NamedTemporaryFile(mode="wt", delete=False, suffix=".bed.gz") as tmp:
    with gzip.open(tmp.name, "wt") as gz:
        for chrom, start, end, p in site_records:
            gz.write(f"{chrom}\t{int(start)}\t{int(end)}\t{p:.6f}\n")
    tmp.flush()
    sites_bt = BedTool(tmp.name)

ov = sites_bt.intersect(markers, wa=True, wb=True)
final = []
for rec in ov:
    try:
        marker_id = rec.fields[7]
    except Exception:
        marker_id = rec.fields[-1]
    p = float(rec.fields[3])
    final.append([marker_id, p, 1.0 - p])

if not final:
    print("⚠️ No overlap found between sites and markers. Exiting.")
    sys.exit(1)

df = pd.DataFrame(final, columns=["marker_id", "p_tumor", "p_normal"])
df.to_csv(output_file, sep="\t", index=False)
print(f"✅ Wrote {output_file} with {len(df):,} rows")

# ----------------- RUN CANCERDETECTOR -----------------
print("\n🚀 Running CancerDetector...")
cmd = ["python", cancerdetector_path, output_file, str(lambda_value), ".", sample_id]
subprocess.run(cmd, check=True)

result_file = f"{sample_id}.tumor_burden.txt"
if not os.path.exists(result_file):
    print("⚠️ No tumor_burden.txt found. CancerDetector may have failed.")
    sys.exit(1)

# ----------------- PLOT -----------------
print(f"\n📄 Reading results from {result_file} ...")
res = pd.read_csv(result_file, sep="\t", header=None)
res.columns = ["tumor", "normal"]
tumor, normal = res.loc[0, "tumor"], res.loc[0, "normal"]

print(f"🧬 cfDNA tumor burden: {tumor:.4f}")
print(f"🧬 Normal cfDNA fraction: {normal:.4f}")

fig, ax = plt.subplots(figsize=(6,5), facecolor="white")
labels, values, colors = ["Tumor cfDNA", "Normal cfDNA"], [tumor, normal], ["firebrick", "steelblue"]

wedges, _ = ax.pie(values, colors=colors, startangle=90, wedgeprops=dict(width=0.4, edgecolor="white"))
for wedge, label, val in zip(wedges, labels, values):
    ang = (wedge.theta2 - wedge.theta1)/2. + wedge.theta1
    x, y = np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang))
    ha = "left" if x > 0 else "right"
    ax.annotate(f"{label}\n{val*100:.1f}%", xy=(x,y), xytext=(1.5*np.sign(x),1.2*y),
                ha=ha, va="center", arrowprops=dict(arrowstyle="-", lw=1.0, color="gray"),
                fontsize=12, fontweight="bold")

ax.set_title(f"cfDNA Composition – {sample_id}", fontsize=14, fontweight="bold", pad=20)
ax.legend(labels, loc="center left", bbox_to_anchor=(1.05, 0.5), frameon=False, fontsize=10)
plt.tight_layout()
plot_path = f"{sample_id}_cfDNA_Composition.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"📊 Plot saved as {plot_path}")
print("✅ Pipeline completed successfully.")

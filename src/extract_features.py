import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import Counter

# ==============================
# INPUT / OUTPUT PATHS
# ==============================
INPUT_FASTA = "data/sample.fasta"
OUTPUT_CSV = "results/low_complexity_features.csv"

# ==============================
# PARAMETERS
# ==============================
LOW_COMPLEXITY_WINDOW = 20
LOW_COMPLEXITY_THRESHOLD = 0.6   # fraction of dominant AA
REPEAT_K = 3                     # k-mer size for repeats

# Disorder propensity (simplified, FASTA-only proxy)
DISORDER_PRONE_AA = set("PEDQKRSGA")

# ==============================
# FUNCTIONS
# ==============================
def low_complexity_fraction(seq):
    """
    Low-complexity defined as windows dominated by one AA
    """
    lc_positions = set()

    for i in range(len(seq) - LOW_COMPLEXITY_WINDOW + 1):
        window = seq[i:i + LOW_COMPLEXITY_WINDOW]
        freq = Counter(window)
        if max(freq.values()) / LOW_COMPLEXITY_WINDOW >= LOW_COMPLEXITY_THRESHOLD:
            lc_positions.update(range(i, i + LOW_COMPLEXITY_WINDOW))

    return len(lc_positions) / len(seq) if seq else 0.0

def disorder_scores(seq):
    """
    FASTA-only disorder proxy:
    fraction of disorder-prone residues per position
    """
    scores = [(1 if aa in DISORDER_PRONE_AA else 0) for aa in seq]
    if not scores:
        return 0.0, 0.0
    return np.mean(scores), np.max(scores)

def repeat_density(seq, k=REPEAT_K):
    """
    Density of repeated k-mers
    """
    if len(seq) < k:
        return 0.0

    kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
    counts = Counter(kmers)

    repeated = sum(v for v in counts.values() if v > 1)
    return repeated / len(kmers)

# ==============================
# MAIN EXTRACTION
# ==============================
records = []

for record in SeqIO.parse(INPUT_FASTA, "fasta"):
    seq = str(record.seq).replace("X", "").replace("U", "")
    if len(seq) < 30:
        continue

    lc_frac = round(low_complexity_fraction(seq), 4)
    mean_dis, max_dis = disorder_scores(seq)
    rep_den = round(repeat_density(seq), 4)

    records.append({
        "Protein_ID": record.id,
        "LowComplexity_Fraction": lc_frac,
        "Mean_Disorder_Score": round(mean_dis, 4),
        "Max_Disorder_Score": round(max_dis, 4),
        "Repeat_Density": rep_den
    })

# ==============================
# EXPORT
# ==============================
df = pd.DataFrame(records)
df.to_csv(OUTPUT_CSV, index=False)

print(f"✅ Low-complexity & repeat features extracted for {len(df)} proteins")
print(f"📁 Output saved to: {OUTPUT_CSV}")

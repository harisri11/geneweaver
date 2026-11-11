"""Generate a mock genome FASTA file for testing GeneWeaver.

Usage:
    python data/make_mock_genome.py --bases 1000000 --out data/mock_1M.fasta
"""
import argparse
import os
import random


def make(path, n_bases, seed=42):
    random.seed(seed)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as f:
        f.write(">chr_mock mock genome\n")
        written = 0
        while written < n_bases:
            n = min(80, n_bases - written)
            f.write("".join(random.choices("ACGT", k=n)) + "\n")
            written += n


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--bases", type=int, default=1_000_000)
    ap.add_argument("--out", default="data/mock_1M.fasta")
    a = ap.parse_args()
    make(a.out, a.bases)
    print(f"Wrote {a.bases:,} bases to {a.out}")

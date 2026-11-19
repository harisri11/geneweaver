# geneweaver/__main__.py
"""Command-line entry point.

Usage:
    python -m geneweaver --fasta data/mock_10000000.fasta \
        --guide GACGCATAAAGATGAGACGC --backend cpu
"""
import argparse
import time

from .pipeline import run_search


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        prog="geneweaver",
        description="Search a genome FASTA file for CRISPR off-target sites.")
    ap.add_argument("--fasta", required=True, help="path to a .fasta genome file")
    ap.add_argument("--guide", required=True, help="20-base guide sequence (A/C/G/T)")
    ap.add_argument("--max-mm", type=int, default=3, help="max mismatches (default 3)")
    ap.add_argument("--backend", choices=["cpu", "gpu"], default="cpu")
    ap.add_argument("--top", type=int, default=20, help="how many hits to print")
    a = ap.parse_args(argv)

    a.guide = a.guide.upper()
    if len(a.guide) != 20 or set(a.guide) - set("ACGT"):
        ap.error("guide must be exactly 20 bases using only A, C, G, T")
    if a.max_mm < 0:
        ap.error("--max-mm must be 0 or more")
    return a


def main(argv=None):
    a = parse_args(argv)
    t0 = time.perf_counter()
    hits = run_search(a.fasta, a.guide, a.max_mm, a.backend)
    secs = time.perf_counter() - t0

    print(f"{len(hits)} hit(s) in {secs:.2f} s (backend={a.backend})")
    if hits:
        print(f"{'chrom':<12}{'position':>12}{'mm':>4}{'score':>8}  target")
        for h in hits[:a.top]:
            print(f"{h.chrom:<12}{h.pos:>12}{h.mismatches:>4}{h.score:>8.1f}  {h.target}")
    return hits


if __name__ == "__main__":
    main()
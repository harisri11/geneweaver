# benchmarks/benchmark.py
"""CPU baseline benchmark for GeneWeaver.

Times the pure Python search on mock genomes of increasing size and
saves the results to benchmarks/results_cpu.csv.

Usage (from the project root):
    python benchmarks/benchmark.py
    python benchmarks/benchmark.py --sizes 1000000 5000000 10000000
"""
import argparse
import csv
import importlib.util
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from geneweaver.cpu_baseline import cpu_search  # noqa: E402
from geneweaver.parser import encode, iter_chunks  # noqa: E402

GUIDE = "GACGCATAAAGATGAGACGC"


def _load_generator():
    path = os.path.join(ROOT, "data", "make_mock_genome.py")
    spec = importlib.util.spec_from_file_location("make_mock_genome", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.make


def ensure_fasta(n_bases):
    path = os.path.join(ROOT, "data", f"mock_{n_bases}.fasta")
    if not os.path.exists(path):
        print(f"Generating {n_bases:,} bases -> {path}")
        _load_generator()(path, n_bases)
    return path


def time_cpu(fasta, max_mm=3, chunk_size=50_000_000):
    guide = encode(GUIDE)
    t0 = time.perf_counter()
    n_hits = 0
    for _chrom, _start, chunk in iter_chunks(fasta, chunk_size):
        n_hits += len(cpu_search(chunk, guide, max_mm))
    return time.perf_counter() - t0, n_hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=int, nargs="+",
                    default=[1_000_000, 5_000_000, 10_000_000])
    ap.add_argument("--max-mm", type=int, default=3)
    ap.add_argument("--out", default=os.path.join(ROOT, "benchmarks", "results_cpu.csv"))
    a = ap.parse_args()

    rows = []
    for n in a.sizes:
        fasta = ensure_fasta(n)
        secs, hits = time_cpu(fasta, a.max_mm)
        rows.append({"bases": n, "cpu_seconds": round(secs, 3),
                     "mbases_per_sec": round(n / secs / 1e6, 3), "hits": hits})
        print(f"{n:>12,} bases | CPU {secs:8.2f} s | {n / secs / 1e6:6.2f} Mbases/s | hits {hits}")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["bases", "cpu_seconds", "mbases_per_sec", "hits"])
        w.writeheader()
        w.writerows(rows)
    print(f"Saved {a.out}")


if __name__ == "__main__":
    main()
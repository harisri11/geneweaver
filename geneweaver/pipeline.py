"""The pipeline that connects parser -> search backend -> scoring."""
import math
import os
import time
from dataclasses import dataclass

from .cpu_baseline import cpu_search
from .parser import encode, iter_chunks
from .scoring import off_target_score

DECODE = "ACGTN"
CHUNK_SIZE = 50_000_000


@dataclass
class Hit:
    chrom: str
    pos: int
    mismatches: int
    target: str
    score: float


def _get_backend(name):
    if name == "cpu":
        return cpu_search
    if name == "gpu":
        from .gpu_kernel import gpu_search  # imported lazily: CPU-only machines still work
        return gpu_search
    raise ValueError(f"unknown backend: {name}")


def estimate_chunks(fasta, chunk_size):
    est_bases = os.path.getsize(fasta) * 80 // 81
    return max(1, math.ceil(est_bases / chunk_size))


def run_search(fasta, guide_str, max_mm=3, backend="cpu",
               chunk_size=CHUNK_SIZE, on_progress=None):
    """Search a FASTA file and return Hit objects sorted by score (high first)."""
    search = _get_backend(backend)
    guide = encode(guide_str)
    total = estimate_chunks(fasta, chunk_size)
    hits, t0 = [], time.perf_counter()

    for k, (chrom, start, chunk) in enumerate(iter_chunks(fasta, chunk_size)):
        new_hits = []
        for p, mm in search(chunk, guide, max_mm):
            if p >= chunk_size:  # overlap zone: the next chunk reports it
                continue
            target = "".join(DECODE[b] for b in chunk[p:p + 20])
            new_hits.append(Hit(chrom, start + p, mm, target,
                                off_target_score(guide_str, target)))
        hits.extend(new_hits)
        if on_progress:
            on_progress({"done": k + 1, "total": total, "new_hits": new_hits,
                         "elapsed": time.perf_counter() - t0})

    hits.sort(key=lambda h: h.score, reverse=True)
    return hits

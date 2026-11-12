"""FASTA parsing, base encoding and overlapping chunking."""
import numpy as np
from Bio import SeqIO

# A=0, C=1, G=2, T=3, anything else (N, etc.) = 4
LUT = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    LUT[_b] = _i
    LUT[_b + 32] = _i  # lowercase


def encode(seq: str) -> np.ndarray:
    """Convert a DNA string to a uint8 array of codes 0-4."""
    return LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]


def iter_chunks(fasta_path, chunk_size=50_000_000, overlap=22):
    """Yield (chrom, global_start, uint8 array).

    overlap = guide(20) + PAM(3) - 1 = 22, so a match that straddles a
    chunk boundary is still fully contained in one of the chunks.
    """
    for rec in SeqIO.parse(fasta_path, "fasta"):
        arr = encode(str(rec.seq))
        start = 0
        while start < len(arr):
            end = min(start + chunk_size + overlap, len(arr))
            yield rec.id, start, arr[start:end]
            start += chunk_size

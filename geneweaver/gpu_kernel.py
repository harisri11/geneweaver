"""GPU (Numba CUDA) CRISPR off-target search.

One GPU thread checks one genome position: it counts mismatches between
the 20-base guide and the genome window, and checks the NGG PAM.
Returns the same (position, mismatches) pairs as cpu_search.
"""
import numpy as np
from numba import cuda

NO_HIT = 255  # sentinel: position has no valid PAM


@cuda.jit
def align_kernel(genome, guide, out_mm):
    i = cuda.grid(1)
    m = guide.shape[0]
    n = genome.shape[0]
    if i < n - m - 2:
        if genome[i + m + 1] == 2 and genome[i + m + 2] == 2:  # NGG PAM
            mm = 0
            for j in range(m):
                if genome[i + j] != guide[j]:
                    mm += 1
            out_mm[i] = mm
        else:
            out_mm[i] = NO_HIT
    elif i < n:
        out_mm[i] = NO_HIT  # last bases: no room for a full window + PAM


def gpu_search(chunk, guide, max_mm=3, threads=256):
    """Return a list of (position, mismatches) for every hit in chunk."""
    n = chunk.shape[0]
    if n < guide.shape[0] + 3:
        return []
    d_genome = cuda.to_device(chunk)                 # Host -> Device
    d_guide = cuda.to_device(guide)
    d_out = cuda.device_array(n, dtype=np.uint8)
    blocks = (n + threads - 1) // threads
    align_kernel[blocks, threads](d_genome, d_guide, d_out)
    cuda.synchronize()
    out = d_out.copy_to_host()                       # Device -> Host
    idx = np.nonzero(out <= max_mm)[0]
    return [(int(i), int(out[i])) for i in idx]
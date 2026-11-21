"""GPU results must match the CPU baseline exactly.

On a laptop without an NVIDIA GPU, run with the CUDA simulator:
    PowerShell:  $env:NUMBA_ENABLE_CUDASIM = "1"; python -m pytest -q
(The simulator is slow, so these tests use small genomes.)
"""
import random

import pytest
from numba import cuda

from geneweaver.cpu_baseline import cpu_search
from geneweaver.parser import encode

pytestmark = pytest.mark.skipif(not cuda.is_available(),
                                reason="needs an NVIDIA GPU or NUMBA_ENABLE_CUDASIM=1")

GUIDE = "GACGCATAAAGATGAGACGC"


def planted_genome(n=3000, seed=2):
    random.seed(seed)
    seq = list("".join(random.choices("ACGT", k=n)))
    for pos, mm_idx in [(100, []), (700, [5]), (1500, [5, 12]), (n - 23, [8])]:
        t = list(GUIDE)
        for i in mm_idx:
            t[i] = "A" if t[i] != "A" else "C"
        seq[pos:pos + 23] = t + list("TGG")
    return "".join(seq)


def test_gpu_matches_cpu_on_planted_genome():
    from geneweaver.gpu_kernel import gpu_search
    genome = encode(planted_genome())
    guide = encode(GUIDE)
    cpu = cpu_search(genome, guide, 3)
    gpu = gpu_search(genome, guide, 3)
    assert gpu == cpu
    assert len(cpu) >= 3


def test_gpu_returns_empty_for_tiny_input():
    from geneweaver.gpu_kernel import gpu_search
    assert gpu_search(encode("ACGT"), encode(GUIDE), 3) == []
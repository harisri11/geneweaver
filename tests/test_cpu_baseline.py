import random

import numpy as np

from geneweaver.cpu_baseline import cpu_search
from geneweaver.parser import encode

GUIDE = "GACGCATAAAGATGAGACGC"


def make_genome(n=20_000, seed=1):
    random.seed(seed)
    return list("".join(random.choices("ACGT", k=n)))


def plant(genome, pos, n_mismatches):
    """Write guide (with n_mismatches changed bases) + TGG PAM at pos."""
    target = list(GUIDE)
    for k in range(n_mismatches):
        target[k * 4] = "A" if target[k * 4] != "A" else "C"
    genome[pos:pos + 23] = target + list("TGG")


def search(genome, max_mm=3):
    return cpu_search(encode("".join(genome)), encode(GUIDE), max_mm)


def test_finds_exact_and_mismatched_targets():
    g = make_genome()
    plant(g, 1000, 0)
    plant(g, 5000, 1)
    plant(g, 9000, 2)
    hits = dict(search(g))
    assert hits[1000] == 0
    assert hits[5000] == 1
    assert hits[9000] == 2


def test_rejects_too_many_mismatches():
    g = make_genome()
    plant(g, 3000, 4)  # 4 mismatches > max_mm=3
    assert 3000 not in dict(search(g, max_mm=3))


def test_requires_ngg_pam():
    g = make_genome()
    plant(g, 4000, 0)
    g[4000 + 21] = "A"  # break the PAM: TGG -> TAG
    assert 4000 not in dict(search(g))


def test_no_hits_in_clean_genome():
    assert search(make_genome()) == []
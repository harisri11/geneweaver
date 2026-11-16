# tests/test_pipeline.py
import random

import pytest

from geneweaver.pipeline import run_search

GUIDE = "GACGCATAAAGATGAGACGC"


def make_fasta(tmp_path, pos, n_mm=1, n=200_000):
    random.seed(1)
    seq = list("".join(random.choices("ACGT", k=n)))
    target = list(GUIDE)
    for k in range(n_mm):
        i = 5 + k * 4  # weighted guide positions (index 0 has weight 0)
        target[i] = "A" if target[i] != "A" else "C"
    seq[pos:pos + 23] = target + list("TGG")
    p = tmp_path / "t.fasta"
    p.write_text(">t\n" + "".join(seq) + "\n")
    return str(p)


def test_finds_planted_hit_with_score(tmp_path):
    f = make_fasta(tmp_path, pos=1000, n_mm=1)
    hits = run_search(f, GUIDE, 3, "cpu", chunk_size=50_000)
    hit = next(h for h in hits if h.pos == 1000)
    assert hit.mismatches == 1
    assert hit.chrom == "t"
    assert 0 < hit.score < 100


def test_hit_on_chunk_boundary_found_exactly_once(tmp_path):
    f = make_fasta(tmp_path, pos=49_990)  # straddles the 50,000 boundary
    hits = run_search(f, GUIDE, 3, "cpu", chunk_size=50_000)
    assert sum(h.pos == 49_990 for h in hits) == 1


def test_results_sorted_by_score(tmp_path):
    f = make_fasta(tmp_path, pos=1000, n_mm=0)
    hits = run_search(f, GUIDE, 3, "cpu", chunk_size=50_000)
    scores = [h.score for h in hits]
    assert scores == sorted(scores, reverse=True)
    assert hits[0].score == 100.0


def test_progress_callback_called_per_chunk(tmp_path):
    f = make_fasta(tmp_path, pos=1000)
    events = []
    run_search(f, GUIDE, 3, "cpu", chunk_size=50_000, on_progress=events.append)
    assert [e["done"] for e in events] == [1, 2, 3, 4]


def test_unknown_backend_raises(tmp_path):
    f = make_fasta(tmp_path, pos=1000)
    with pytest.raises(ValueError):
        run_search(f, GUIDE, 3, "tpu")
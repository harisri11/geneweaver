import numpy as np

from geneweaver.parser import encode, iter_chunks


def test_encode_basic():
    assert encode("ACGTN").tolist() == [0, 1, 2, 3, 4]


def test_encode_lowercase():
    assert encode("acgt").tolist() == [0, 1, 2, 3]


def test_chunks_cover_genome_with_overlap(tmp_path):
    seq = "ACGT" * 250  # 1000 bases
    f = tmp_path / "t.fasta"
    f.write_text(">t\n" + seq + "\n")
    chunks = list(iter_chunks(str(f), chunk_size=400, overlap=22))
    assert [s for _, s, _ in chunks] == [0, 400, 800]
    assert [len(a) for _, _, a in chunks] == [422, 422, 200]
    # overlap region is identical in consecutive chunks
    assert np.array_equal(chunks[0][2][400:422], chunks[1][2][:22])

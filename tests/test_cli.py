import random

import pytest

from geneweaver.__main__ import main, parse_args

GUIDE = "GACGCATAAAGATGAGACGC"


def planted_fasta(tmp_path):
    random.seed(1)
    seq = list("".join(random.choices("ACGT", k=50_000)))
    seq[2000:2023] = list(GUIDE + "TGG")
    p = tmp_path / "g.fasta"
    p.write_text(">chrT\n" + "".join(seq) + "\n")
    return str(p)


def test_cli_prints_planted_hit(tmp_path, capsys):
    hits = main(["--fasta", planted_fasta(tmp_path), "--guide", GUIDE])
    out = capsys.readouterr().out
    assert "chrT" in out and "2000" in out
    assert hits[0].pos == 2000 and hits[0].score == 100.0


def test_cli_lowercase_guide_is_accepted():
    assert parse_args(["--fasta", "x", "--guide", GUIDE.lower()]).guide == GUIDE


@pytest.mark.parametrize("bad", ["ACGT", GUIDE + "A", "GACGCATAAAGATGAGACGN"])
def test_cli_rejects_bad_guide(bad):
    with pytest.raises(SystemExit):
        parse_args(["--fasta", "x", "--guide", bad])